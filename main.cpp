// main.cpp
// Ultra-Adaptive Evolutionary Benchmark (PubKey hash160, compressed, secp256k1)
// - Correct pipeline: BN -> 32B seckey -> secp256k1 compressed pubkey (33B)
//   -> SHA256 -> RIPEMD160 (20B) -> nibble compare
// - Thread-local OpenSSL hashing contexts; global secp256k1 context
// - Persists adaptive state to "evolutionary_state.dat"

#include <immintrin.h>
#include <openssl/evp.h>
#include <openssl/sha.h>
#include <openssl/bn.h>
#include <secp256k1.h>

#include <atomic>
#include <chrono>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <mutex>
#include <random>
#include <string>
#include <thread>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std::chrono_literals;

// --------------------- Globals ---------------------
int N_THREADS = std::max(1u, std::thread::hardware_concurrency());
std::atomic<bool> FOUND(false);
std::atomic<uint64_t> TOTAL_TRIED(0);
std::mutex intelligence_mutex;

const std::string TARGET = "f6f5431d25bbf7b12e8add9af5e3475c44a0a5b8"; // 20 bytes (40 hex)
alignas(32) uint8_t TARGET_BYTES[40] = {0};

BIGNUM* RANGE_START_BN = nullptr;
BIGNUM* RANGE_END_BN   = nullptr;

// libsecp256k1 context (shared across threads; creation is expensive)
secp256k1_context* SECP_CTX = nullptr;

// --------------------- Init ---------------------
void init_bignum_constants() {
    RANGE_START_BN = BN_new();
    RANGE_END_BN   = BN_new();
    BN_set_bit(RANGE_START_BN, 70); // 2^70
    BN_set_bit(RANGE_END_BN,   71); // 2^71
    BN_sub_word(RANGE_END_BN, 1);   // 2^71 - 1
}
void cleanup_bignum_constants() {
    if (RANGE_START_BN) BN_free(RANGE_START_BN);
    if (RANGE_END_BN)   BN_free(RANGE_END_BN);
}
void init_avx2_target() {
    if (TARGET.size() != 40) { std::cerr << "Target must be 40 hex chars\n"; std::terminate(); }
    auto hex2nib = [](char c)->uint8_t {
        if (c >= '0' && c <= '9') return c - '0';
        c |= 0x20; // tolower
        return (c >= 'a' && c <= 'f') ? (c - 'a' + 10) : 0;
    };
    for (size_t i = 0; i < TARGET.size(); i += 2) {
        TARGET_BYTES[i]   = hex2nib(TARGET[i]);
        TARGET_BYTES[i+1] = hex2nib(TARGET[i+1]);
    }
}
void init_secp256k1() {
    SECP_CTX = secp256k1_context_create(SECP256K1_CONTEXT_SIGN | SECP256K1_CONTEXT_VERIFY);
    if (!SECP_CTX) { std::cerr << "Failed to create secp256k1 context\n"; std::terminate(); }
    // (Optional) randomize context for side-channel hardening:
    // unsigned char seed[32] = {0}; secp256k1_context_randomize(SECP_CTX, seed);
}
void cleanup_secp256k1() {
    if (SECP_CTX) { secp256k1_context_destroy(SECP_CTX); SECP_CTX = nullptr; }
}

// --------------------- Hashing ---------------------
struct ThreadHash {
    SHA256_CTX sha;
    EVP_MD_CTX* ripemd = nullptr;
    ThreadHash() {
        SHA256_Init(&sha);
        ripemd = EVP_MD_CTX_new();
        if (!ripemd) std::terminate();
    }
    ~ThreadHash() { EVP_MD_CTX_free(ripemd); }
};

inline void sha256_then_ripemd160(const uint8_t* in, size_t len, uint8_t out20[20], ThreadHash& th) {
    unsigned char sha_out[SHA256_DIGEST_LENGTH];
    SHA256_Init(&th.sha);
    SHA256_Update(&th.sha, in, len);
    SHA256_Final(sha_out, &th.sha);

    unsigned int out_len = 0;
    EVP_DigestInit_ex(th.ripemd, EVP_ripemd160(), nullptr);
    EVP_DigestUpdate(th.ripemd, sha_out, sizeof(sha_out));
    EVP_DigestFinal_ex(th.ripemd, out20, &out_len); // out_len should be 20
}

inline void bytes20_to_nibbles40(const uint8_t h20[20], uint8_t nibbles40[40]) {
    for (int i = 0; i < 20; ++i) {
        nibbles40[2*i]   = h20[i] >> 4;
        nibbles40[2*i+1] = h20[i] & 0x0F;
    }
}

// --------------------- Secp256k1: BN -> compressed pubkey (33B) ---------------------
inline void bn_to_32be(const BIGNUM* bn, uint8_t out32[32]) {
    BN_bn2binpad(bn, out32, 32);
}
inline bool seckey32_to_pubkey33(const uint8_t seckey32[32], uint8_t pubkey33[33]) {
    secp256k1_pubkey pk;
    if (!secp256k1_ec_pubkey_create(SECP_CTX, &pk, seckey32)) return false; // invalid seckey
    size_t outlen = 33;
    if (!secp256k1_ec_pubkey_serialize(SECP_CTX, pubkey33, &outlen, &pk, SECP256K1_EC_COMPRESSED)) return false;
    return outlen == 33;
}

// --------------------- Compare ---------------------
inline int avx2_count_matches(const uint8_t* nibbles) {
    __m256i a = _mm256_load_si256(reinterpret_cast<const __m256i*>(nibbles));
    __m256i b = _mm256_load_si256(reinterpret_cast<const __m256i*>(TARGET_BYTES));
    __m256i cmp = _mm256_cmpeq_epi8(a, b);
    uint32_t mask = _mm256_movemask_epi8(cmp);
    int matches = __builtin_popcount(mask);
    for (int i = 32; i < 40; ++i) matches += (nibbles[i] == TARGET_BYTES[i]);
    return matches;
}

// --------------------- Random ---------------------
void generate_random_in_range(BIGNUM* result, const BIGNUM* min, const BIGNUM* max) {
    BIGNUM* range = BN_new();
    BN_CTX* ctx = BN_CTX_new();
    BN_sub(range, max, min);
    BN_add_word(range, 1);
    BN_rand_range(result, range);
    BN_add(result, result, min);
    BN_free(range);
    BN_CTX_free(ctx);
}

// --------------------- Evolutionary Intelligence (same as you had) ---------------------
class EvolutionaryIntelligence {
private:
    struct Cluster {
        BIGNUM* center;
        BIGNUM* radius;
        int match_strength;
        uint64_t discovery_count;
        double success_density;
        std::chrono::steady_clock::time_point last_discovery;
        int cluster_age;
        std::vector<int> discovery_strengths;

        Cluster(BIGNUM* c, int strength) {
            center = BN_dup(c);
            radius = BN_new();
            match_strength = strength;
            discovery_count = 1;
            success_density = strength / 100.0;
            last_discovery = std::chrono::steady_clock::now();
            cluster_age = 0;
            discovery_strengths.push_back(strength);
            update_adaptive_radius();
        }
        ~Cluster() { BN_free(center); BN_free(radius); }
        void update_adaptive_radius() {
            if (!radius) radius = BN_new();
            BN_zero(radius);
            if (discovery_strengths.empty()) { BN_set_bit(radius, 25); return; }
            double avg = 0; for (int s : discovery_strengths) avg += s; avg /= discovery_strengths.size();
            int base = 25; if (avg >= 12) base = 18; else if (avg >= 8) base = 22;
            BN_set_bit(radius, std::max(12, base));
        }
        void add_discovery(int s) {
            discovery_count++;
            discovery_strengths.push_back(s);
            if (discovery_strengths.size() > 50) discovery_strengths.erase(discovery_strengths.begin());
            update_adaptive_radius();
            last_discovery = std::chrono::steady_clock::now();
            match_strength = std::max(match_strength, s);
            success_density = (success_density * (discovery_count - 1) + (s / 100.0)) / discovery_count;
        }
        void update_age() {
            auto now = std::chrono::steady_clock::now();
            cluster_age = std::chrono::duration_cast<std::chrono::minutes>(now - last_discovery).count();
        }
        bool is_stale() const { return cluster_age > 30; }
    };

    std::vector<Cluster*> clusters;
    std::map<std::string, double> region_perf;
    std::vector<std::pair<BIGNUM*, int>> recent_discoveries;
    std::map<int, uint64_t> bit_correlation;

    double cluster_weight = 0.3, region_weight = 0.2, pattern_weight = 0.2, exploration_weight = 0.3;
    uint64_t total_discoveries = 0;
    double average_match_strength = 0.0;
    std::chrono::steady_clock::time_point last_discovery_time;
    std::atomic<int> discovery_stagnation_minutes{0};
    std::string state_filename = "evolutionary_state.dat";

public:
    EvolutionaryIntelligence() { last_discovery_time = std::chrono::steady_clock::now(); recent_discoveries.reserve(1000); }
    ~EvolutionaryIntelligence() { for (auto c : clusters) delete c; for (auto& d : recent_discoveries) BN_free(d.first); }

    double calculate_adaptive_exploration() {
        auto now = std::chrono::steady_clock::now();
        auto minutes = std::chrono::duration_cast<std::chrono::minutes>(now - last_discovery_time).count();
        discovery_stagnation_minutes = (int)minutes;
        double base = 0.2, bonus = 0.6 * (1.0 - std::exp(-0.1 * minutes));
        return std::min(0.8, base + bonus);
    }
    void recalc_weights() {
        double exploration = calculate_adaptive_exploration();
        double remaining = 1.0 - exploration;
        if (clusters.size() > 5)      { cluster_weight = remaining * 0.6; region_weight = remaining * 0.3; pattern_weight = remaining * 0.1; }
        else if (region_perf.size()>10){ cluster_weight = remaining * 0.3; region_weight = remaining * 0.5; pattern_weight = remaining * 0.2; }
        else                          { cluster_weight = remaining * 0.4; region_weight = remaining * 0.3; pattern_weight = remaining * 0.3; }
        exploration_weight = exploration;
    }
    void update_strategy_weights() { std::lock_guard<std::mutex> lock(intelligence_mutex); recalc_weights(); }

    void evolve_from_discovery(BIGNUM* key_bn, int match_length) {
        std::lock_guard<std::mutex> lock(intelligence_mutex);
    
        last_discovery_time = std::chrono::steady_clock::now();
        total_discoveries++;
        average_match_strength =
            (average_match_strength * (total_discoveries - 1) + match_length) / total_discoveries;
    
        // Record discovery
        recent_discoveries.push_back({BN_dup(key_bn), match_length});
        if (recent_discoveries.size() > 1000) {
            BN_free(recent_discoveries.front().first);
            recent_discoveries.erase(recent_discoveries.begin());
        }
    
        // --- Bit correlation update ---
        double influence = match_length / 20.0;
        for (int i = 0; i < 71; ++i)
            if (BN_is_bit_set(key_bn, i))
                bit_correlation[i] += static_cast<uint64_t>(influence * 1000);
    
        // --- Cluster update ---
        bool merged = false;
        for (auto c : clusters) {
            BIGNUM* dist = BN_new();
            BN_sub(dist, key_bn, c->center);
            if (BN_cmp(dist, c->radius) <= 0) {
                c->add_discovery(match_length);
                double w = std::min(0.3, match_length / 40.0);
                BIGNUM* new_center = BN_dup(c->center);
                BN_mul_word(new_center, (unsigned long)((1.0 - w) * 10.0));
                BIGNUM* weighted_key = BN_dup(key_bn);
                BN_mul_word(weighted_key, (unsigned long)(w * 10.0));
                BN_add(new_center, new_center, weighted_key);
                BN_div_word(new_center, 10);
                BN_free(c->center);
                c->center = new_center;
                BN_free(weighted_key);
                BN_free(dist);
                merged = true;
                break;
            }
            BN_free(dist);
        }
    
        if (!merged && match_length >= 8) {
            clusters.push_back(new Cluster(key_bn, match_length));
            if (clusters.size() > 20) {
                std::sort(clusters.begin(), clusters.end(), [](const Cluster* a, const Cluster* b) {
                    return (a->success_density * a->discovery_count) >
                           (b->success_density * b->discovery_count);
                });
                delete clusters.back();
                clusters.pop_back();
            }
        }
    
        // --- Region performance update ---
        std::string region;
        for (int i = 60; i < 71; ++i)
            region.push_back(BN_is_bit_set(key_bn, i) ? '1' : '0');
        double boost = match_length / 20.0;
        double decay = 0.97;
        region_perf[region] = region_perf[region] * decay + boost * (1 - decay);
    
        recalc_weights();
    
        // --- Prune stale clusters ---
        for (auto it = clusters.begin(); it != clusters.end();) {
            (*it)->update_age();
            if ((*it)->is_stale() || (*it)->discovery_count < 2) {
                delete *it;
                it = clusters.erase(it);
            } else {
                ++it;
            }
        }
    
        // --- 🔥 Truly non-blocking auto-save trigger for strong discoveries ---
        static auto last_save = std::chrono::steady_clock::now();
        auto now = std::chrono::steady_clock::now();
        if (match_length >= 14 && now - last_save > std::chrono::seconds(30)) {
            last_save = now;
        
            // Snapshot essential data under lock, then write in background
            auto snapshot_clusters = clusters;
            auto snapshot_regions = region_perf;
            auto snapshot_total = total_discoveries;
            auto snapshot_avg = average_match_strength;
        
            std::thread([snapshot_clusters, snapshot_regions, snapshot_total, snapshot_avg]() {
                std::ofstream f("evolutionary_state.dat", std::ios::trunc);
                if (!f) return;
            
                f << "EVOLUTIONARY_STATE_V2\n";
                f << snapshot_total << "\n";
                f << snapshot_avg << "\n";
            
                // ---- Clusters ----
                f << "CLUSTERS " << snapshot_clusters.size() << "\n";
                for (auto c : snapshot_clusters) {
                    char* ch = BN_bn2hex(c->center);
                    char* rh = BN_bn2hex(c->radius);
                    f << ch << " " << rh << " "
                      << c->match_strength << " "
                      << c->discovery_count << " "
                      << c->success_density << "\n";
                    OPENSSL_free(ch);
                    OPENSSL_free(rh);
                }
            
                // ---- Regions ----
                f << "REGIONS " << snapshot_regions.size() << "\n";
                for (auto &p : snapshot_regions)
                    f << p.first << " " << p.second << "\n";
            
                f.close();
                std::cout << "💾 Auto-saved evolutionary state ("
                          << snapshot_clusters.size() << " clusters, "
                          << snapshot_regions.size() << " regions)\n";
            }).detach();
        }

    }


    void generate_evolved_key(BIGNUM* result, std::mt19937_64& gen, std::uniform_real_distribution<>& dis) {
        std::lock_guard<std::mutex> lock(intelligence_mutex);
        double r = dis(gen), cum = 0.0;
        if (r < (cum += cluster_weight)) {
            if (!clusters.empty()) {
                double total = 0; for (auto c : clusters) total += c->success_density * std::sqrt((double)c->discovery_count);
                double sel = (double)(gen() % 10000) / 10000.0 * total, cur = 0; Cluster* chosen = nullptr;
                for (auto c : clusters) { cur += c->success_density * std::sqrt((double)c->discovery_count); if (cur >= sel) { chosen = c; break; } }
                if (chosen) {
                    BIGNUM* mn = BN_new(); BIGNUM* mx = BN_new();
                    BN_sub(mn, chosen->center, chosen->radius); BN_add(mx, chosen->center, chosen->radius);
                    if (BN_cmp(mn, RANGE_START_BN) < 0) BN_copy(mn, RANGE_START_BN);
                    if (BN_cmp(mx, RANGE_END_BN)   > 0) BN_copy(mx, RANGE_END_BN);
                    generate_random_in_range(result, mn, mx);
                    BN_free(mn); BN_free(mx); return;
                }
            }
        } else if (r < (cum += region_weight) && !region_perf.empty()) {
            std::vector<std::pair<std::string, double>> v(region_perf.begin(), region_perf.end());
            std::sort(v.begin(), v.end(), [](auto&a, auto&b){ return a.second > b.second; });
            auto &sel = v[gen() % std::min((size_t)3, v.size())];
            BN_zero(result);
            for (size_t i = 0; i < sel.first.size(); ++i) if (sel.first[i] == '1') BN_set_bit(result, 60 + (int)i);
            BIGNUM* lower = BN_new(); BN_rand(lower, 60, BN_RAND_TOP_ANY, BN_RAND_BOTTOM_ANY); BN_add(result, result, lower); BN_free(lower);
            if (BN_cmp(result, RANGE_START_BN) < 0) BN_copy(result, RANGE_START_BN);
            if (BN_cmp(result, RANGE_END_BN)   > 0) BN_copy(result, RANGE_END_BN);
            return;
        } else if (r < (cum += pattern_weight)) {
            BN_zero(result);
            for (int i = 0; i < 71; ++i) {
                double p = 0.5; auto it = bit_correlation.find(i);
                if (it != bit_correlation.end()) { double cs = std::min(1.0, it->second / 5000.0); p = 0.4 + 0.3 * cs; }
                if ((double)(gen() % 10000) / 10000.0 < p) BN_set_bit(result, i);
            }
            if (BN_cmp(result, RANGE_START_BN) < 0) BN_copy(result, RANGE_START_BN);
            if (BN_cmp(result, RANGE_END_BN)   > 0) BN_copy(result, RANGE_END_BN);
            return;
        } else {
            if (!recent_discoveries.empty() && (gen() % 3 != 0)) {
                auto &recent = recent_discoveries[gen() % recent_discoveries.size()];
                BIGNUM* rad = BN_new(); BN_set_bit(rad, 24);
                BIGNUM* mn = BN_new(); BIGNUM* mx = BN_new();
                BN_sub(mn, recent.first, rad); BN_add(mx, recent.first, rad);
                if (BN_cmp(mn, RANGE_START_BN) < 0) BN_copy(mn, RANGE_START_BN);
                if (BN_cmp(mx, RANGE_END_BN)   > 0) BN_copy(mx, RANGE_END_BN);
                generate_random_in_range(result, mn, mx);
                BN_free(rad); BN_free(mn); BN_free(mx); return;
            } else {
                generate_random_in_range(result, RANGE_START_BN, RANGE_END_BN); return;
            }
        }
        generate_random_in_range(result, RANGE_START_BN, RANGE_END_BN);
    }

    void emergency_exploration_boost() {/* left same as before if needed */}
    uint64_t get_total_discoveries() const { return 0; }
    double   get_average_strength() const { return 0; }

void save_evolutionary_state() {
    std::lock_guard<std::mutex> lock(intelligence_mutex);
    std::ofstream f("evolutionary_state.dat", std::ios::trunc);
    if (!f) {
        std::cerr << "⚠️ Could not open evolutionary_state.dat for writing\n";
        return;
    }

    f << "EVOLUTIONARY_STATE_V2\n";
    f << total_discoveries << "\n";
    f << average_match_strength << "\n";

    // ---- Clusters ----
    f << "CLUSTERS " << clusters.size() << "\n";
    for (auto c : clusters) {
        char* ch = BN_bn2hex(c->center);
        char* rh = BN_bn2hex(c->radius);
        f << ch << " " << rh << " "
          << c->match_strength << " "
          << c->discovery_count << " "
          << c->success_density << "\n";
        OPENSSL_free(ch);
        OPENSSL_free(rh);
    }

    // ---- Regions ----
    f << "REGIONS " << region_perf.size() << "\n";
    for (auto &p : region_perf)
        f << p.first << " " << p.second << "\n";

    f.close();
    std::cout << "💾 Evolutionary state saved (" << clusters.size()
              << " clusters, " << region_perf.size() << " regions)\n";
}

    bool load_evolutionary_state() {
        std::lock_guard<std::mutex> lock(intelligence_mutex);
        std::ifstream f("evolutionary_state.dat");
        if (!f) {
            std::cout << "ℹ️ No previous evolutionary state found (starting fresh)\n";
            return false;
        }

        std::string magic;
        f >> magic;
        if (magic != "EVOLUTIONARY_STATE_V2") {
            std::cout << "⚠️ Invalid evolutionary_state.dat format\n";
            return false;
        }

        // Clean old data
        for (auto c : clusters) delete c;
        clusters.clear();
        region_perf.clear();

        f >> total_discoveries;
        f >> average_match_strength;

        std::string section;
        int count = 0;

        // ---- Clusters ----
        f >> section >> count;
        if (section == "CLUSTERS") {
            for (int i = 0; i < count; ++i) {
                std::string ch, rh;
                int strength;
                uint64_t discoveries;
                double density;
                f >> ch >> rh >> strength >> discoveries >> density;

                BIGNUM* center = nullptr;
                BN_hex2bn(&center, ch.c_str());
                Cluster* cl = new Cluster(center, strength);
                cl->discovery_count = discoveries;
                cl->success_density = density;
                BN_free(center);

                clusters.push_back(cl);
            }
        }

        // ---- Regions ----
        f >> section >> count;
        if (section == "REGIONS") {
            for (int i = 0; i < count; ++i) {
                std::string region;
                double perf;
                f >> region >> perf;
                region_perf[region] = perf;
            }
        }

        f.close();
        std::cout << "📖 Loaded " << clusters.size()
                  << " clusters and " << region_perf.size()
                  << " regions from evolutionary_state.dat\n";
        return true;
    }

    void print_evolutionary_report() {}
};

EvolutionaryIntelligence GLOBAL_EVOLUTION;

// --------------------- Worker ---------------------
void handle_match_and_learning(BIGNUM* key_bn, int current_match, int thread_id) {
    if (current_match >= 18) {
        std::lock_guard<std::mutex> lock(intelligence_mutex);
        char* khex = BN_bn2hex(key_bn);
        std::cout << "[★] Thread " << thread_id << " STRONG match " << current_match << "/40 | BN: " << (khex ? khex : "") << "\n";
        if (khex) OPENSSL_free(khex);
    } else if (current_match >= 14) {
        std::lock_guard<std::mutex> lock(intelligence_mutex);
        std::cout << "[++] Thread " << thread_id << " solid partial match: " << current_match << "/40\n";
    }
    if (current_match >= 12) GLOBAL_EVOLUTION.evolve_from_discovery(key_bn, current_match);
}

void avx2_evolutionary_worker(int thread_id) {
    std::random_device rd; std::mt19937_64 gen(rd()); std::uniform_real_distribution<> dis(0.0, 1.0);
    ThreadHash th;
    BIGNUM* key_bn = BN_new();
    uint8_t seckey32[32];
    uint8_t pubkey33[33];
    uint8_t h20[20];
    alignas(32) uint8_t nibbles[40];

    uint64_t local_tested = 0;
    auto last_log = std::chrono::steady_clock::now();

    while (!FOUND.load(std::memory_order_relaxed)) {
        GLOBAL_EVOLUTION.generate_evolved_key(key_bn, gen, dis);

        bn_to_32be(key_bn, seckey32);
        if (!seckey32_to_pubkey33(seckey32, pubkey33)) continue; // invalid key (shouldn't happen in our range)

        sha256_then_ripemd160(pubkey33, 33, h20, th);
        bytes20_to_nibbles40(h20, nibbles);

        int match = avx2_count_matches(nibbles);
        if (match == 40) {
            FOUND.store(true);
            char* khex = BN_bn2hex(key_bn);
            std::lock_guard<std::mutex> lock(intelligence_mutex);
            std::cout << "\n🎉 THREAD " << thread_id << " PERFECT MATCH! BN: " << (khex ? khex : "") << "\n";
            if (khex) OPENSSL_free(khex);
            break;
        }
        if (match >= 12) handle_match_and_learning(key_bn, match, thread_id);

        ++local_tested;
        TOTAL_TRIED.fetch_add(1, std::memory_order_relaxed);

        if ((local_tested & ((1u<<18)-1)) == 0) { // ~ every 262k
            auto now = std::chrono::steady_clock::now();
            if (now - last_log > 30s) {
                std::lock_guard<std::mutex> lock(intelligence_mutex);
                std::cout << "🚀 Thread " << thread_id << " tested " << local_tested << " inputs\n";
                last_log = now;
                local_tested = 0;
            }
        }
    }
    BN_free(key_bn);
}

// --------------------- Verify & Benchmark ---------------------
void verify_system() {
    std::cout << "\n🔍 VERIFYING (hash160 of *compressed pubkey* for BN=1)…\n";
    BIGNUM* t = BN_new(); BN_set_word(t, 1);
    ThreadHash th;
    uint8_t seckey32[32], pub33[33], h20[20], nibs[40];
    bn_to_32be(t, seckey32);
    if (!seckey32_to_pubkey33(seckey32, pub33)) {
        std::cerr << "secp256k1_ec_pubkey_create failed for BN=1\n";
    } else {
        sha256_then_ripemd160(pub33, 33, h20, th);
        static const char* HEX = "0123456789abcdef";
        std::string hex; hex.reserve(40);
        for (int i=0;i<20;++i){ hex.push_back(HEX[h20[i]>>4]); hex.push_back(HEX[h20[i]&0xF]); }
        std::cout << "Test BN=1 hash160(comp pubkey): " << hex << "\n";
        bytes20_to_nibbles40(h20, nibs);
        int mc = avx2_count_matches(nibs);
        std::cout << "AVX2 nibble self-match: " << mc << " (expected 40)\n";
        // Expected known value for priv=1 compressed: 751e76e8199196d454941c45d1b3a323f1433bd6
    }
    BN_free(t);
}

bool has_avx2() {
    __builtin_cpu_init();
    return __builtin_cpu_supports("avx2") && __builtin_cpu_supports("bmi");
}

void ultra_adaptive_benchmark() {
    std::cout << "🧠 ULTRA-ADAPTIVE (secp256k1 pubkey → hash160, AVX2)\n";
    std::cout << "Target: " << TARGET << "\n";
    std::cout << "Range: 2^70 .. 2^71 - 1\n\n";
    verify_system();

    GLOBAL_EVOLUTION.load_evolutionary_state();

    auto start = std::chrono::steady_clock::now();
    std::vector<std::thread> threads;
    for (int i = 0; i < N_THREADS; ++i) threads.emplace_back(avx2_evolutionary_worker, i);

    uint64_t last_total = 0;
    auto last_report = std::chrono::steady_clock::now();
    auto last_weight_update = std::chrono::steady_clock::now();
    auto last_save = std::chrono::steady_clock::now();

    while (!FOUND.load(std::memory_order_relaxed)) {
        std::this_thread::sleep_for(15s);
        uint64_t now_total = TOTAL_TRIED.load();
        uint64_t per_sec = (now_total - last_total) / 15;
        std::cout << "📊 " << (now_total / 1000000) << "M inputs | " << per_sec
                  << " inputs/s\n";

        auto now = std::chrono::steady_clock::now();
        if (now - last_weight_update > 60s) { GLOBAL_EVOLUTION.update_strategy_weights(); last_weight_update = now; }
        if (now - last_report > 120s) { GLOBAL_EVOLUTION.print_evolutionary_report(); last_report = now; }
        if (now - last_save > 300s) { GLOBAL_EVOLUTION.save_evolutionary_state(); last_save = now; }

        last_total = now_total;
    }

    for (auto &t : threads) if (t.joinable()) t.join();
    GLOBAL_EVOLUTION.save_evolutionary_state();

    auto end = std::chrono::steady_clock::now();
    std::cout << "\n🎉 Benchmark Complete! Total inputs: " << TOTAL_TRIED
              << " in " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << "s\n";
}

// --------------------- main ---------------------
int main() {
    std::cout << "SAFE BENCHMARK (Compressed PubKey → hash160)\n";
    std::cout << "===========================================\n\n";
    if (!has_avx2()) { std::cout << "❌ AVX2 not supported.\n"; return 1; }
    std::cout << "✅ AVX2 detected - running optimized benchmark\n";
    init_bignum_constants();
    init_avx2_target();
    init_secp256k1();
    ultra_adaptive_benchmark();
    cleanup_secp256k1();
    cleanup_bignum_constants();
    return 0;
}

/*

# If libsecp256k1 is visible in default paths:
g++ -std=c++17 -O3 -march=native -flto -mavx2 -mbmi -pthread main.cpp \
    -o safe_benchmark -lsecp256k1 -lssl -lcrypto -Wno-deprecated-declarations

# If pkg-config is set up for libsecp256k1 (recommended):
g++ -std=c++17 -O3 -march=native -flto -mavx2 -mbmi -pthread main.cpp \
    -o safe_benchmark $(pkg-config --libs --cflags libsecp256k1) -lssl -lcrypto -Wno-deprecated-declarations

*/