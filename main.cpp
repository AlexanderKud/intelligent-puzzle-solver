#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <chrono>
#include <cmath>
#include <thread>
#include <atomic>
#include <mutex>
#include <map>
#include <algorithm>
#include <queue>
#include <bitset>
#include <fstream>
#include <openssl/sha.h>
#include <openssl/ripemd.h>
#include <openssl/ec.h>
#include <openssl/obj_mac.h>
#include <openssl/bn.h>

// 🧠 ULTRA-ADAPTIVE EVOLUTIONARY BITCOIN PUZZLE SOLVER
int N_THREADS = std::thread::hardware_concurrency();
constexpr int ANALYSIS_WINDOW = 10000;

std::atomic<bool> FOUND(false);
std::mutex intelligence_mutex;
std::atomic<uint64_t> TOTAL_TRIED(0);

const std::string TARGET = "f6f5431d25bbf7b12e8add9af5e3475c44a0a5b8";

BIGNUM* RANGE_START_BN = nullptr;
BIGNUM* RANGE_END_BN = nullptr;

void init_bignum_constants() {
    RANGE_START_BN = BN_new();
    RANGE_END_BN = BN_new();
    
    BN_set_bit(RANGE_START_BN, 70);
    BN_set_bit(RANGE_END_BN, 71);
    BN_sub_word(RANGE_END_BN, 1);
}

void cleanup_bignum_constants() {
    if (RANGE_START_BN) BN_free(RANGE_START_BN);
    if (RANGE_END_BN) BN_free(RANGE_END_BN);
}

// 🧠 ULTRA-ADAPTIVE EVOLUTIONARY INTELLIGENCE
class EvolutionaryIntelligence {
private:
    struct EvolutionaryCluster {
        BIGNUM* center;
        BIGNUM* radius;
        int match_strength;
        uint64_t discovery_count;
        double success_density;
        std::chrono::steady_clock::time_point last_discovery;
        int cluster_age;
        std::vector<int> discovery_strengths;
        
        EvolutionaryCluster(BIGNUM* center_bn, int strength) {
            center = BN_dup(center_bn);
            radius = BN_new();
            match_strength = strength;
            discovery_count = 1;
            success_density = strength / 100.0;
            last_discovery = std::chrono::steady_clock::now();
            cluster_age = 0;
            discovery_strengths.push_back(strength);
            update_adaptive_radius();
        }
        
        ~EvolutionaryCluster() {
            BN_free(center);
            BN_free(radius);
        }
        
        void update_age() {
            auto now = std::chrono::steady_clock::now();
            cluster_age = std::chrono::duration_cast<std::chrono::minutes>(now - last_discovery).count();
        }
        
        bool is_stale() const {
            return cluster_age > 45;
        }
        
        void update_adaptive_radius() {
            BN_free(radius);
            radius = BN_new();
            
            if (discovery_strengths.empty()) {
                BN_set_bit(radius, 25);
                return;
            }
            
            double avg_strength = 0.0;
            for (int strength : discovery_strengths) {
                avg_strength += strength;
            }
            avg_strength /= discovery_strengths.size();
            
            int base_radius_bits = 25;
            
            if (avg_strength >= 14) {
                base_radius_bits = 18 - std::min(static_cast<unsigned long long>(discovery_count / 5), 5ULL);
            } else if (avg_strength >= 10) {
                base_radius_bits = 22 - std::min(static_cast<unsigned long long>(discovery_count / 10), 3ULL);
            } else {
                base_radius_bits = 25;
            }
            
            BN_set_bit(radius, std::max(10, base_radius_bits));
        }
        
        void add_discovery(int strength) {
            discovery_count++;
            discovery_strengths.push_back(strength);
            if (discovery_strengths.size() > 50) {
                discovery_strengths.erase(discovery_strengths.begin());
            }
            update_adaptive_radius();
            last_discovery = std::chrono::steady_clock::now();
            match_strength = std::max(match_strength, strength);
            success_density = (success_density * (discovery_count - 1) + (strength / 100.0)) / discovery_count;
        }
    };
    
    std::vector<EvolutionaryCluster*> active_clusters;
    std::map<std::string, double> region_performance;
    std::vector<std::pair<BIGNUM*, int>> recent_discoveries;
    std::map<int, uint64_t> bit_correlation;
    
    double cluster_weight = 0.4;
    double region_weight = 0.3;
    double pattern_weight = 0.2;
    double exploration_weight = 0.1;
    
    uint64_t total_discoveries = 0;
    double average_match_strength = 0.0;
    
    std::chrono::steady_clock::time_point last_discovery_time;
    std::atomic<int> discovery_stagnation_minutes{0};
    std::string state_filename = "evolutionary_state.dat";

public:
    // 🚨 EMERGENCY EXPLORATION BOOST
    void emergency_exploration_boost() {
        std::lock_guard<std::mutex> lock(intelligence_mutex);
        
        if (total_discoveries == 0 && discovery_stagnation_minutes > 2) {
            std::cout << "🚨 EMERGENCY: No discoveries after " << discovery_stagnation_minutes 
                      << " minutes and " << TOTAL_TRIED << " keys. Boosting exploration to 80%!\n";
            
            exploration_weight = 0.8;
            cluster_weight = 0.1;
            region_weight = 0.05;
            pattern_weight = 0.05;
            
            // Also add some random clusters to kickstart learning
            if (active_clusters.empty()) {
                std::cout << "🧠 Seeding random exploration clusters...\n";
                std::random_device rd;
                std::mt19937_64 gen(rd());
                
                for (int i = 0; i < 5; i++) {
                    BIGNUM* random_center = BN_new();
                    generate_random_in_range(random_center, RANGE_START_BN, RANGE_END_BN, gen);
                    
                    // Create cluster with simulated 6-char match to start learning
                    EvolutionaryCluster* cluster = new EvolutionaryCluster(random_center, 6);
                    active_clusters.push_back(cluster);
                    
                    BN_free(random_center);
                }
            }
        }
    }
    
    uint64_t get_total_discoveries() const {
        return total_discoveries;
    }

    EvolutionaryIntelligence() {
        last_discovery_time = std::chrono::steady_clock::now();
        recent_discoveries.reserve(1000);
    }
    
    ~EvolutionaryIntelligence() {
        for (auto cluster : active_clusters) {
            delete cluster;
        }
    }
    
    double calculate_adaptive_exploration() {
        auto now = std::chrono::steady_clock::now();
        auto minutes_since_discovery = std::chrono::duration_cast<std::chrono::minutes>(
            now - last_discovery_time).count();
        
        discovery_stagnation_minutes = minutes_since_discovery;
        
        double base_exploration = 0.1;
        double stagnation_bonus = 0.4 * (1.0 - exp(-0.05 * minutes_since_discovery));
        
        return std::min(0.5, base_exploration + stagnation_bonus);
    }
    
    void recalculate_strategy_weights() {
        double exploration = calculate_adaptive_exploration();
        double remaining = 1.0 - exploration;
        
        cluster_weight = remaining * 0.5;
        region_weight = remaining * 0.3;
        pattern_weight = remaining * 0.2;
        
        exploration_weight = exploration;
    }
    
    void evolve_from_discovery(BIGNUM* key_bn, int match_length, const std::string& hash_result) {
        std::lock_guard<std::mutex> lock(intelligence_mutex);
        
        last_discovery_time = std::chrono::steady_clock::now();
        total_discoveries++;
        average_match_strength = (average_match_strength * (total_discoveries - 1) + match_length) / total_discoveries;
        
        recent_discoveries.push_back({BN_dup(key_bn), match_length});
        if (recent_discoveries.size() > 1000) {
            BN_free(recent_discoveries[0].first);
            recent_discoveries.erase(recent_discoveries.begin());
        }
        
        evolve_bit_patterns(key_bn, match_length);
        evolve_clusters(key_bn, match_length);
        evolve_region_performance(key_bn, match_length);
        recalculate_strategy_weights();
        prune_stale_clusters();
    }
    
    void evolve_bit_patterns(BIGNUM* key_bn, int match_length) {
        std::bitset<256> bit_pattern;
        for (int i = 0; i < 256; i++) {
            if (BN_is_bit_set(key_bn, i)) {
                bit_pattern.set(i);
            }
        }
        
        double influence = match_length / 20.0;
        for (int i = 0; i < 71; i++) {
            if (bit_pattern[i]) {
                bit_correlation[i] += static_cast<uint64_t>(influence * 1000);
            }
        }
    }
    
    void evolve_clusters(BIGNUM* key_bn, int match_length) {
        if (match_length < 10) return;
        
        bool merged = false;
        
        for (auto cluster : active_clusters) {
            BIGNUM* distance = BN_new();
            BN_CTX* ctx = BN_CTX_new();
            
            BN_sub(distance, key_bn, cluster->center);
            
            if (BN_cmp(distance, cluster->radius) <= 0) {
                cluster->add_discovery(match_length);
                
                BIGNUM* new_center = BN_new();
                BN_CTX* ctx2 = BN_CTX_new();
                
                double weight = match_length / 20.0;
                BN_mul_word(cluster->center, static_cast<int>((1 - weight) * 10));
                BN_mul_word(key_bn, static_cast<int>(weight * 10));
                BN_add(cluster->center, cluster->center, key_bn);
                BN_div_word(cluster->center, 10);
                
                BN_free(new_center);
                BN_CTX_free(ctx2);
                merged = true;
            }
            
            BN_free(distance);
            BN_CTX_free(ctx);
        }
        
        if (!merged && match_length >= 12) {
            active_clusters.push_back(new EvolutionaryCluster(key_bn, match_length));
            
            if (active_clusters.size() > 25) {
                std::sort(active_clusters.begin(), active_clusters.end(),
                         [](const EvolutionaryCluster* a, const EvolutionaryCluster* b) {
                             return (a->success_density * a->discovery_count) > 
                                    (b->success_density * b->discovery_count);
                         });
                delete active_clusters.back();
                active_clusters.pop_back();
            }
        }
    }
    
    void evolve_region_performance(BIGNUM* key_bn, int match_length) {
        std::string region_id = get_evolutionary_region(key_bn);
        double performance_boost = match_length / 20.0;
        double decay_factor = 0.95;
        
        region_performance[region_id] = 
            region_performance[region_id] * decay_factor + performance_boost * (1 - decay_factor);
    }
    
    void prune_stale_clusters() {
        for (auto it = active_clusters.begin(); it != active_clusters.end(); ) {
            (*it)->update_age();
            if ((*it)->is_stale()) {
                delete *it;
                it = active_clusters.erase(it);
            } else {
                ++it;
            }
        }
    }
    
    std::string get_evolutionary_region(BIGNUM* key_bn) {
        std::string region_id;
        for (int i = 63; i < 71; i++) {
            region_id += BN_is_bit_set(key_bn, i) ? '1' : '0';
        }
        for (int i = 31; i < 39; i++) {
            region_id += BN_is_bit_set(key_bn, i) ? '1' : '0';
        }
        return region_id;
    }
    
    void generate_evolved_key(BIGNUM* result, std::mt19937_64& gen, std::uniform_real_distribution<>& dis) {
        std::lock_guard<std::mutex> lock(intelligence_mutex);
        
        double strategy_selector = dis(gen);
        double cumulative = 0.0;
        
        if (strategy_selector < (cumulative += cluster_weight) && !active_clusters.empty()) {
            generate_from_evolved_clusters(result, gen);
        } else if (strategy_selector < (cumulative += region_weight) && !region_performance.empty()) {
            generate_from_high_performance_regions(result, gen);
        } else if (strategy_selector < (cumulative += pattern_weight)) {
            generate_from_evolved_patterns(result, gen);
        } else {
            generate_evolutionary_exploration(result, gen);
        }
    }
    
    void generate_from_evolved_clusters(BIGNUM* result, std::mt19937_64& gen) {
        double total_fitness = 0;
        for (auto cluster : active_clusters) {
            total_fitness += cluster->success_density * cluster->discovery_count;
        }
        
        double selection = static_cast<double>(gen() % 10000) / 10000.0 * total_fitness;
        double current = 0;
        
        EvolutionaryCluster* selected_cluster = nullptr;
        for (auto cluster : active_clusters) {
            current += cluster->success_density * cluster->discovery_count;
            if (current >= selection) {
                selected_cluster = cluster;
                break;
            }
        }
        
        if (selected_cluster) {
            BIGNUM* min_val = BN_new();
            BIGNUM* max_val = BN_new();
            
            BN_sub(min_val, selected_cluster->center, selected_cluster->radius);
            BN_add(max_val, selected_cluster->center, selected_cluster->radius);
            
            if (BN_cmp(min_val, RANGE_START_BN) < 0) BN_copy(min_val, RANGE_START_BN);
            if (BN_cmp(max_val, RANGE_END_BN) > 0) BN_copy(max_val, RANGE_END_BN);
            
            generate_random_in_range(result, min_val, max_val, gen);
            
            BN_free(min_val);
            BN_free(max_val);
        } else {
            generate_random_in_range(result, RANGE_START_BN, RANGE_END_BN, gen);
        }
    }
    
    void generate_from_high_performance_regions(BIGNUM* result, std::mt19937_64& gen) {
        std::vector<std::pair<std::string, double>> sorted_regions(
            region_performance.begin(), region_performance.end());
        std::sort(sorted_regions.begin(), sorted_regions.end(),
                 [](const auto& a, const auto& b) { return a.second > b.second; });
        
        if (sorted_regions.empty()) return;
        
        const auto& selected_region = sorted_regions[gen() % std::min(5, (int)sorted_regions.size())];
        
        BN_zero(result);
        std::string region_pattern = selected_region.first;
        
        for (size_t i = 0; i < 8; i++) {
            if (region_pattern[i] == '1') {
                BN_set_bit(result, 63 + i);
            }
        }
        
        for (size_t i = 8; i < 16; i++) {
            if (region_pattern[i] == '1') {
                BN_set_bit(result, 31 + (i - 8));
            }
        }
        
        BIGNUM* random_lower = BN_new();
        BN_rand(random_lower, 31, BN_RAND_TOP_ANY, BN_RAND_BOTTOM_ANY);
        BN_add(result, result, random_lower);
        BN_free(random_lower);
        
        if (BN_cmp(result, RANGE_START_BN) < 0) BN_copy(result, RANGE_START_BN);
        if (BN_cmp(result, RANGE_END_BN) > 0) BN_copy(result, RANGE_END_BN);
    }
    
    void generate_from_evolved_patterns(BIGNUM* result, std::mt19937_64& gen) {
        BN_zero(result);
        
        for (int i = 0; i < 71; i++) {
            double bit_probability = 0.5;
            
            if (bit_correlation.find(i) != bit_correlation.end()) {
                double correlation_strength = std::min(1.0, bit_correlation[i] / 10000.0);
                bit_probability = 0.3 + 0.4 * correlation_strength;
            }
            
            if (static_cast<double>(gen() % 10000) / 10000.0 < bit_probability) {
                BN_set_bit(result, i);
            }
        }
        
        if (BN_cmp(result, RANGE_START_BN) < 0) BN_copy(result, RANGE_START_BN);
        if (BN_cmp(result, RANGE_END_BN) > 0) BN_copy(result, RANGE_END_BN);
    }
    
    void generate_evolutionary_exploration(BIGNUM* result, std::mt19937_64& gen) {
        if (!recent_discoveries.empty() && (gen() % 4 != 0)) {
            const auto& recent = recent_discoveries[gen() % recent_discoveries.size()];
            BIGNUM* exploration_center = recent.first;
            
            BIGNUM* radius = BN_new();
            BN_set_bit(radius, 25);
            
            BIGNUM* min_val = BN_new();
            BIGNUM* max_val = BN_new();
            
            BN_sub(min_val, exploration_center, radius);
            BN_add(max_val, exploration_center, radius);
            
            if (BN_cmp(min_val, RANGE_START_BN) < 0) BN_copy(min_val, RANGE_START_BN);
            if (BN_cmp(max_val, RANGE_END_BN) > 0) BN_copy(max_val, RANGE_END_BN);
            
            generate_random_in_range(result, min_val, max_val, gen);
            
            BN_free(radius);
            BN_free(min_val);
            BN_free(max_val);
        } else {
            generate_random_in_range(result, RANGE_START_BN, RANGE_END_BN, gen);
        }
    }
    
    void generate_random_in_range(BIGNUM* result, const BIGNUM* min, const BIGNUM* max, std::mt19937_64& gen) {
        BIGNUM* range = BN_new();
        BN_CTX* ctx = BN_CTX_new();
        
        BN_sub(range, max, min);
        BN_add_word(range, 1);
        
        BN_rand_range(result, range);
        BN_add(result, result, min);
        
        BN_free(range);
        BN_CTX_free(ctx);
    }
    
    void save_evolutionary_state() {
        std::lock_guard<std::mutex> lock(intelligence_mutex);
        
        std::ofstream state_file(state_filename);
        if (!state_file) {
            std::cout << "❌ Could not save state to " << state_filename << "\n";
            return;
        }
        
        state_file << "EVOLUTIONARY_STATE_V1\n";
        state_file << total_discoveries << "\n";
        state_file << average_match_strength << "\n";
        
        state_file << "CLUSTERS " << active_clusters.size() << "\n";
        for (auto cluster : active_clusters) {
            char* center_hex = BN_bn2hex(cluster->center);
            char* radius_hex = BN_bn2hex(cluster->radius);
            state_file << center_hex << " " << radius_hex << " "
                      << cluster->match_strength << " " << cluster->discovery_count << " "
                      << cluster->success_density << " " << cluster->cluster_age << "\n";
            OPENSSL_free(center_hex);
            OPENSSL_free(radius_hex);
        }
        
        state_file << "REGIONS " << region_performance.size() << "\n";
        for (const auto& [region, performance] : region_performance) {
            state_file << region << " " << performance << "\n";
        }
        
        state_file << "BITS " << bit_correlation.size() << "\n";
        for (const auto& [bit, correlation] : bit_correlation) {
            state_file << bit << " " << correlation << "\n";
        }
        
        state_file << "STRATEGY_WEIGHTS " << cluster_weight << " " << region_weight 
                  << " " << pattern_weight << " " << exploration_weight << "\n";
        
        state_file.close();
        std::cout << "💾 Evolutionary state saved to " << state_filename << "\n";
    }
    
    bool load_evolutionary_state() {
        std::lock_guard<std::mutex> lock(intelligence_mutex);
        
        std::ifstream state_file(state_filename);
        if (!state_file) {
            std::cout << "📝 No previous state found, starting fresh evolution\n";
            return false;
        }
        
        std::string magic;
        state_file >> magic;
        if (magic != "EVOLUTIONARY_STATE_V1") {
            std::cout << "❌ Invalid state file format\n";
            return false;
        }
        
        for (auto cluster : active_clusters) {
            delete cluster;
        }
        active_clusters.clear();
        region_performance.clear();
        bit_correlation.clear();
        
        state_file >> total_discoveries >> average_match_strength;
        
        std::string section;
        int count;
        
        state_file >> section >> count;
        for (int i = 0; i < count; i++) {
            std::string center_hex, radius_hex;
            int strength, discoveries, age;
            double density;
            
            state_file >> center_hex >> radius_hex >> strength >> discoveries >> density >> age;
            
            BIGNUM* center = nullptr;
            BN_hex2bn(&center, center_hex.c_str());
            
            EvolutionaryCluster* cluster = new EvolutionaryCluster(center, strength);
            cluster->discovery_count = discoveries;
            cluster->success_density = density;
            cluster->cluster_age = age;
            
            BN_free(center);
            active_clusters.push_back(cluster);
        }
        
        state_file >> section >> count;
        for (int i = 0; i < count; i++) {
            std::string region;
            double performance;
            state_file >> region >> performance;
            region_performance[region] = performance;
        }
        
        state_file >> section >> count;
        for (int i = 0; i < count; i++) {
            int bit;
            uint64_t correlation;
            state_file >> bit >> correlation;
            bit_correlation[bit] = correlation;
        }
        
        state_file >> section >> cluster_weight >> region_weight >> pattern_weight >> exploration_weight;
        
        state_file.close();
        std::cout << "📖 Evolutionary state loaded: " << active_clusters.size() 
                  << " clusters, " << region_performance.size() << " regions\n";
        return true;
    }
    
    void print_evolutionary_report() {
        std::lock_guard<std::mutex> lock(intelligence_mutex);
        
        std::cout << "\n🧠 ULTRA-ADAPTIVE EVOLUTIONARY REPORT:\n";
        std::cout << "Total Discoveries: " << total_discoveries << "\n";
        std::cout << "Average Match Strength: " << average_match_strength << "\n";
        std::cout << "Discovery Stagnation: " << discovery_stagnation_minutes << " minutes\n";
        std::cout << "Active Clusters: " << active_clusters.size() << "\n";
        std::cout << "Tracked Regions: " << region_performance.size() << "\n";
        
        std::cout << "Adaptive Strategy Weights:\n";
        std::cout << "  Cluster: " << (cluster_weight * 100) << "%\n";
        std::cout << "  Region: " << (region_weight * 100) << "%\n";
        std::cout << "  Pattern: " << (pattern_weight * 100) << "%\n";
        std::cout << "  Exploration: " << (exploration_weight * 100) << "%\n";
        
        if (!active_clusters.empty()) {
            std::cout << "Top Evolutionary Clusters:\n";
            std::sort(active_clusters.begin(), active_clusters.end(),
                     [](const auto& a, const auto& b) { 
                         return a->success_density > b->success_density; 
                     });
            for (size_t i = 0; i < std::min(active_clusters.size(), size_t(3)); i++) {
                char* center_hex = BN_bn2hex(active_clusters[i]->center);
                char* radius_hex = BN_bn2hex(active_clusters[i]->radius);
                std::cout << "  Cluster " << i << ": strength=" << active_clusters[i]->match_strength
                          << ", density=" << active_clusters[i]->success_density 
                          << ", discoveries=" << active_clusters[i]->discovery_count 
                          << ", radius_bits=" << BN_num_bits(active_clusters[i]->radius)
                          << ", center=" << std::string(center_hex).substr(0, 16) << "...\n";
                OPENSSL_free(center_hex);
                OPENSSL_free(radius_hex);
            }
        }
    }
};

EvolutionaryIntelligence GLOBAL_EVOLUTION;

// 🧠 TURBO HASHING FUNCTION
std::string turbo_hash_bn(BIGNUM* key_bn) {
    thread_local EC_KEY* kk = nullptr;
    thread_local const EC_GROUP* g = nullptr;
    thread_local EC_POINT* p = nullptr;
    
    if (!kk) {
        kk = EC_KEY_new_by_curve_name(NID_secp256k1);
        g = EC_KEY_get0_group(kk);
        p = EC_POINT_new(g);
    }
    
    EC_KEY_set_private_key(kk, key_bn);
    EC_POINT_mul(g, p, key_bn, NULL, NULL, NULL);
    
    unsigned char pub[65];
    size_t pl = EC_POINT_point2oct(g, p, POINT_CONVERSION_COMPRESSED, pub, 65, NULL);
    
    unsigned char sha[SHA256_DIGEST_LENGTH];
    SHA256(pub, pl, sha);
    
    unsigned char rm[RIPEMD160_DIGEST_LENGTH];
    RIPEMD160(sha, SHA256_DIGEST_LENGTH, rm);
    
    thread_local char out[41];
    static const char* hex = "0123456789abcdef";
    for (int i = 0; i < 20; i++) {
        out[i * 2] = hex[rm[i] >> 4];
        out[i * 2 + 1] = hex[rm[i] & 0xF];
    }
    
    return std::string(out, 40);
}

// 🧠 EVOLUTIONARY WORKER
void evolutionary_worker(int thread_id) {
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    
    EC_KEY* kk = EC_KEY_new_by_curve_name(NID_secp256k1);
    BIGNUM* key_bn = BN_new();
    const EC_GROUP* g = EC_KEY_get0_group(kk);
    EC_POINT* p = EC_POINT_new(g);
    unsigned char pub[65];
    unsigned char sha[SHA256_DIGEST_LENGTH];
    unsigned char rm[RIPEMD160_DIGEST_LENGTH];
    char out[41];
    static const char* hex_chars = "0123456789abcdef";
    
    uint64_t local_tested = 0;
    auto last_report = std::chrono::steady_clock::now();
    
    while (!FOUND) {
        GLOBAL_EVOLUTION.generate_evolved_key(key_bn, gen, dis);
        
        EC_KEY_set_private_key(kk, key_bn);
        EC_POINT_mul(g, p, key_bn, NULL, NULL, NULL);
        size_t pl = EC_POINT_point2oct(g, p, POINT_CONVERSION_COMPRESSED, pub, 65, NULL);
        SHA256(pub, pl, sha);
        RIPEMD160(sha, SHA256_DIGEST_LENGTH, rm);
        
        bool full_match = true;
        int current_match = 0;
        
        for (int i = 0; i < 20; i++) {
            char c1 = hex_chars[rm[i] >> 4];
            char c2 = hex_chars[rm[i] & 0xF];
            
            if (c1 != TARGET[i * 2] || c2 != TARGET[i * 2 + 1]) {
                full_match = false;
                current_match = i * 2;
                if (c1 == TARGET[i * 2]) {
                    current_match++;
                }
                break;
            }
        }
        
        if (full_match) {
            FOUND = true;
            char* found_key_hex = BN_bn2hex(key_bn);
            std::lock_guard<std::mutex> lock(intelligence_mutex);
            std::cout << "\n🎉 EVOLUTIONARY WORKER " << thread_id << " SOLVED THE PUZZLE!\n";
            std::cout << "Key: " << found_key_hex << "\n";
            OPENSSL_free(found_key_hex);
            break;
        }
        
        // 🚨 LOWERED THRESHOLD FROM 8→6 AND 14→10
        if (current_match >= 6) {
            GLOBAL_EVOLUTION.evolve_from_discovery(key_bn, current_match, "");
            
            if (current_match >= 10) {
                char* key_hex = BN_bn2hex(key_bn);
                std::cout << "🔥 Evolutionary breakthrough: Thread " << thread_id 
                          << " found " << current_match << " char match!\n";
                OPENSSL_free(key_hex);
            }
        }
        
        local_tested++;
        TOTAL_TRIED++;
        
        if (local_tested % 50000 == 0) {
            auto now = std::chrono::steady_clock::now();
            if (now - last_report > std::chrono::seconds(30)) {
                std::lock_guard<std::mutex> lock(intelligence_mutex);
                std::cout << "🧬 Evolutionary Thread " << thread_id << " tested " << local_tested 
                          << " keys (total: " << TOTAL_TRIED << ")\n";
                last_report = now;
                local_tested = 0;
            }
        }
    }
    
    EC_POINT_free(p);
    BN_free(key_bn);
    EC_KEY_free(kk);
}

// 🧠 PARALLEL PATTERN DISCOVERY WORKER
void parallel_pattern_discovery_worker(int worker_id, int total_workers) {
    std::cout << "🔍 Continuous Pattern Discovery Worker " << worker_id << " started...\n";
    
    std::random_device rd;
    std::mt19937_64 gen(rd());
    
    EC_KEY* kk = EC_KEY_new_by_curve_name(NID_secp256k1);
    BIGNUM* test_bn = BN_new();
    
    uint64_t base_offset = worker_id * (1ULL << 30);
    uint64_t batch_counter = 0;
    
    while (!FOUND) {
        std::vector<BIGNUM*> evolutionary_patterns;
        batch_counter++;
        
        // 🧠 DYNAMIC PATTERN GENERATION - DIFFERENT STRATEGIES PER WORKER
        if (worker_id == 0) {
            // Powers of 2 with random offsets
            for (int exp = 70; exp < 72; exp++) {
                BIGNUM* base = BN_new();
                BN_set_bit(base, exp);
                
                // Generate 100 random offsets around each power of 2
                for (int i = 0; i < 100; i++) {
                    BIGNUM* pattern = BN_dup(base);
                    int64_t offset = (gen() % 2000000) - 1000000;
                    if (offset >= 0) {
                        BN_add_word(pattern, offset);
                    } else {
                        BN_sub_word(pattern, -offset);
                    }
                    
                    if (BN_cmp(pattern, RANGE_START_BN) >= 0 && BN_cmp(pattern, RANGE_END_BN) <= 0) {
                        evolutionary_patterns.push_back(pattern);
                    } else {
                        BN_free(pattern);
                    }
                }
                BN_free(base);
            }
        } 
        else if (worker_id == 1) {
            // Arithmetic sequences from different starting points
            BIGNUM* start_point = BN_new();
            BN_set_bit(start_point, 70);
            BN_add_word(start_point, base_offset + batch_counter * 1000000);
            
            for (int step = 1; step <= 10000; step *= 10) {
                for (int mult = 1; mult <= 100; mult++) {
                    BIGNUM* pattern = BN_dup(start_point);
                    BN_add_word(pattern, mult * step);
                    
                    if (BN_cmp(pattern, RANGE_START_BN) >= 0 && BN_cmp(pattern, RANGE_END_BN) <= 0) {
                        evolutionary_patterns.push_back(pattern);
                    } else {
                        BN_free(pattern);
                        break;
                    }
                }
            }
            BN_free(start_point);
        }
        else if (worker_id == 2) {
            // Geometric-like patterns
            for (int base_val = 2; base_val <= 5; base_val++) {
                for (int multiplier = 1; multiplier <= 100; multiplier++) {
                    BIGNUM* pattern = BN_new();
                    BN_set_bit(pattern, 70);
                    BN_add_word(pattern, base_val * multiplier * 1000000 + worker_id * 100000);
                    
                    if (BN_cmp(pattern, RANGE_START_BN) >= 0 && BN_cmp(pattern, RANGE_END_BN) <= 0) {
                        evolutionary_patterns.push_back(pattern);
                    } else {
                        BN_free(pattern);
                    }
                }
            }
        }
        else {
            // Prime-like offsets and random mathematical structures
            for (int i = 0; i < 500; i++) {
                BIGNUM* pattern = BN_new();
                BN_set_bit(pattern, 70);
                
                // Mix of structured and random offsets
                uint64_t offset = base_offset + 
                                 (gen() % (1ULL << 25)) + 
                                 (batch_counter * 1234567) +
                                 (worker_id * 987654321);
                
                BN_add_word(pattern, offset);
                
                if (BN_cmp(pattern, RANGE_START_BN) >= 0 && BN_cmp(pattern, RANGE_END_BN) <= 0) {
                    evolutionary_patterns.push_back(pattern);
                } else {
                    BN_free(pattern);
                }
            }
        }
        
        // 🧠 TEST CURRENT BATCH
        for (BIGNUM* pattern : evolutionary_patterns) {
            if (FOUND) break;
            
            std::string hash = turbo_hash_bn(pattern);
            int match_length = 0;
            for (size_t i = 0; i < hash.size() && i < TARGET.size(); i++) {
                if (hash[i] == TARGET[i]) {
                    match_length++;
                } else {
                    break;
                }
            }
            
            // 🚨 LOWERED THRESHOLD FROM 10→6
            if (match_length >= 6) {
                GLOBAL_EVOLUTION.evolve_from_discovery(pattern, match_length, hash);
                if (match_length >= 8) {
                    char* pattern_hex = BN_bn2hex(pattern);
                    std::cout << "🔍 Worker " << worker_id << " pattern: " << match_length 
                              << " chars at " << std::string(pattern_hex).substr(0, 16) << "...\n";
                    OPENSSL_free(pattern_hex);
                }
            }
            
            if (hash == TARGET) {
                FOUND = true;
                char* found_key_hex = BN_bn2hex(pattern);
                std::cout << "\n🎯 PATTERN WORKER " << worker_id << " SOLVED THE PUZZLE!\n";
                std::cout << "Key: " << found_key_hex << "\n";
                OPENSSL_free(found_key_hex);
            }
            
            TOTAL_TRIED++;
        }
        
        // 🧠 CLEANUP BATCH
        for (BIGNUM* pattern : evolutionary_patterns) {
            BN_free(pattern);
        }
        
        // 🧠 BRIEF PAUSE TO ALLOW OTHER THREADS
        if (!FOUND && batch_counter % 10 == 0) {
            std::this_thread::sleep_for(std::chrono::milliseconds(50));
        }
    }
    
    BN_free(test_bn);
    EC_KEY_free(kk);
}

void generate_random_in_range(BIGNUM* result, const BIGNUM* min, const BIGNUM* max, std::mt19937_64& gen) {
    BIGNUM* range = BN_new();
    BN_CTX* ctx = BN_CTX_new();
    
    BN_sub(range, max, min);
    BN_add_word(range, 1);
    
    BN_rand_range(result, range);
    BN_add(result, result, min);
    
    BN_free(range);
    BN_CTX_free(ctx);
}

// 🧠 DEEP EXPLORATION WORKER
void deep_exploration_worker(int worker_id) {
    std::random_device rd;
    std::mt19937_64 gen(rd());
    
    EC_KEY* kk = EC_KEY_new_by_curve_name(NID_secp256k1);
    BIGNUM* key_bn = BN_new();
    const EC_GROUP* g = EC_KEY_get0_group(kk);
    EC_POINT* p = EC_POINT_new(g);
    unsigned char pub[65];
    unsigned char sha[SHA256_DIGEST_LENGTH];
    unsigned char rm[RIPEMD160_DIGEST_LENGTH];
    char out[41];
    static const char* hex = "0123456789abcdef";
    
    int exploration_cycle = 0;
    
    while (!FOUND) {
        // 🧠 DIVERSE STARTING POINTS ACROSS ENTIRE RANGE
        std::vector<BIGNUM*> exploration_bases;
        
        // Cover different regions: start, middle, end of 2^70-2^71 range
        for (int region = 0; region < 8; region++) {
            BIGNUM* base = BN_new();
            BN_set_bit(base, 70);
            
            // Spread bases across the entire range
            uint64_t region_size = (1ULL << 60); // Divide range into 8 regions
            uint64_t region_start = region * region_size;
            uint64_t worker_offset = worker_id * (region_size / 4);
            
            BN_add_word(base, region_start + worker_offset + (exploration_cycle * 1000000));
            
            // Ensure we're within range
            if (BN_cmp(base, RANGE_START_BN) >= 0 && BN_cmp(base, RANGE_END_BN) <= 0) {
                exploration_bases.push_back(base);
            } else {
                BN_free(base);
            }
        }
        
        // Add some completely random starting points
        for (int i = 0; i < 2; i++) {
            BIGNUM* random_base = BN_new();
            generate_random_in_range(random_base, RANGE_START_BN, RANGE_END_BN, gen);
            exploration_bases.push_back(random_base);
        }
        
        BIGNUM* exploration_center = exploration_bases[gen() % exploration_bases.size()];
        char* center_hex = BN_bn2hex(exploration_center);
        
        std::cout << "🏔️  Deep Exploration " << worker_id << " at " 
                  << std::string(center_hex).substr(0, 20) << "...\n";
        OPENSSL_free(center_hex);
        
        BIGNUM* exploration_radius = BN_new();
        BN_set_bit(exploration_radius, 18); // Focused deep search
        
        BIGNUM* min_explore = BN_new();
        BIGNUM* max_explore = BN_new();
        
        BN_sub(min_explore, exploration_center, exploration_radius);
        BN_add(max_explore, exploration_center, exploration_radius);
        
        if (BN_cmp(min_explore, RANGE_START_BN) < 0) BN_copy(min_explore, RANGE_START_BN);
        if (BN_cmp(max_explore, RANGE_END_BN) > 0) BN_copy(max_explore, RANGE_END_BN);
        
        // 🧠 INTENSIVE LOCAL SEARCH
        for (int i = 0; i < 50000 && !FOUND; i++) {
            BIGNUM* range = BN_new();
            BN_CTX* ctx = BN_CTX_new();
            
            BN_sub(range, max_explore, min_explore);
            BN_add_word(range, 1);
            
            BN_rand_range(key_bn, range);
            BN_add(key_bn, key_bn, min_explore);
            
            EC_KEY_set_private_key(kk, key_bn);
            EC_POINT_mul(g, p, key_bn, NULL, NULL, NULL);
            size_t pl = EC_POINT_point2oct(g, p, POINT_CONVERSION_COMPRESSED, pub, 65, NULL);
            SHA256(pub, pl, sha);
            RIPEMD160(sha, SHA256_DIGEST_LENGTH, rm);
            
            bool full_match = true;
            int match_length = 0;
            for (int j = 0; j < 20; j++) {
                char c1 = hex[rm[j] >> 4];
                char c2 = hex[rm[j] & 0xF];
                if (c1 != TARGET[j * 2] || c2 != TARGET[j * 2 + 1]) {
                    full_match = false;
                    match_length = j * 2;
                    if (c1 == TARGET[j * 2]) match_length++;
                    break;
                }
            }
            
            if (full_match) {
                FOUND = true;
                char* found_key_hex = BN_bn2hex(key_bn);
                std::cout << "\n🎉 DEEP EXPLORATION " << worker_id << " SOLVED THE PUZZLE!\n";
                std::cout << "Key: " << found_key_hex << "\n";
                OPENSSL_free(found_key_hex);
                break;
            }
            
            // 🚨 LOWERED THRESHOLD FROM 12→8
            if (match_length >= 8) {
                GLOBAL_EVOLUTION.evolve_from_discovery(key_bn, match_length, "");
                if (match_length >= 10) {
                    char* key_hex = BN_bn2hex(key_bn);
                    std::cout << "💎 Deep exploration gem: " << match_length << " char match!\n";
                    OPENSSL_free(key_hex);
                }
            }
            
            TOTAL_TRIED++;
            
            BN_free(range);
            BN_CTX_free(ctx);
        }
        
        // 🧠 CLEANUP
        for (BIGNUM* base : exploration_bases) {
            BN_free(base);
        }
        BN_free(exploration_radius);
        BN_free(min_explore);
        BN_free(max_explore);
        
        exploration_cycle++;
        
        if (!FOUND) {
            std::this_thread::sleep_for(std::chrono::seconds(3));
        }
    }
    
    EC_POINT_free(p);
    BN_free(key_bn);
    EC_KEY_free(kk);
}


// 🧠 MAIN EVOLUTIONARY SYSTEM
void ultra_adaptive_evolutionary_solver() {
    std::cout << "🧠 ULTRA-ADAPTIVE EVOLUTIONARY BITCOIN PUZZLE SOLVER 🧠\n";
    std::cout << "Target: " << TARGET << "\n";
    
    // 🧠 VERIFY HASHING IS WORKING
    std::cout << "🔍 Verifying hash function...\n";
    BIGNUM* test_key = BN_new();
    BN_set_bit(test_key, 70);
    BN_add_word(test_key, 12345);
    std::string test_hash = turbo_hash_bn(test_key);
    std::cout << "Test hash (first 20 chars): " << test_hash.substr(0, 20) << "...\n";
    std::cout << "Target hash (first 20 chars): " << TARGET.substr(0, 20) << "...\n";
    BN_free(test_key);
    
    auto start_time = std::chrono::steady_clock::now();
    
    std::vector<std::thread> evolutionary_threads;
    
    std::cout << "🚀 Starting " << N_THREADS << " evolutionary workers...\n";
    for (int i = 0; i < N_THREADS; i++) {
        evolutionary_threads.emplace_back(evolutionary_worker, i);
    }
    
    std::cout << "🔍 Starting 4 continuous pattern discovery workers...\n";
    for (int i = 0; i < 4; i++) {
        evolutionary_threads.emplace_back(parallel_pattern_discovery_worker, i, 4);
    }
    
    std::cout << "🏔️  Starting 2 improved deep exploration workers...\n";
    for (int i = 0; i < 2; i++) {
        evolutionary_threads.emplace_back(deep_exploration_worker, i);
    }
    
    uint64_t last_total = 0;
    auto last_evolution_report = std::chrono::steady_clock::now();
    auto last_save_time = std::chrono::steady_clock::now();
    
    while (!FOUND) {
        std::this_thread::sleep_for(std::chrono::seconds(30));
        
        auto now = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - start_time);
        uint64_t current_total = TOTAL_TRIED;
        uint64_t keys_per_sec = (current_total - last_total) / 30;
        
        std::cout << "\n📊 EVOLUTIONARY PROGRESS:\n";
        std::cout << "⏱️  Elapsed: " << elapsed.count() << "s | ";
        std::cout << "Rate: " << keys_per_sec << " keys/s | ";
        std::cout << "Total: " << current_total << " keys | ";
        std::cout << "Discoveries: " << GLOBAL_EVOLUTION.get_total_discoveries() << "\n";
        
        // 🚨 EMERGENCY EXPLORATION BOOST
        if (current_total > 1000000 && GLOBAL_EVOLUTION.get_total_discoveries() == 0) {
            GLOBAL_EVOLUTION.emergency_exploration_boost();
        }
        
        if (now - last_evolution_report > std::chrono::minutes(2)) {
            GLOBAL_EVOLUTION.print_evolutionary_report();
            last_evolution_report = now;
        }
        
        if (now - last_save_time > std::chrono::minutes(5)) {
            GLOBAL_EVOLUTION.save_evolutionary_state();
            last_save_time = now;
        }
        
        last_total = current_total;
    }
    
    for (auto& t : evolutionary_threads) {
        if (t.joinable()) t.join();
    }
    
    GLOBAL_EVOLUTION.save_evolutionary_state();
    std::cout << "\n🎉 EVOLUTION COMPLETED! Total keys: " << TOTAL_TRIED << "\n";
}

int main() {
    std::cout << "🧠 ULTRA-ADAPTIVE EVOLUTIONARY BITCOIN PUZZLE SOLVER 🧠\n";
    std::cout << "======================================================\n\n";
    
    init_bignum_constants();
    GLOBAL_EVOLUTION.load_evolutionary_state();
    
    auto start = std::chrono::steady_clock::now();
    ultra_adaptive_evolutionary_solver();
    auto end = std::chrono::steady_clock::now();
    
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    std::cout << "\n======================================================\n";
    std::cout << "Total evolutionary time: " << elapsed.count() << " seconds\n";
    
    if (FOUND) {
        std::cout << "🏆 ULTRA-ADAPTIVE EVOLUTION TRIUMPHS! 🏆\n";
    } else {
        std::cout << "🌱 Evolutionary journey continues... 🌱\n";
    }
    
    cleanup_bignum_constants();
    return 0;
}