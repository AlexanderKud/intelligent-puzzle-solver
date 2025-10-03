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
#include <openssl/sha.h>
#include <openssl/ripemd.h>
#include <openssl/ec.h>
#include <openssl/obj_mac.h>
#include <openssl/bn.h>

// 🧠 INTELLIGENT ADAPTIVE SEARCH WITH SELF-EVOLVING STRATEGY 🧠
// Using BIGNUM for 71-bit range since uint64_t can't handle 2^70+
constexpr int N_THREADS = 12;
constexpr int CHUNK_SIZE = 10000;

std::atomic<bool> FOUND(false);
std::mutex intelligence_mutex;
std::atomic<uint64_t> TOTAL_TRIED(0);

const std::string TARGET = "f6f5431d25bbf7b12e8add9af5e3475c44a0a5b8";

// Big number constants for 2^70 and 2^71
BIGNUM* HACKY_RANGE_START_BN = nullptr;
BIGNUM* HACKY_RANGE_END_BN = nullptr;

void init_bignum_constants() {
    HACKY_RANGE_START_BN = BN_new();
    HACKY_RANGE_END_BN = BN_new();
    
    // Set 2^70
    BN_zero(HACKY_RANGE_START_BN);
    BN_set_bit(HACKY_RANGE_START_BN, 70);
    
    // Set 2^71 - 1
    BN_zero(HACKY_RANGE_END_BN);
    BN_set_bit(HACKY_RANGE_END_BN, 71);
    BN_sub_word(HACKY_RANGE_END_BN, 1);
}

void cleanup_bignum_constants() {
    if (HACKY_RANGE_START_BN) BN_free(HACKY_RANGE_START_BN);
    if (HACKY_RANGE_END_BN) BN_free(HACKY_RANGE_END_BN);
}

// 🧠 INTELLIGENCE SYSTEM - TRACKS PATTERNS AND ADAPTS
struct SearchIntelligence {
    std::map<int, uint64_t> prefix_matches; // prefix_length -> count
    std::map<std::string, int> promising_ranges; // key_hex -> match_length
    std::vector<std::pair<std::string, int>> best_matches; // (key_hex, match_length)
    int adaptive_threshold = 8; // Start focusing on 8+ character matches
    double focus_intensity = 0.3; // % of threads focusing on promising areas
    
    void record_match(const std::string& key_hex, int match_length, const std::string& hash) {
        std::lock_guard<std::mutex> lock(intelligence_mutex);
        
        prefix_matches[match_length]++;
        
        // Add to promising ranges if it's a good match
        if (match_length >= adaptive_threshold) {
            promising_ranges[key_hex] = match_length;
            
            // Keep best matches sorted
            best_matches.push_back({key_hex, match_length});
            std::sort(best_matches.begin(), best_matches.end(), 
                     [](const auto& a, const auto& b) { return a.second > b.second; });
            
            // Keep only top 20 best matches
            if (best_matches.size() > 20) {
                best_matches.pop_back();
            }
            
            // Adapt threshold based on findings
            if (match_length > adaptive_threshold + 2) {
                adaptive_threshold = match_length - 1;
                std::cout << "🧠 INTELLIGENCE: Raised threshold to " << adaptive_threshold 
                          << " (found match of length " << match_length << ")\n";
            }
        }
        
        // Increase focus on promising areas as we learn
        if (prefix_matches[match_length] % 100 == 0) {
            focus_intensity = std::min(0.8, 0.3 + (prefix_matches[match_length] / 1000.0) * 0.1);
        }
    }
    
    std::vector<std::string> get_promising_centers() {
        std::lock_guard<std::mutex> lock(intelligence_mutex);
        std::vector<std::string> centers;
        for (const auto& [center, length] : promising_ranges) {
            centers.push_back(center);
        }
        return centers;
    }
    
    void print_intelligence() {
        std::lock_guard<std::mutex> lock(intelligence_mutex);
        std::cout << "\n🧠 SEARCH INTELLIGENCE REPORT:\n";
        std::cout << "Adaptive Threshold: " << adaptive_threshold << "\n";
        std::cout << "Focus Intensity: " << (focus_intensity * 100) << "%\n";
        std::cout << "Promising Ranges: " << promising_ranges.size() << "\n";
        std::cout << "Prefix Match Distribution:\n";
        for (const auto& [length, count] : prefix_matches) {
            if (count > 0) {
                std::cout << "  " << length << " chars: " << count << " matches\n";
            }
        }
        if (!best_matches.empty()) {
            std::cout << "Best Matches:\n";
            for (size_t i = 0; i < std::min(best_matches.size(), size_t(5)); i++) {
                std::cout << "  " << best_matches[i].second << " chars at key " 
                          << best_matches[i].first.substr(0, 16) << "...\n";
            }
        }
    }
};

SearchIntelligence GLOBAL_INTELLIGENCE;

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

// Generate random BIGNUM in range [min, max]
void generate_random_bignum_in_range(BIGNUM* result, const BIGNUM* min, const BIGNUM* max, std::mt19937_64& gen) {
    BIGNUM* range = BN_new();
    BN_CTX* ctx = BN_CTX_new();
    
    // Calculate range = max - min + 1
    BN_sub(range, max, min);
    BN_add_word(range, 1);
    
    // Generate random number in range
    BN_rand_range(result, range);
    
    // Add min to get number in [min, max]
    BN_add(result, result, min);
    
    BN_free(range);
    BN_CTX_free(ctx);
}

// Convert BIGNUM to hex string
std::string bn_to_hex(const BIGNUM* bn) {
    char* hex = BN_bn2hex(bn);
    std::string result(hex);
    OPENSSL_free(hex);
    return result;
}

// Convert hex string to BIGNUM
BIGNUM* hex_to_bn(const std::string& hex) {
    BIGNUM* bn = BN_new();
    BN_hex2bn(&bn, hex.c_str());
    return bn;
}

// 🧠 ADAPTIVE WORKER - CAN SWITCH BETWEEN RANDOM AND FOCUSED SEARCH
void adaptive_worker(int thread_id) {
    std::random_device rd;
    std::mt19937_64 gen(rd());
    
    // Thread-local crypto context
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
    auto last_switch = std::chrono::steady_clock::now();
    bool in_focused_mode = false;
    std::string focus_center_hex;
    BIGNUM* focus_radius_bn = BN_new();
    BN_set_word(focus_radius_bn, 1ULL << 30); // Start with 2^30 radius
    
    while (!FOUND) {
        // 🧠 DECIDE SEARCH MODE BASED ON INTELLIGENCE
        auto now = std::chrono::steady_clock::now();
        bool should_switch = (now - last_switch > std::chrono::seconds(30));
        
        if (should_switch || !in_focused_mode) {
            double focus_chance = GLOBAL_INTELLIGENCE.focus_intensity;
            auto promising_centers = GLOBAL_INTELLIGENCE.get_promising_centers();
            
            if (!promising_centers.empty() && (gen() % 1000 < static_cast<int>(focus_chance * 1000))) {
                // Switch to focused mode
                in_focused_mode = true;
                focus_center_hex = promising_centers[gen() % promising_centers.size()];
                std::lock_guard<std::mutex> lock(intelligence_mutex);
                std::cout << "🎯 Thread " << thread_id << " focusing around key " 
                          << focus_center_hex.substr(0, 16) << "...\n";
            } else {
                // Switch to exploration mode
                in_focused_mode = false;
            }
            last_switch = now;
        }
        
        // 🧠 GENERATE KEY BASED ON CURRENT MODE
        if (in_focused_mode) {
            // Generate in focused region around promising key
            BIGNUM* center_bn = hex_to_bn(focus_center_hex);
            BIGNUM* min_bn = BN_new();
            BIGNUM* max_bn = BN_new();
            
            // Calculate min = center - radius, max = center + radius
            BN_sub(min_bn, center_bn, focus_radius_bn);
            BN_add(max_bn, center_bn, focus_radius_bn);
            
            // Clamp to global range
            if (BN_cmp(min_bn, HACKY_RANGE_START_BN) < 0) {
                BN_copy(min_bn, HACKY_RANGE_START_BN);
            }
            if (BN_cmp(max_bn, HACKY_RANGE_END_BN) > 0) {
                BN_copy(max_bn, HACKY_RANGE_END_BN);
            }
            
            generate_random_bignum_in_range(key_bn, min_bn, max_bn, gen);
            
            BN_free(center_bn);
            BN_free(min_bn);
            BN_free(max_bn);
        } else {
            // Random exploration in full range
            generate_random_bignum_in_range(key_bn, HACKY_RANGE_START_BN, HACKY_RANGE_END_BN, gen);
        }
        
        // 🧠 TURBO HASH
        EC_KEY_set_private_key(kk, key_bn);
        EC_POINT_mul(g, p, key_bn, NULL, NULL, NULL);
        size_t pl = EC_POINT_point2oct(g, p, POINT_CONVERSION_COMPRESSED, pub, 65, NULL);
        SHA256(pub, pl, sha);
        RIPEMD160(sha, SHA256_DIGEST_LENGTH, rm);
        
        // 🧠 SMART COMPARE WITH PATTERN DETECTION
        bool full_match = true;
        int current_match = 0;
        
        for (int i = 0; i < 20; i++) {
            char c1 = hex_chars[rm[i] >> 4];
            char c2 = hex_chars[rm[i] & 0xF];
            
            if (c1 != TARGET[i * 2] || c2 != TARGET[i * 2 + 1]) {
                full_match = false;
                current_match = i * 2; // Each byte gives 2 hex chars
                if (c1 == TARGET[i * 2]) {
                    current_match++; // First char matched
                }
                break;
            }
        }
        
        if (full_match) {
            FOUND = true;
            std::string found_key_hex = bn_to_hex(key_bn);
            std::lock_guard<std::mutex> lock(intelligence_mutex);
            std::cout << "\n🎉 FULL MATCH FOUND BY ADAPTIVE THREAD " << thread_id << "!\n";
            std::cout << "Key: " << found_key_hex << "\n";
            break;
        }
        
        // 🧠 RECORD INTELLIGENCE ABOUT PARTIAL MATCHES
        if (current_match >= 6) { // Record matches of 6+ characters
            std::string key_hex = bn_to_hex(key_bn);
            GLOBAL_INTELLIGENCE.record_match(key_hex, current_match, "");
        }
        
        local_tested++;
        TOTAL_TRIED++;
        
        // 🧠 OCCASIONAL INTELLIGENCE REPORTING
        if (local_tested % 500000 == 0) {
            std::lock_guard<std::mutex> lock(intelligence_mutex);
            std::cout << "🧠 Thread " << thread_id << " tested " << local_tested 
                      << " keys (" << (in_focused_mode ? "FOCUSED" : "EXPLORING") << ")\n";
            if (current_match >= 8) {
                std::string key_hex = bn_to_hex(key_bn);
                std::cout << "   Found promising match: " << current_match 
                          << " characters at key " << key_hex.substr(0, 16) << "...\n";
            }
            local_tested = 0;
        }
    }
    
    EC_POINT_free(p);
    BN_free(key_bn);
    BN_free(focus_radius_bn);
    EC_KEY_free(kk);
}

// 🧠 PATTERN-ANALYSIS WORKER - LOOKS FOR MATHEMATICAL PATTERNS
void pattern_analysis_worker() {
    std::cout << "🔍 Pattern Analysis Worker started...\n";
    
    EC_KEY* kk = EC_KEY_new_by_curve_name(NID_secp256k1);
    BIGNUM* key_bn = BN_new();
    
    // Test common mathematical patterns within our range
    std::vector<BIGNUM*> test_patterns;
    
    // Powers of 2 near our range
    for (int exp = 70; exp < 72; exp++) {
        BIGNUM* pattern = BN_new();
        BN_set_bit(pattern, exp);
        test_patterns.push_back(pattern);
        
        // pattern - 1
        BIGNUM* pattern_minus = BN_dup(pattern);
        BN_sub_word(pattern_minus, 1);
        test_patterns.push_back(pattern_minus);
        
        // pattern + 1
        BIGNUM* pattern_plus = BN_dup(pattern);
        BN_add_word(pattern_plus, 1);
        test_patterns.push_back(pattern_plus);
    }
    
    // Test all patterns
    for (BIGNUM* pattern : test_patterns) {
        if (FOUND) break;
        
        // Ensure pattern is within range
        if (BN_cmp(pattern, HACKY_RANGE_START_BN) < 0 || BN_cmp(pattern, HACKY_RANGE_END_BN) > 0) {
            continue;
        }
        
        std::string hash = turbo_hash_bn(pattern);
        int match_length = 0;
        for (size_t i = 0; i < hash.size() && i < TARGET.size(); i++) {
            if (hash[i] == TARGET[i]) {
                match_length++;
            } else {
                break;
            }
        }
        
        if (match_length >= 10) {
            std::string pattern_hex = bn_to_hex(pattern);
            GLOBAL_INTELLIGENCE.record_match(pattern_hex, match_length, hash);
            std::cout << "🔍 Pattern match: " << match_length 
                      << " chars at mathematical pattern " << pattern_hex << "\n";
        }
        
        if (hash == TARGET) {
            FOUND = true;
            std::string found_key_hex = bn_to_hex(pattern);
            std::cout << "\n🎯 PATTERN SOLVED THE PUZZLE!\n";
            std::cout << "Key: " << found_key_hex << "\n";
            break;
        }
    }
    
    // Cleanup patterns
    for (BIGNUM* pattern : test_patterns) {
        BN_free(pattern);
    }
    
    BN_free(key_bn);
    EC_KEY_free(kk);
}

// 🧠 NEIGHBORHOOD EXPLORER - DEEPLY EXPLORES AROUND PROMISING KEYS
void neighborhood_explorer(int worker_id) {
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
    
    while (!FOUND) {
        // Get promising centers to explore
        auto centers = GLOBAL_INTELLIGENCE.get_promising_centers();
        if (centers.empty()) {
            std::this_thread::sleep_for(std::chrono::seconds(5));
            continue;
        }
        
        std::string center_hex = centers[gen() % centers.size()];
        BIGNUM* center_bn = hex_to_bn(center_hex);
        
        BIGNUM* radius_bn = BN_new();
        BN_set_word(radius_bn, 1ULL << 20); // Explore 2^20 range around center
        
        BIGNUM* min_bn = BN_new();
        BIGNUM* max_bn = BN_new();
        
        // Calculate min = center - radius, max = center + radius
        BN_sub(min_bn, center_bn, radius_bn);
        BN_add(max_bn, center_bn, radius_bn);
        
        // Clamp to global range
        if (BN_cmp(min_bn, HACKY_RANGE_START_BN) < 0) {
            BN_copy(min_bn, HACKY_RANGE_START_BN);
        }
        if (BN_cmp(max_bn, HACKY_RANGE_END_BN) > 0) {
            BN_copy(max_bn, HACKY_RANGE_END_BN);
        }
        
        std::cout << "🏘️  Neighborhood Explorer " << worker_id 
                  << " deeply exploring around key " << center_hex.substr(0, 16) << "...\n";
        
        for (int i = 0; i < 100000 && !FOUND; i++) {
            generate_random_bignum_in_range(key_bn, min_bn, max_bn, gen);
            
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
                std::string found_key_hex = bn_to_hex(key_bn);
                std::cout << "\n🎉 NEIGHBORHOOD EXPLORER FOUND IT!\n";
                std::cout << "Key: " << found_key_hex << "\n";
                break;
            }
            
            if (match_length >= 12) {
                std::string key_hex = bn_to_hex(key_bn);
                GLOBAL_INTELLIGENCE.record_match(key_hex, match_length, "");
                if (match_length >= 16) {
                    std::cout << "🔥 HOT ZONE: " << match_length 
                              << " char match in neighborhood exploration!\n";
                }
            }
            
            TOTAL_TRIED++;
        }
        
        BN_free(center_bn);
        BN_free(radius_bn);
        BN_free(min_bn);
        BN_free(max_bn);
    }
    
    EC_POINT_free(p);
    BN_free(key_bn);
    EC_KEY_free(kk);
}

// 🧠 SELF-EVOLVING ATTACK SYSTEM
void intelligent_71bit_attack() {
    std::cout << "🧠 LAUNCHING INTELLIGENT SELF-EVOLVING 71-BIT ATTACK 🧠\n";
    std::cout << "Target: " << TARGET << "\n";
    std::cout << "Range: 2^70 to 2^71\n";
    std::cout << "Threads: " << N_THREADS << "\n";
    std::cout << "This system will ADAPT and EVOLVE its strategy based on partial successes!\n\n";
    
    auto start_time = std::chrono::steady_clock::now();
    
    std::vector<std::thread> threads;
    
    // Start adaptive workers (can switch between exploration and focused search)
    std::cout << "🚀 Starting " << N_THREADS << " adaptive workers...\n";
    for (int i = 0; i < N_THREADS; i++) {
        threads.emplace_back(adaptive_worker, i);
    }
    
    // Start pattern analysis
    std::cout << "🔍 Starting pattern analysis...\n";
    threads.emplace_back(pattern_analysis_worker);
    
    // Start neighborhood explorers
    std::cout << "🏘️  Starting neighborhood explorers...\n";
    for (int i = 0; i < 2; i++) {
        threads.emplace_back(neighborhood_explorer, i);
    }
    
    // Intelligence monitoring and reporting
    uint64_t last_total = 0;
    auto last_intelligence_report = std::chrono::steady_clock::now();
    
    while (!FOUND) {
        std::this_thread::sleep_for(std::chrono::seconds(30));
        
        auto now = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - start_time);
        uint64_t current_total = TOTAL_TRIED;
        uint64_t keys_per_sec = (current_total - last_total) / 30;
        
        std::cout << "\n📊 PROGRESS UPDATE:\n";
        std::cout << "⏱️  Elapsed: " << elapsed.count() << "s | ";
        std::cout << "Rate: " << keys_per_sec << " keys/s | ";
        std::cout << "Total: " << current_total << " keys\n";
        
        // Regular intelligence reports
        if (now - last_intelligence_report > std::chrono::minutes(5)) {
            GLOBAL_INTELLIGENCE.print_intelligence();
            last_intelligence_report = now;
        }
        
        last_total = current_total;
        
        // Adaptive timeout based on progress
        if (elapsed.count() > 3600 && keys_per_sec < 1000) { // 1 hour with low rate
            std::cout << "💤 Low detection rate - the puzzle might be in a sparse region.\n";
            std::cout << "🧠 Increasing focus intensity to find more promising areas...\n";
            GLOBAL_INTELLIGENCE.focus_intensity = std::min(0.9, GLOBAL_INTELLIGENCE.focus_intensity + 0.1);
        }
        
        if (elapsed.count() > 86400) { // 24 hours
            std::cout << "🧠 24H Analysis: This is a VERY hard puzzle!\n";
            std::cout << "Keys tested: " << current_total << " / ~1.2e21 possible\n";
            std::cout << "At current rate, estimated time: " 
                      << (1.2e21 / keys_per_sec / 86400 / 365) << " years\n";
            break;
        }
    }
    
    // Wait for all threads
    for (auto& t : threads) {
        if (t.joinable()) t.join();
    }
    
    auto end_time = std::chrono::steady_clock::now();
    auto total_elapsed = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    
    if (FOUND) {
        std::cout << "\n🎉 INTELLIGENT SYSTEM SOLVED THE 71-BIT PUZZLE! 🎉\n";
        std::cout << "Time: " << total_elapsed.count() << " seconds\n";
        std::cout << "Total Keys Tested: " << TOTAL_TRIED << "\n";
        std::cout << "Final Intelligence Report:\n";
        GLOBAL_INTELLIGENCE.print_intelligence();
    } else {
        std::cout << "\n💀 Intelligent system couldn't crack it this time...\n";
        std::cout << "But we learned a lot! Final intelligence:\n";
        GLOBAL_INTELLIGENCE.print_intelligence();
    }
}

int main() {
    std::cout << "🧠 SELF-EVOLVING INTELLIGENT BITCOIN PUZZLE SOLVER 🧠\n";
    std::cout << "====================================================\n\n";
    
    // Initialize OpenSSL and big number constants
    init_bignum_constants();
    
    auto start = std::chrono::steady_clock::now();
    
    intelligent_71bit_attack();
    
    auto end = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    
    std::cout << "\n====================================================\n";
    std::cout << "Total execution time: " << elapsed.count() << " seconds\n";
    
    if (FOUND) {
        std::cout << "🏆 INTELLIGENCE TRIUMPHS! 🏆\n";
    } else {
        std::cout << "🤖 The AI will continue learning... 🤖\n";
    }
    
    // Cleanup
    cleanup_bignum_constants();
    
    return 0;
}