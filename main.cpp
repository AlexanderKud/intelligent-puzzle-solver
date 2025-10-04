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
#include <openssl/sha.h>
#include <openssl/ripemd.h>
#include <openssl/ec.h>
#include <openssl/obj_mac.h>
#include <openssl/bn.h>

// 🧠 REAL-TIME INTELLIGENT FILTERING WITH CORNER DETECTION 🧠
// constexpr int N_THREADS = 16;
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

// 🧠 REAL-TIME PATTERN ANALYSIS ENGINE
class RealTimeIntelligence {
private:
    struct PatternCluster {
        BIGNUM* center;
        int match_strength;
        uint64_t sample_count;
        double success_density;
        std::vector<int> bit_patterns;
        
        PatternCluster(BIGNUM* center_bn, int strength) {
            center = BN_dup(center_bn);
            match_strength = strength;
            sample_count = 1;
            success_density = strength / 100.0;
        }
        
        ~PatternCluster() {
            BN_free(center);
        }
    };
    
    std::vector<PatternCluster*> active_clusters;
    std::map<std::string, double> corner_performance; // corner_hash -> success_rate
    std::vector<std::pair<BIGNUM*, int>> recent_matches; // (key, match_length)
    std::map<int, uint64_t> bit_position_correlation; // bit_pos -> match_correlation
    
    // Bit pattern analysis
    static const int BIT_ANALYSIS_DEPTH = 256;
    std::vector<std::vector<double>> bit_transition_matrix; // [current_pattern][next_bit] -> probability
    
public:
    RealTimeIntelligence() {
        // Initialize bit transition matrix
        bit_transition_matrix.resize(BIT_ANALYSIS_DEPTH, std::vector<double>(2, 0.5));
    }
    
    ~RealTimeIntelligence() {
        for (auto cluster : active_clusters) {
            delete cluster;
        }
    }
    
    // 🧠 ANALYZE EACH KEY IN REAL-TIME AND UPDATE FILTERING STRATEGY
    void analyze_key_pattern(BIGNUM* key_bn, int match_length, const std::string& hash_result) {
        std::lock_guard<std::mutex> lock(intelligence_mutex);
        
        // 1. Update recent matches queue
        recent_matches.push_back({BN_dup(key_bn), match_length});
        if (recent_matches.size() > 1000) {
            BN_free(recent_matches[0].first);
            recent_matches.erase(recent_matches.begin());
        }
        
        // 2. Bit-level pattern analysis
        analyze_bit_patterns(key_bn, match_length);
        
        // 3. Cluster analysis for promising regions
        update_clusters(key_bn, match_length);
        
        // 4. Corner detection in the search space
        detect_promising_corners();
        
        // 5. Adaptive bit transition probabilities
        update_transition_matrix(key_bn, match_length);
    }
    
    // 🧠 BIT-LEVEL PATTERN ANALYSIS
    void analyze_bit_patterns(BIGNUM* key_bn, int match_length) {
        // Convert BIGNUM to bit pattern
        std::bitset<256> bit_pattern;
        for (int i = 0; i < 256; i++) {
            if (BN_is_bit_set(key_bn, i)) {
                bit_pattern.set(i);
            }
        }
        
        // Analyze correlation between bit positions and match success
        if (match_length >= 12) {
            for (int i = 0; i < 71; i++) { // Only analyze relevant bits for 71-bit range
                if (bit_pattern[i]) {
                    bit_position_correlation[i]++;
                }
            }
        }
    }
    
    // 🧠 CLUSTER PROMISING REGIONS
    void update_clusters(BIGNUM* key_bn, int match_length) {
        if (match_length < 10) return;
        
        bool added_to_existing = false;
        
        // Try to add to existing cluster
        for (auto cluster : active_clusters) {
            BIGNUM* distance = BN_new();
            BN_CTX* ctx = BN_CTX_new();
            
            BN_sub(distance, key_bn, cluster->center);
            
            // If close to cluster center, update cluster
            if (BN_num_bits(distance) < 20) { // Within 2^20 range
                cluster->sample_count++;
                cluster->match_strength = std::max(cluster->match_strength, match_length);
                cluster->success_density = (cluster->success_density * (cluster->sample_count - 1) + 
                                          (match_length / 100.0)) / cluster->sample_count;
                added_to_existing = true;
            }
            
            BN_free(distance);
            BN_CTX_free(ctx);
        }
        
        // Create new cluster if no suitable existing cluster
        if (!added_to_existing && match_length >= 12) {
            active_clusters.push_back(new PatternCluster(key_bn, match_length));
            
            // Keep only top clusters
            if (active_clusters.size() > 20) {
                std::sort(active_clusters.begin(), active_clusters.end(),
                         [](const PatternCluster* a, const PatternCluster* b) {
                             return a->success_density > b->success_density;
                         });
                delete active_clusters.back();
                active_clusters.pop_back();
            }
        }
    }
    
    // 🧠 DETECT PROMISING CORNERS OF THE SEARCH SPACE
    void detect_promising_corners() {
        if (recent_matches.size() < 100) return;
        
        // Analyze recent matches to find promising corners
        std::map<std::string, std::pair<uint64_t, uint64_t>> corner_stats; // corner_id -> (total, successes)
        
        for (const auto& [key_bn, match_length] : recent_matches) {
            // Define "corners" based on bit patterns in upper bits
            std::string corner_id = get_corner_identifier(key_bn);
            
            corner_stats[corner_id].first++;
            if (match_length >= 12) {
                corner_stats[corner_id].second++;
            }
        }
        
        // Update corner performance
        for (const auto& [corner, stats] : corner_stats) {
            double success_rate = static_cast<double>(stats.second) / stats.first;
            corner_performance[corner] = success_rate;
        }
    }
    
    // 🧠 GET CORNER IDENTIFIER BASED ON BIT PATTERNS
    std::string get_corner_identifier(BIGNUM* key_bn) {
        // Use upper 16 bits to define "corners" of the search space
        std::string corner_id;
        for (int i = 55; i < 71; i++) { // Upper 16 bits of 71-bit range
            corner_id += BN_is_bit_set(key_bn, i) ? '1' : '0';
        }
        return corner_id;
    }
    
    // 🧠 UPDATE BIT TRANSITION PROBABILITIES
    void update_transition_matrix(BIGNUM* key_bn, int match_length) {
        // Analyze transitions between successful bit patterns
        // This is simplified - in reality you'd track sequences of successful keys
        double learning_rate = 0.01 * (match_length - 8); // Learn more from better matches
        
        // Update probabilities based on successful keys
        for (int i = 0; i < 70; i++) {
            int current_bit = BN_is_bit_set(key_bn, i);
            int next_bit = BN_is_bit_set(key_bn, i + 1);
            
            // Smooth update of transition probabilities
            bit_transition_matrix[current_bit][next_bit] = 
                (1 - learning_rate) * bit_transition_matrix[current_bit][next_bit] + 
                learning_rate * 1.0;
        }
    }
    
    // 🧠 GENERATE INTELLIGENT NEXT KEY BASED ON LEARNED PATTERNS
    void generate_intelligent_key(BIGNUM* result, std::mt19937_64& gen, 
                                 std::uniform_real_distribution<>& dis) {
        std::lock_guard<std::mutex> lock(intelligence_mutex);
        
        // Strategy selection based on learned patterns
        double strategy = dis(gen);
        
        if (strategy < 0.4 && !active_clusters.empty()) {
            // Cluster-based generation (40% probability)
            generate_from_clusters(result, gen);
        } else if (strategy < 0.7 && !corner_performance.empty()) {
            // Corner-based generation (30% probability)
            generate_from_corners(result, gen);
        } else if (strategy < 0.9) {
            // Bit-pattern-based generation (20% probability)
            generate_from_bit_patterns(result, gen);
        } else {
            // Pure exploration (10% probability)
            generate_random_in_range(result, RANGE_START_BN, RANGE_END_BN, gen);
        }
    }
    
    // 🧠 GENERATE KEY FROM PROMISING CLUSTERS
    void generate_from_clusters(BIGNUM* result, std::mt19937_64& gen) {
        // Weighted selection by cluster success density
        double total_weight = 0;
        for (auto cluster : active_clusters) {
            total_weight += cluster->success_density;
        }
        
        double selection = static_cast<double>(gen() % 10000) / 10000.0 * total_weight;
        double current = 0;
        
        PatternCluster* selected_cluster = nullptr;
        for (auto cluster : active_clusters) {
            current += cluster->success_density;
            if (current >= selection) {
                selected_cluster = cluster;
                break;
            }
        }
        
        if (selected_cluster) {
            // Generate near cluster center with adaptive radius
            BIGNUM* radius = BN_new();
            BN_set_bit(radius, 20 - selected_cluster->match_strength); // Smaller radius for better clusters
            
            BIGNUM* min_val = BN_new();
            BIGNUM* max_val = BN_new();
            
            BN_sub(min_val, selected_cluster->center, radius);
            BN_add(max_val, selected_cluster->center, radius);
            
            // Clamp to range
            if (BN_cmp(min_val, RANGE_START_BN) < 0) BN_copy(min_val, RANGE_START_BN);
            if (BN_cmp(max_val, RANGE_END_BN) > 0) BN_copy(max_val, RANGE_END_BN);
            
            generate_random_in_range(result, min_val, max_val, gen);
            
            BN_free(radius);
            BN_free(min_val);
            BN_free(max_val);
        }
    }
    
    // 🧠 GENERATE KEY FROM HIGH-PERFORMANCE CORNERS
    void generate_from_corners(BIGNUM* result, std::mt19937_64& gen) {
        // Find best performing corners
        std::vector<std::pair<std::string, double>> sorted_corners(
            corner_performance.begin(), corner_performance.end());
        std::sort(sorted_corners.begin(), sorted_corners.end(),
                 [](const auto& a, const auto& b) { return a.second > b.second; });
        
        if (sorted_corners.empty()) return;
        
        // Select from top 3 corners
        const auto& selected_corner = sorted_corners[gen() % std::min(3, (int)sorted_corners.size())];
        
        // Create base key from corner pattern
        BN_zero(result);
        for (size_t i = 0; i < selected_corner.first.size(); i++) {
            if (selected_corner.first[i] == '1') {
                BN_set_bit(result, 55 + i); // Set upper bits according to corner pattern
            }
        }
        
        // Add random lower bits
        BIGNUM* random_lower = BN_new();
        BN_rand(random_lower, 55, BN_RAND_TOP_ANY, BN_RAND_BOTTOM_ANY); // 55 lower bits
        BN_add(result, result, random_lower);
        BN_free(random_lower);
    }
    
    // 🧠 GENERATE KEY USING LEARNED BIT PATTERNS
    void generate_from_bit_patterns(BIGNUM* result, std::mt19937_64& gen) {
        BN_zero(result);
        
        int current_bit = gen() % 2;
        for (int i = 0; i < 71; i++) {
            // Use learned transition probabilities
            double prob = bit_transition_matrix[current_bit][1];
            if (static_cast<double>(gen() % 10000) / 10000.0 < prob) {
                BN_set_bit(result, i);
                current_bit = 1;
            } else {
                current_bit = 0;
            }
        }
        
        // Ensure within range
        if (BN_cmp(result, RANGE_START_BN) < 0) BN_copy(result, RANGE_START_BN);
        if (BN_cmp(result, RANGE_END_BN) > 0) BN_copy(result, RANGE_END_BN);
    }
    
    // 🧠 UTILITY: GENERATE RANDOM IN RANGE
    void generate_random_in_range(BIGNUM* result, const BIGNUM* min, const BIGNUM* max, 
                                 std::mt19937_64& gen) {
        BIGNUM* range = BN_new();
        BN_CTX* ctx = BN_CTX_new();
        
        BN_sub(range, max, min);
        BN_add_word(range, 1);
        
        BN_rand_range(result, range);
        BN_add(result, result, min);
        
        BN_free(range);
        BN_CTX_free(ctx);
    }
    
    // 🧠 GET INTELLIGENCE REPORT
    void print_intelligence_report() {
        std::lock_guard<std::mutex> lock(intelligence_mutex);
        
        std::cout << "\n🧠 REAL-TIME INTELLIGENCE REPORT:\n";
        std::cout << "Active Clusters: " << active_clusters.size() << "\n";
        std::cout << "Tracked Corners: " << corner_performance.size() << "\n";
        std::cout << "Recent Matches Analyzed: " << recent_matches.size() << "\n";
        
        if (!active_clusters.empty()) {
            std::cout << "Top Clusters:\n";
            for (size_t i = 0; i < std::min(active_clusters.size(), size_t(3)); i++) {
                std::cout << "  Cluster " << i << ": strength=" << active_clusters[i]->match_strength
                          << ", density=" << active_clusters[i]->success_density 
                          << ", samples=" << active_clusters[i]->sample_count << "\n";
            }
        }
        
        if (!corner_performance.empty()) {
            std::cout << "Top Performing Corners:\n";
            std::vector<std::pair<std::string, double>> sorted_corners(
                corner_performance.begin(), corner_performance.end());
            std::sort(sorted_corners.begin(), sorted_corners.end(),
                     [](const auto& a, const auto& b) { return a.second > b.second; });
            
            for (size_t i = 0; i < std::min(sorted_corners.size(), size_t(3)); i++) {
                std::cout << "  Corner " << sorted_corners[i].first << ": success_rate=" 
                          << sorted_corners[i].second << "\n";
            }
        }
        
        // Show bit correlation for top positions
        if (!bit_position_correlation.empty()) {
            std::cout << "Top Bit Correlations:\n";
            std::vector<std::pair<int, uint64_t>> sorted_bits(
                bit_position_correlation.begin(), bit_position_correlation.end());
            std::sort(sorted_bits.begin(), sorted_bits.end(),
                     [](const auto& a, const auto& b) { return a.second > b.second; });
            
            for (size_t i = 0; i < std::min(sorted_bits.size(), size_t(5)); i++) {
                std::cout << "  Bit " << sorted_bits[i].first << ": correlation=" 
                          << sorted_bits[i].second << "\n";
            }
        }
    }
};

RealTimeIntelligence GLOBAL_INTELLIGENCE;

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

// 🧠 INTELLIGENT WORKER WITH REAL-TIME FILTERING
void intelligent_worker(int thread_id) {
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
        // 🧠 GENERATE KEY USING REAL-TIME INTELLIGENCE
        GLOBAL_INTELLIGENCE.generate_intelligent_key(key_bn, gen, dis);
        
        // 🧠 TURBO HASH
        EC_KEY_set_private_key(kk, key_bn);
        EC_POINT_mul(g, p, key_bn, NULL, NULL, NULL);
        size_t pl = EC_POINT_point2oct(g, p, POINT_CONVERSION_COMPRESSED, pub, 65, NULL);
        SHA256(pub, pl, sha);
        RIPEMD160(sha, SHA256_DIGEST_LENGTH, rm);
        
        // 🧠 SMART COMPARE WITH PATTERN ANALYSIS
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
            std::string found_key_hex = BN_bn2hex(key_bn);
            std::lock_guard<std::mutex> lock(intelligence_mutex);
            std::cout << "\n🎉 FULL MATCH FOUND BY INTELLIGENT WORKER " << thread_id << "!\n";
            std::cout << "Key: " << found_key_hex << "\n";
            break;
        }
        
        // 🧠 REAL-TIME PATTERN ANALYSIS AND FILTERING UPDATE
        if (current_match >= 8) {
            GLOBAL_INTELLIGENCE.analyze_key_pattern(key_bn, current_match, "");
            
            if (current_match >= 14) {
                std::string key_hex = BN_bn2hex(key_bn);
                std::cout << "🔥 Thread " << thread_id << " found strong pattern: " 
                          << current_match << " chars at " << key_hex.substr(0, 20) << "...\n";
            }
        }
        
        local_tested++;
        TOTAL_TRIED++;
        
        // 🧠 PERIODIC REPORTING
        if (local_tested % 100000 == 0) {
            auto now = std::chrono::steady_clock::now();
            if (now - last_report > std::chrono::seconds(30)) {
                std::lock_guard<std::mutex> lock(intelligence_mutex);
                std::cout << "🧠 Thread " << thread_id << " tested " << local_tested 
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

// 🧠 PATTERN EXPLORATION WORKER - TESTS MATHEMATICAL STRUCTURES
void pattern_exploration_worker() {
    std::cout << "🔬 Pattern Exploration Worker started...\n";
    
    std::random_device rd;
    std::mt19937_64 gen(rd());
    
    EC_KEY* kk = EC_KEY_new_by_curve_name(NID_secp256k1);
    BIGNUM* test_bn = BN_new();
    
    // Test various mathematical structures
    std::vector<BIGNUM*> test_patterns;
    
    // Powers of 2 and nearby values
    for (int exp = 70; exp < 72; exp++) {
        BIGNUM* pattern = BN_new();
        BN_set_bit(pattern, exp);
        test_patterns.push_back(pattern);
        
        // Nearby values
        for (int offset = -1000; offset <= 1000; offset += 100) {
            BIGNUM* nearby = BN_dup(pattern);
            if (offset > 0) {
                BN_add_word(nearby, offset);
            } else {
                BN_sub_word(nearby, -offset);
            }
            if (BN_cmp(nearby, RANGE_START_BN) >= 0 && BN_cmp(nearby, RANGE_END_BN) <= 0) {
                test_patterns.push_back(nearby);
            } else {
                BN_free(nearby);
            }
        }
    }
    
    // Arithmetic progressions from promising areas
    auto promising_centers = []() -> std::vector<BIGNUM*> {
        // This would come from global intelligence in a real implementation
        std::vector<BIGNUM*> centers;
        BIGNUM* center = BN_new();
        BN_set_bit(center, 70);
        BN_add_word(center, 123456789);
        centers.push_back(center);
        return centers;
    }();
    
    for (BIGNUM* center : promising_centers) {
        for (int step = 1; step <= 10000; step *= 10) {
            for (int dir = -1; dir <= 1; dir += 2) {
                BIGNUM* progression = BN_dup(center);
                BN_add_word(progression, dir * step);
                if (BN_cmp(progression, RANGE_START_BN) >= 0 && BN_cmp(progression, RANGE_END_BN) <= 0) {
                    test_patterns.push_back(progression);
                } else {
                    BN_free(progression);
                }
            }
        }
    }
    
    // Test all patterns
    for (BIGNUM* pattern : test_patterns) {
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
        
        if (match_length >= 10) {
            GLOBAL_INTELLIGENCE.analyze_key_pattern(pattern, match_length, hash);
            std::string pattern_hex = BN_bn2hex(pattern);
            std::cout << "🔬 Pattern match: " << match_length 
                      << " chars at " << pattern_hex.substr(0, 20) << "...\n";
        }
        
        if (hash == TARGET) {
            FOUND = true;
            std::string found_key_hex = BN_bn2hex(pattern);
            std::cout << "\n🎯 PATTERN EXPLORATION SOLVED THE PUZZLE!\n";
            std::cout << "Key: " << found_key_hex << "\n";
        }
        
        TOTAL_TRIED++;
    }
    
    // Cleanup
    for (BIGNUM* pattern : test_patterns) {
        BN_free(pattern);
    }
    BN_free(test_bn);
    EC_KEY_free(kk);
}

// 🧠 MAIN INTELLIGENT SEARCH SYSTEM
void real_time_intelligent_search() {
    std::cout << "🧠 LAUNCHING REAL-TIME INTELLIGENT FILTERING SYSTEM 🧠\n";
    std::cout << "Target: " << TARGET << "\n";
    std::cout << "Range: 2^70 to 2^71\n";
    std::cout << "Threads: " << N_THREADS << "\n";
    std::cout << "This system does CONTINUOUS PATTERN ANALYSIS and ADAPTIVE FILTERING!\n\n";
    
    auto start_time = std::chrono::steady_clock::now();
    
    std::vector<std::thread> threads;
    
    // Start intelligent workers
    std::cout << "🚀 Starting " << N_THREADS << " intelligent workers...\n";
    for (int i = 0; i < N_THREADS; i++) {
        threads.emplace_back(intelligent_worker, i);
    }
    
    // Start pattern exploration
    std::cout << "🔬 Starting pattern exploration...\n";
    threads.emplace_back(pattern_exploration_worker);
    
    // Intelligence monitoring
    uint64_t last_total = 0;
    auto last_intelligence_report = std::chrono::steady_clock::now();
    
    while (!FOUND) {
        std::this_thread::sleep_for(std::chrono::seconds(30));
        
        auto now = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - start_time);
        uint64_t current_total = TOTAL_TRIED;
        uint64_t keys_per_sec = (current_total - last_total) / 30;
        
        std::cout << "\n📊 REAL-TIME PROGRESS:\n";
        std::cout << "⏱️  Elapsed: " << elapsed.count() << "s | ";
        std::cout << "Rate: " << keys_per_sec << " keys/s | ";
        std::cout << "Total: " << current_total << " keys\n";
        
        // Regular intelligence reports
        if (now - last_intelligence_report > std::chrono::minutes(2)) {
            GLOBAL_INTELLIGENCE.print_intelligence_report();
            last_intelligence_report = now;
        }
        
        last_total = current_total;
        
        // Adaptive timeout
        if (elapsed.count() > 3600 && keys_per_sec < 1000) {
            std::cout << "💤 Low detection rate - the system is focusing on promising areas.\n";
        }
    }
    
    // Wait for all threads
    for (auto& t : threads) {
        if (t.joinable()) t.join();
    }
    
    auto end_time = std::chrono::steady_clock::now();
    auto total_elapsed = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    
    std::cout << "\n🎉 REAL-TIME INTELLIGENT SYSTEM COMPLETED!\n";
    std::cout << "Time: " << total_elapsed.count() << " seconds\n";
    std::cout << "Total Keys Tested: " << TOTAL_TRIED << "\n";
    std::cout << "Final Intelligence Report:\n";
    GLOBAL_INTELLIGENCE.print_intelligence_report();
}

int main() {
    std::cout << "🧠 REAL-TIME INTELLIGENT BITCOIN PUZZLE SOLVER 🧠\n";
    std::cout << "=================================================\n\n";
    
    // Initialize
    init_bignum_constants();
    
    auto start = std::chrono::steady_clock::now();
    real_time_intelligent_search();
    auto end = std::chrono::steady_clock::now();
    
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    std::cout << "\n=================================================\n";
    std::cout << "Total execution time: " << elapsed.count() << " seconds\n";
    
    if (FOUND) {
        std::cout << "🏆 REAL-TIME INTELLIGENCE TRIUMPHS! 🏆\n";
    } else {
        std::cout << "🤖 System continues to learn and adapt... 🤖\n";
    }
    
    cleanup_bignum_constants();
    return 0;
}