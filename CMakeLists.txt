cmake_minimum_required(VERSION 3.16)
project(SAFE_Benchmark
    VERSION 1.0
    DESCRIPTION "Ultra-Adaptive Evolutionary Bitcoin Benchmark (secp256k1 + AVX2 + OpenSSL)"
    LANGUAGES CXX)

# ---- Compiler configuration ----
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Optimize for native CPU
add_compile_options(-O3 -march=native -flto -mavx2 -mbmi -pthread -Wno-deprecated-declarations)
add_link_options(-flto)

# ---- Create the main executable ----
add_executable(safe_benchmark main.cpp)

# ---- Threads ----
find_package(Threads REQUIRED)
target_link_libraries(safe_benchmark PRIVATE Threads::Threads)

# ---- OpenSSL ----
find_package(OpenSSL REQUIRED)
if (OpenSSL_FOUND)
    message(STATUS "Found OpenSSL: ${OPENSSL_VERSION}")
    target_include_directories(safe_benchmark PRIVATE ${OPENSSL_INCLUDE_DIR})
    target_link_libraries(safe_benchmark PRIVATE OpenSSL::SSL OpenSSL::Crypto)
else()
    message(FATAL_ERROR "❌ OpenSSL not found. Please install libssl-dev or openssl-devel.")
endif()

# ---- Try to find existing libsecp256k1 ----
find_path(SECP256K1_INCLUDE_DIR secp256k1.h)
find_library(SECP256K1_LIBRARY secp256k1)

if (SECP256K1_INCLUDE_DIR AND SECP256K1_LIBRARY)
    message(STATUS "✅ Found system libsecp256k1: ${SECP256K1_LIBRARY}")
    target_include_directories(safe_benchmark PRIVATE ${SECP256K1_INCLUDE_DIR})
    target_link_libraries(safe_benchmark PRIVATE ${SECP256K1_LIBRARY})
else()
    message(WARNING "⚠️ libsecp256k1 not found — it will be built automatically from GitHub")

    include(FetchContent)
    FetchContent_Declare(
        secp256k1
        GIT_REPOSITORY https://github.com/bitcoin-core/secp256k1.git
        GIT_TAG master
    )

    # Configure minimal build
    FetchContent_GetProperties(secp256k1)
    if (NOT secp256k1_POPULATED)
        FetchContent_Populate(secp256k1)
        set(SECP256K1_BUILD_TESTS OFF CACHE BOOL "" FORCE)
        set(SECP256K1_BUILD_BENCHMARK OFF CACHE BOOL "" FORCE)
        set(SECP256K1_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
        set(SECP256K1_ENABLE_MODULE_ECDH ON CACHE BOOL "" FORCE)
        set(SECP256K1_ENABLE_MODULE_RECOVERY ON CACHE BOOL "" FORCE)

        add_subdirectory(${secp256k1_SOURCE_DIR} ${secp256k1_BINARY_DIR} EXCLUDE_FROM_ALL)
    endif()

    target_link_libraries(safe_benchmark PRIVATE secp256k1)
    target_include_directories(safe_benchmark PRIVATE ${secp256k1_SOURCE_DIR}/include)
endif()

# ---- Platform-specific (macOS OpenSSL from Homebrew) ----
if(APPLE)
    include_directories(/usr/local/opt/openssl/include /opt/homebrew/opt/openssl/include)
    link_directories(/usr/local/opt/openssl/lib /opt/homebrew/opt/openssl/lib)
endif()

# ---- Output directory ----
set_target_properties(safe_benchmark PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "----------------------------------------")
message(STATUS "✅ Configuration complete!")
message(STATUS "Run:")
message(STATUS "   cd build && make -j$(nproc)")
message(STATUS "Binary will be in: ${CMAKE_BINARY_DIR}/bin/safe_benchmark")
