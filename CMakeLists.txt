cmake_minimum_required(VERSION 3.10)
project(IntelligentPuzzleSolver)

set(CMAKE_CXX_STANDARD 17)

# OPTIMIZATION FOR PERFORMANCE
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -flto -fomit-frame-pointer -funroll-loops -pthread")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -O3 -flto -pthread")

find_package(OpenSSL REQUIRED)
add_executable(intelligent_solver main.cpp)

target_link_libraries(intelligent_solver OpenSSL::SSL OpenSSL::Crypto)

# SUPPRESS DEPRECATION WARNINGS FOR EDUCATIONAL CODE
target_compile_options(intelligent_solver PRIVATE -w)