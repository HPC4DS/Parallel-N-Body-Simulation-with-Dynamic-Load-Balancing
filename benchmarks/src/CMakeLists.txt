add_library(benchmark benchmark.c)

# Public headers
target_include_directories(benchmark PUBLIC ${CMAKE_SOURCE_DIR}/benchmarks/include/)

find_and_link_mpi_for_target(benchmark PRIVATE)
find_and_link_openmp_for_target(benchmark PRIVATE)

target_link_libraries(benchmark PRIVATE utils)
target_link_libraries(benchmark PRIVATE print_debug)

target_include_directories(benchmark PUBLIC ${CMAKE_SOURCE_DIR}/include)

set(BENCHMARK_PRINT_DEBUG_LEVEL ${GLOBAL_PRINT_DEBUG_LEVEL} CACHE STRING "Debug level for benchmark")
define_print_debug_level_for_target(benchmark ${BENCHMARK_PRINT_DEBUG_LEVEL})