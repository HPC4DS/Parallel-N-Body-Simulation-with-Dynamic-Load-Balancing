//
// Created by Matteo Ranzi on 08/12/25.
//

//TODO: examine the real usage of memory in GB


#include <cstring>
#include <mpi.h>

#include <unistd.h>
#include <vector>

#include "debug/print_debug.h"
#include "utils/build_info.h"
#include "utils/random_values.hpp"

#include "benchmark.hpp"
#include "linear_algorithms.hpp"

constexpr int defaultRadixBits = 10;
constexpr int min_vector_size_exp = 10;
constexpr int max_vector_size_exp = 32;

static size_t parse_size_safe(const char* s);

int main(int argc, char *argv[]) {
    int comm_size;
    int my_rank;
    MPI_File benchmark_log_file;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    print_build_info(my_rank);

    const size_t radix_bits = parse_size_safe(getenv("PBS_ARRAY_INDEX")) == static_cast<size_t>(-1) ?
                               defaultRadixBits : parse_size_safe(getenv("PBS_ARRAY_INDEX"));

    PRINT_DEBUG_INFO_R(my_rank, "PBS array index: %zu\n", radix_bits);

    //=============================================================================
    BenchmarkConfig benchmark_config;
    default_benchmark_config(&benchmark_config);
    benchmark_config.repetitions = 15;
    benchmark_config.mpi_log_file = &benchmark_log_file;
    snprintf(benchmark_config.description, sizeof(benchmark_config.description), "linear radix_sort<%zu> (<0> is std::sort)", radix_bits);
    std::strcpy(benchmark_config.sweep_name, "vector_size");
    std::strcpy(benchmark_config.const_name, "radix_bits");
    benchmark_config.const_value = static_cast<size_t>(radix_bits);

    benchmark_init(my_rank, &benchmark_config);


    // Use size_t for sizes to avoid overflow and match random_fill signature
    size_t vector_size;
    std::vector<uint64_t> numbers;
    std::function<void()> pre = [&]() {
        try {
            random_fill(numbers, vector_size, 0LU, 1ULL << 63);
        } catch (const std::length_error &e) {
            PRINT_DEBUG_ERROR_R(my_rank, "Allocation failed in random_fill: %s", e.what());
            // abort the MPI job to avoid continuing with invalid state
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    };
    std::function<void()> app;
    switch (radix_bits) {
        case 0: app = [&]{std::sort(numbers.begin(), numbers.end());}; break;
        case 1: app = [&] {radix_sort<1>(numbers.begin(), numbers.end());}; break;
        case 2: app = [&] {radix_sort<2>(numbers.begin(), numbers.end());}; break;
        case 3: app = [&] {radix_sort<3>(numbers.begin(), numbers.end());}; break;
        case 4: app = [&] {radix_sort<4>(numbers.begin(), numbers.end());}; break;
        case 5: app = [&] {radix_sort<5>(numbers.begin(), numbers.end());}; break;
        case 6: app = [&] {radix_sort<6>(numbers.begin(), numbers.end());}; break;
        case 7: app = [&] {radix_sort<7>(numbers.begin(), numbers.end());}; break;
        case 8: app = [&] {radix_sort<8>(numbers.begin(), numbers.end());}; break;
        case 9: app = [&] {radix_sort<9>(numbers.begin(), numbers.end());}; break;
        case 10: app = [&] {radix_sort<10>(numbers.begin(), numbers.end());}; break;
        case 11: app = [&] {radix_sort<11>(numbers.begin(), numbers.end());}; break;
        case 12: app = [&] {radix_sort<12>(numbers.begin(), numbers.end());}; break;
        case 13: app = [&] {radix_sort<13>(numbers.begin(), numbers.end());}; break;
        case 14: app = [&] {radix_sort<14>(numbers.begin(), numbers.end());}; break;
        case 15: app = [&] {radix_sort<15>(numbers.begin(), numbers.end());}; break;
        case 16: app = [&] {radix_sort<16>(numbers.begin(), numbers.end());}; break;
        default: break;
    }
    std::function<void()> post = [&]() {
        // Verify sorting
        for (size_t i = 1; i < numbers.size(); ++i) {
            if (numbers[i - 1] > numbers[i]) {
                PRINT_DEBUG_ERROR_R(my_rank, "SORTING FAILED");
                MPI_Abort(comm_size, MPI_COMM_WORLD);
            }
        }
    };

    // Build sweep values carefully: use size_t and avoid shifts that overflow or exceed limits
    std::vector<size_t> sweep_values;
    for (int i = min_vector_size_exp; i <= max_vector_size_exp; i++) {
        const size_t candidate = (1ULL << i);
        // stop if candidate exceeds what a vector can allocate
        if (candidate == 0) break;
        if (candidate > numbers.max_size()) break;
        sweep_values.push_back(candidate);
    }

    // Run benchmarks for each sweep value
    for (const auto& sv : sweep_values) {
        vector_size = sv;
        benchmark_config.sweep_value = sv;
        benchmark_run(my_rank, REPETITION_STRATEGY::MIN_TIME, &benchmark_config, pre, app, post, nullptr);
    }

    benchmark_finalize(my_rank, &benchmark_config);
    //=============================================================================


    MPI_Finalize();

    return 0;
}

static size_t parse_size_safe(const char* s) {
    if (!s) return -1;
    try {
        const unsigned long long v = std::stoull(std::string(s));
        return v;
    } catch (...) {
        return -1;
    }
}