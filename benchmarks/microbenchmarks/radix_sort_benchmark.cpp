//
// Created by Matteo Ranzi on 08/12/25.
//

#include <mpi.h>

#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <vector>

#include "debug/print_debug.h"
#include "utils/build_info.h"
#include "utils/random_values.hpp"

#include "benchmark.hpp"
#include "linear_algorithms.hpp"

constexpr int radixBits = 10;
constexpr int max_vector_size_exp = 10; // up to 2^32 elements

int main(int argc, char *argv[]) {
    int comm_size;
    int my_rank;
    MPI_File benchmark_log_file;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    print_build_info(my_rank);


    //=============================================================================
    BenchmarkConfig benchmark_config;
    default_benchmark_config(&benchmark_config);
    benchmark_config.mpi_log_file = &benchmark_log_file;
    snprintf(benchmark_config.description, sizeof(benchmark_config.description), "linear radix_sort<%d>", radixBits);
    strcpy(benchmark_config.sweep_name, "vector_size");

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
    std::function<void()> app = [&]() {
        static_assert(radixBits >= 1 && radixBits <= 16, "radixBits must be between 1 and 16");
        // direct compile-time call removes unreachable-case warnings
        radix_sort<radixBits>(numbers.begin(), numbers.end());
    };
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
    for (int i = 0; i <= max_vector_size_exp; i++) {
        const size_t candidate = (1ULL << i);
        // stop if candidate exceeds what a vector can allocate
        if (candidate == 0) break;
        if (candidate > numbers.max_size()) break;
        sweep_values.push_back(candidate);
    }

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
