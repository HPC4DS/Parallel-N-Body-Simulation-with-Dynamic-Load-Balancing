//
// Created by Matteo Ranzi on 08/12/25.
//

#include <omp.h>
#include <mpi.h>

#include <cstdlib>
#include <cstring>
#include <unistd.h>

#include "debug/print_debug.h"
#include "debug/unique_print_debug.h"
#include "utils/build_info.h"
#include "utils/random_values.hpp"

#include "benchmark.hpp"
#include "BuildInfo.h"
#include "linear_algorithms.hpp"

constexpr int radixBits = 10;
constexpr int max_vector_size_exp = 32;

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

    int vector_size;
    std::vector<uint64_t> numbers;
    std::function<void()> pre = [&]() {
        random_fill(numbers, vector_size, 0LU, 1LU << 63);
    };
    std::function<void()> app = [&]() {
        switch (radixBits) {
        case 1:  radix_sort<1>(numbers.begin(), numbers.end()); break;
        case 2:  radix_sort<2>(numbers.begin(), numbers.end()); break;
        case 3:  radix_sort<3>(numbers.begin(), numbers.end()); break;
        case 4:  radix_sort<4>(numbers.begin(), numbers.end()); break;
        case 5:  radix_sort<5>(numbers.begin(), numbers.end()); break;
        case 6:  radix_sort<6>(numbers.begin(), numbers.end()); break;
        case 7:  radix_sort<7>(numbers.begin(), numbers.end()); break;
        case 8:  radix_sort<8>(numbers.begin(), numbers.end()); break;
        case 9:  radix_sort<9>(numbers.begin(), numbers.end()); break;
        case 10: radix_sort<10>(numbers.begin(), numbers.end()); break;
        case 11: radix_sort<11>(numbers.begin(), numbers.end()); break;
        case 12: radix_sort<12>(numbers.begin(), numbers.end()); break;
        case 13: radix_sort<13>(numbers.begin(), numbers.end()); break;
        case 14: radix_sort<14>(numbers.begin(), numbers.end()); break;
        case 15: radix_sort<15>(numbers.begin(), numbers.end()); break;
        default: break;
        }
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

    std::vector<int> sweep_values;
    for (int i = 0; i <= max_vector_size_exp; i++) {
        sweep_values.push_back(1<<i);
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
