//
// Created by Matteo Ranzi on 20/11/25.
//

#ifndef CLION_OPENMP_MPI_BENCHMARK_HPP
#define CLION_OPENMP_MPI_BENCHMARK_HPP

#include "benchmark.h"
#include <functional>

extern "C" inline void benchmark_cpp_trampoline(void* userdata) {
    const auto* fn = static_cast<std::function<void()>*>(userdata);
    (*fn)();
}

inline void benchmark_run(const int my_rank, const BenchmarkConfig* benchmark_config, std::function<void()>& fn, BenchmarkResult* benchmark_result) {
    return benchmark_run_c(my_rank,benchmark_config, benchmark_cpp_trampoline, &fn, benchmark_result);
}

template <typename Lambda>
void benchmark_run(const int my_rank, const BenchmarkConfig* benchmark_config, Lambda&& lambda, BenchmarkResult* benchmark_result) {
    std::function<void()> fn = lambda;
    return benchmark_run_c(my_rank, benchmark_config, benchmark_cpp_trampoline, &fn, benchmark_result);
}

#endif //CLION_OPENMP_MPI_BENCHMARK_HPP