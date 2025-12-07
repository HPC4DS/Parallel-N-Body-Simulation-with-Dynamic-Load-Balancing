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

inline void benchmark_run(const int my_rank, const BenchmarkConfig* benchmark_config, std::function<void()>& preHook, std::function<void()>& app, std::function<void()>& postHook, BenchmarkResult* benchmark_result) {
    return benchmark_run_c(my_rank,benchmark_config, benchmark_cpp_trampoline, &preHook, benchmark_cpp_trampoline, &app, benchmark_cpp_trampoline, &postHook, benchmark_result);
}

template <typename Lambda>
void benchmark_run(const int my_rank, const BenchmarkConfig* benchmark_config, Lambda&& preHook, Lambda&& app, Lambda&& postHook, BenchmarkResult* benchmark_result) {
    std::function<void()> fn_app = app;
    std::function<void()> fn_preHook = preHook;
    std::function<void()> fn_postHook = postHook;
    return benchmark_run_c(my_rank, benchmark_config, benchmark_cpp_trampoline, &fn_preHook, benchmark_cpp_trampoline, &fn_app, benchmark_cpp_trampoline, &fn_postHook, benchmark_result);
}

#endif //CLION_OPENMP_MPI_BENCHMARK_HPP