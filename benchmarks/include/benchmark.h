//
// Created by Matteo Ranzi on 02/11/25.
//

#ifndef CLION_OPENMP_MPI_BENCHMARK_H
#define CLION_OPENMP_MPI_BENCHMARK_H

#include <mpi.h>
#include <stdio.h>

//TODO currently local benchmark is overwritten with previous ones
//TODO convert benchmark to C++

#ifdef HPC_RUN
#define BENCHMARK_DIR getenv("HPC_BENCHMARK_DIR")
#define APPLICATION_ID getenv("JOB_ID")
#else
#define BENCHMARK_DIR getenv("LOCAL_BENCHMARK_DIR")
#define APPLICATION_ID getenv("APPLICATION_ID")
#endif

#ifdef __cplusplus
extern "C" {
#endif

    typedef void (*Application)(void *arguments);

    typedef struct {
        char description[256];
        int max_iterations; // Maximum number of iterations for the benchmark
        int repetitions; // How many times to repeat the entire benchmark
        double min_time;
        const char *logs_dir;
        MPI_File *mpi_log_file;

        size_t const_value;
        char const_name[64];

        size_t sweep_value;
        char sweep_name[64];
    } BenchmarkConfig;

    typedef struct {
        double min_time;
        double max_time;
        double avg_time;
    } BenchmarkResult;

    typedef enum REPETITION_STRATEGY {
        FIXED_ITERATIONS,
        MIN_TIME
    } REPETITION_STRATEGY;

    void default_benchmark_config(BenchmarkConfig* benchmark_config);
    void benchmark_init(int my_rank, const BenchmarkConfig* benchmark_config);
    void benchmark_finalize(int my_rank, const BenchmarkConfig* benchmark_config);

    void benchmark_run(int my_rank, REPETITION_STRATEGY repetition_strategy, const BenchmarkConfig* benchmark_config, Application preHook, void *preHookArgs, Application app, void *appArgs, Application postHook, void *postHookArgs, BenchmarkResult* benchmark_result);

#ifdef __cplusplus
}
#endif

#endif //CLION_OPENMP_MPI_BENCHMARK_H