//
// Created by Matteo Ranzi on 02/11/25.
//

#ifndef CLION_OPENMP_MPI_BENCHMARK_H
#define CLION_OPENMP_MPI_BENCHMARK_H

#include <mpi.h>
#include <stdio.h>
#include <string.h>

//TODO currently local benchmark is overwritten with previous ones

#ifdef __cplusplus
extern "C" {
#endif

    typedef void (*Application)(int rank);

    typedef struct {
        char name[256];
        int n_iterations;
        const char *logs_dir;
        MPI_File *mpi_log_file;
    } BenchmarkConfig;

    typedef struct {
        double min_time;
        double max_time;
        double avg_time;
    } BenchmarkResult;

    void benchmark_init(int my_rank, const BenchmarkConfig* benchmark_config);
    void benchmark_finalize(int my_rank, const BenchmarkConfig* benchmark_config);

    void benchmark_run(int my_rank, const BenchmarkConfig* benchmark_config, Application app, BenchmarkResult* benchmark_result);

#ifdef __cplusplus
}
#endif

#endif //CLION_OPENMP_MPI_BENCHMARK_H