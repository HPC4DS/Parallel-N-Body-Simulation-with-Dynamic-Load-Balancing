//
// Created by Matteo Ranzi on 02/11/25.
//

#ifndef CLION_OPENMP_MPI_BENCHMARK_H
#define CLION_OPENMP_MPI_BENCHMARK_H

#include <mpi.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

    typedef void (*Application)(int rank);

    /*typedef struct {
        const char* name;
        double avg_time;
        size_t problem_size;
    } BenchmarkResult;*/

    void benchmark_init(int my_rank, const char *logs_dir, MPI_File *benchmark_log_file);
    void benchmark_finalize(int my_rank, MPI_File *benchmark_log_file);

    void benchmark_run(int my_rank, const char* benchmark_name, Application app, const MPI_File *benchmark_log_file);

#ifdef __cplusplus
}
#endif

#endif //CLION_OPENMP_MPI_BENCHMARK_H