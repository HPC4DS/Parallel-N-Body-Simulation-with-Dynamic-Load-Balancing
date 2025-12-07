/***********************************************************************************************************************
* Project: Parallel N-Body Simulation with Dynamic Load Balancing
* Author : Matteo Ranzi
* Description: Implement a parallel gravitational N-body simulation that runs efficiently on a cluster using MPI for distributed domain decomposition and OpenMP for intra-node parallelism.
***********************************************************************************************************************/

#include <omp.h>
#include <mpi.h>

#include <cstdlib>
#include <unistd.h>
#include <cstring>
#include <unistd.h>

#include "debug/print_debug.h"
#include "debug/unique_print_debug.h"
#include "utils/info_utils.h"

#include "benchmark.hpp"
#include "BuildInfo.h"

#ifdef HPC_RUN
#define BENCHMARK_DIR getenv("HPC_BENCHMARK_DIR")
#else
#define BENCHMARK_DIR getenv("LOCAL_BENCHMARK_DIR")
#endif

#define APPLICATION_ID getenv("APPLICATION_ID")


void wait_for_attach_debugger(int rank);

int main(int argc, char *argv[]) {
    int comm_size;
    int my_rank;
    MPI_File benchmark_log_file;


    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    print_build_info(my_rank);

    wait_for_attach_debugger(my_rank);

    //=============================================================================
    BenchmarkConfig benchmark_config;
    const std::string logs_dir_str = std::string(BENCHMARK_DIR ? BENCHMARK_DIR : "") + "/" + (APPLICATION_ID ? APPLICATION_ID : "default_app");
    benchmark_config.logs_dir = strdup(logs_dir_str.c_str());
    benchmark_config.mpi_log_file = &benchmark_log_file;
    strcpy(benchmark_config.description, "Main app testing benchmark");
    benchmark_config.max_iterations = 50;
    benchmark_config.min_time = 2.0; // seconds

    benchmark_init(my_rank, &benchmark_config);
    std::function<void()> pre = [&](){};
    std::function<void()> app = [&](){PRINT_DEBUG_INFO_R(my_rank, "Hello Parallel N-Body Simulation with Dynamic Load Balancing\n");};
    std::function<void()> post = [&](){};
    benchmark_run(my_rank, REPETITION_STRATEGY::MIN_TIME, &benchmark_config, pre, app, post, nullptr);
    benchmark_finalize(my_rank, &benchmark_config);
    //=============================================================================

    
    MPI_Finalize();

    return 0;
}


/**
 * <a href="https://www.jetbrains.com/help/clion/openmpi.html#debug-procedure">
 *      How to debug OpenMPI projects in CLion
 * </a>
 *
 * @param rank rank of the process.
 * @brief If environment variable WAIT_FOR_DEBUGGER is set to 1, Process 0 blocks until 'breakpoint' value is changed from debugger.<br/>
 *        All the other MPI processes (if exist) wait via MPI_Barrier.
 */
void wait_for_attach_debugger(const int rank) {
    const char* flag = std::getenv("WAIT_FOR_DEBUGGER");
    if (flag && strcmp(flag, "1") == 0) {
        if (rank == 0) {
            printf("\n*** Waiting for debugger to attach to process %d ***\n", getpid());
            // ReSharper disable once CppLocalVariableMayBeConst
            int breakpoint = 0;
            // ReSharper disable once CppDFALoopConditionNotUpdated // ReSharper disable once CppDFAConstantConditions // ReSharper disable once CppDFAEndlessLoop
            while (!breakpoint) {
                sleep(3); // >>> A breakpoint must be here <<<
            }
        }

        //If MPI execution, all processes wait for Process 0
#ifdef MPI_VERSION
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }
}

// =================================================================
//  Generated with Matteo Ranzi’s CLion Project Template
//  Kickstart C/C++ projects with MPI, OpenMP and HPC@UNITN remote deployment
// =================================================================