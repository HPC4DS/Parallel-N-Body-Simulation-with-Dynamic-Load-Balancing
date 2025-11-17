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
#include <iostream>

#include "debug/print_debug.h"
#include "debug/unique_print_debug.h"

#include "benchmark/benchmark.h"
#include "BuildInfo.h"

#ifdef HPC_RUN
#define LOGS_DIR getenv("HPC_JOB_LOGS_DIR")
#else
#define LOGS_DIR getenv("LOCAL_LOGS_DIR")
#endif


void wait_for_attach_debugger(int rank);
void print_build_info(int rank);

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
    benchmark_init(my_rank, LOGS_DIR, &benchmark_log_file);
    benchmark_run(my_rank, "Void Application", [](const int rank){PRINT_DEBUG_INFO(rank, "Hello Parallel N-Body Simulation with Dynamic Load Balancing\n");}, &benchmark_log_file);
    benchmark_finalize(my_rank, &benchmark_log_file);
    //=============================================================================

    
    MPI_Finalize();

    return 0;
}

void print_build_info(int rank) {
    UNIQUE_PRINT(rank, "Build Version: %s\n", BUILD_VERSION);
    UNIQUE_PRINT(rank, "Git Commit Hash: %s\n", GIT_COMMIT_HASH);
    UNIQUE_PRINT(rank, "Git Dirty: %s\n", GIT_DIRTY);
    UNIQUE_PRINT(rank, "Build Type: %s\n", BUILD_TYPE);
    UNIQUE_PRINT(rank, "Compiler: %s\n", CXX_COMPILER);
    UNIQUE_PRINT(rank, "CXX Flags: %s\n", CXX_FLAGS);
    UNIQUE_PRINT(rank, "Build Timestamp: %s\n", BUILD_TIMESTAMP);
    UNIQUE_PRINT(rank, "System: %s (%s)\n", SYSTEM_NAME, SYSTEM_PROCESSOR);
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