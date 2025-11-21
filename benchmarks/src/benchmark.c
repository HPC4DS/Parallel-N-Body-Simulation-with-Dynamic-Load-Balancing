//
// Created by Matteo Ranzi on 02/11/25.
//

#include "benchmark.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <stdarg.h>
#include <unistd.h>
#include <errno.h>

#include <float.h>
#include "debug/print_debug.h"
#include "debug/unique_print_debug.h"
#include "utils/info_utils.h"

#define MAX_LOG_FILENAME_LENGTH 1024
#define BENCHMARK_FILE_NAME "benchmark.log"

/**
 * Aborts the MPI program with the given error code.
 *
 * @param err The error code to exit with.
 */
static void benchmark_abort(const int err){
    int mpi_initialized = 0;
    MPI_Initialized(&mpi_initialized);
    if (mpi_initialized) {
        /* MPI_Abort does not return; use a non-zero error code. */
        MPI_Abort(MPI_COMM_WORLD, err == 0 ? 1 : err);
    } else {
        exit(err == 0 ? 1 : err);
    }
}

/**
 * Handles MPI errors by printing an error message and exiting the program.
 *
 * @param my_rank The rank of the current process.
 * @param ret_value The MPI error code.
 * @param message A custom error message to display.
 */
static void handle_mpi_error_and_abort(const int my_rank, const int ret_value, const char *message) {
    char error_string[MPI_MAX_ERROR_STRING];
    int error_string_len;
    MPI_Error_string(ret_value, error_string, &error_string_len);
    PRINT_DEBUG_ERROR_R(my_rank, "%s: %s\n", message, error_string);

    benchmark_abort(ret_value);
}

/**
 * Creates the benchmark log filename by combining the logs directory and the benchmark file name.
 *
 * @param my_rank The rank of the current process.
 * @param benchmark_log_filename The buffer to store the generated filename.
 * @param logs_dir The directory where the log file will be stored.
 */
static void create_benchmark_filename(const int my_rank, char *benchmark_log_filename, const char *logs_dir) {
    const int len = snprintf(benchmark_log_filename, MAX_LOG_FILENAME_LENGTH, "%s/%s", logs_dir, BENCHMARK_FILE_NAME);
    if (len < 0 || len >= MAX_LOG_FILENAME_LENGTH) {
        PRINT_DEBUG_ERROR_R(my_rank, "[BENCHMARK] Benchmark log filename too long (truncated)\n");
        benchmark_abort(1);
    }
}

/**
 * Opens the benchmark log file for writing. If the file already exists, it is removed.
 *
 * @param my_rank The rank of the current process.
 * @param filename The name of the file to open.
 * @param mpi_file A pointer to the MPI file handle.
 */
static void open_benchmark_file(const int my_rank, const char *filename, MPI_File *mpi_file) {
    if (my_rank == 0) {
        if (remove(filename) != 0 && errno != ENOENT) {
            PRINT_DEBUG_ERROR_R(my_rank, "[BENCHMARK] Failed to remove existing file '%s': %s\n", filename, strerror(errno));
            benchmark_abort(1);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    const int ret = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_APPEND, MPI_INFO_NULL, mpi_file);
    if (ret != MPI_SUCCESS) {
        handle_mpi_error_and_abort(my_rank, ret, "[BENCHMARK] Cannot open benchmark log file");
    }
}

/**
 * Logs a unique message to the benchmark file from rank 0.
 *
 * @param my_rank The rank of the current process.
 * @param file A pointer to the MPI file handle.
 * @param format The format string for the log message.
 * @param ... Additional arguments for the format string.
 */
static void benchmark_unique_log(const int my_rank, const MPI_File *file, const char *format, ...) {
    if (my_rank != 0) return;

    char buffer[4096];
    va_list args;
    va_start(args, format);
    int len = vsnprintf(buffer, sizeof(buffer), format, args);
    va_end(args);

    if (len < 0) return;

    if (len < sizeof(buffer)) {
        if (MPI_File_write(*file, buffer, len, MPI_CHAR, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
            handle_mpi_error_and_abort(my_rank, MPI_ERR_OTHER, "[BENCHMARK] Error writing to benchmark file");
        }
    } else {
        PRINT_DEBUG_WARN_R(my_rank, "[BENCHMARK] Log message too long (truncated)\n");
    }
}

/**
 * Writes the benchmark log header, including the start time, benchmark name, and processor names.
 *
 * @param my_rank The rank of the current process.
 * @param file A pointer to the MPI file handle.
 * @param benchmark_name The name of the benchmark.
 */
static void write_benchmark_log_header(const int my_rank, const MPI_File *file, const char *benchmark_name) {
    // Print date and time
    if (my_rank == 0) {
        const time_t t = time(NULL);
        struct tm tm;
        localtime_r(&t, &tm);

        benchmark_unique_log(my_rank, file, "==================================================\n");
        benchmark_unique_log(my_rank, file, "Benchmark started at: %d-%02d-%02d %02d:%02d:%02d\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
        benchmark_unique_log(my_rank, file, "Benchmark Name: %s\n", benchmark_name);
    }

    // Gather fixed-size processor names to rank 0
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    char node_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(node_name, &name_len);

    char *all_names = NULL;
    if (my_rank == 0) {
        all_names = malloc((size_t)world_size * MPI_MAX_PROCESSOR_NAME);
        if (!all_names) {
            handle_mpi_error_and_abort(my_rank, MPI_ERR_OTHER, "[BENCHMARK] malloc failed");
        }
    }

    const int ret = MPI_Gather(node_name, MPI_MAX_PROCESSOR_NAME, MPI_CHAR,
               all_names, MPI_MAX_PROCESSOR_NAME, MPI_CHAR,
               0, MPI_COMM_WORLD);

    if (ret != MPI_SUCCESS) {
        handle_mpi_error_and_abort(my_rank, ret, "[BENCHMARK] MPI_Gather failed");
    }

    benchmark_unique_log(my_rank, file, "Number of MPI processes: %d\n", world_size);
    if (my_rank == 0) {
        for (int r = 0; r < world_size; r++) {
            char *recv_name = all_names + r * MPI_MAX_PROCESSOR_NAME;
            benchmark_unique_log(my_rank, file, "Rank: %d > node: %s\n", r, recv_name);
        }
        free(all_names);
    }

    benchmark_unique_log(my_rank, file, "------------------\n");


    char build_info[512];
    build_info_to_string(build_info, sizeof(build_info));
    benchmark_unique_log(my_rank, file, "Build Information:\n");
    benchmark_unique_log(my_rank, file, "%s", build_info);
    benchmark_unique_log(my_rank, file, "==================================================\n\n");
}

/**
 * Initializes the benchmark by creating and opening the log file.
 * benchmark_init is a collective operation
 *
 * @param my_rank The rank of the current process.
 * @param benchmark_config The configuration for the benchmark.
 */
void benchmark_init(const int my_rank, const BenchmarkConfig* benchmark_config) {
    char benchmark_log_filename[MAX_LOG_FILENAME_LENGTH];

    create_benchmark_filename(my_rank, benchmark_log_filename, benchmark_config->logs_dir);
    open_benchmark_file(my_rank, benchmark_log_filename, benchmark_config->mpi_log_file);

    UNIQUE_PRINT_DEBUG_INFO(my_rank, "[BENCHMARK] Output log to: %s\n", benchmark_log_filename);

    write_benchmark_log_header(my_rank, benchmark_config->mpi_log_file, benchmark_config->description);
}

/**
 * Finalizes the benchmark by closing the log file.
 * benchmark_finalize is a collective operation
 *
 * @param my_rank The rank of the current process.
 * @param benchmark_config The configuration for the benchmark.
 */
void benchmark_finalize(const int my_rank, const BenchmarkConfig* benchmark_config) {
    const int ret = MPI_File_close(benchmark_config->mpi_log_file);
    if (ret != MPI_SUCCESS) {
        char error_string[MPI_MAX_ERROR_STRING];
        int error_string_len;
        MPI_Error_string(ret, error_string, &error_string_len);
        PRINT_DEBUG_WARN_R(my_rank, "[BENCHMARK] Error closing benchmark log file: %s\n", error_string);
    }
}


//FIXME: if application calls any collective MPI function, imbalance statistics may be incorrect (since benchmark cannot know which process waited more)
/**
 * Runs the benchmark for the given application and logs the execution time.
 * benchmark_run_c is a collective operation
 *
 * @param my_rank The rank of the current process.
 * @param benchmark_config The configuration for the benchmark.
 * @param app The application function to benchmark.
 * @param arguments The arguments to pass to the application function.
 * @param benchmark_result The structure to store benchmark results.
 */
void benchmark_run_c(const int my_rank, const BenchmarkConfig* benchmark_config, const Application app, void *arguments, BenchmarkResult* benchmark_result) {
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    double *local_times = malloc(benchmark_config->n_iterations * sizeof(double));
    double *global_times = malloc(benchmark_config->n_iterations * sizeof(double));
    double *all_local_times = NULL;

    if (my_rank == 0) {
        all_local_times = malloc(world_size * benchmark_config->n_iterations * sizeof(double));
    }

    UNIQUE_PRINT_DEBUG_INFO(my_rank, "[BENCHMARK] ***STARTING BENCHMARK***\n");

    for (int i = 0; i < benchmark_config->n_iterations; i++) {
        MPI_Barrier(MPI_COMM_WORLD);
        const double start_time = MPI_Wtime();
        app(arguments);
        const double local_end = MPI_Wtime();

        MPI_Barrier(MPI_COMM_WORLD);
        const double global_end = MPI_Wtime();

        local_times[i] = local_end - start_time;
        global_times[i] = global_end - start_time;
    }

    UNIQUE_PRINT_DEBUG_INFO(my_rank, "[BENCHMARK] ***ENDING BENCHMARK***\n");

    // Gather all local times to rank 0
    MPI_Gather(local_times, benchmark_config->n_iterations, MPI_DOUBLE,
               all_local_times, benchmark_config->n_iterations, MPI_DOUBLE,
               0, MPI_COMM_WORLD);

    if (my_rank == 0) {
        benchmark_unique_log(my_rank, benchmark_config->mpi_log_file, "\n***\n");
        benchmark_unique_log(my_rank, benchmark_config->mpi_log_file, "Results over %d iterations:\n\n", benchmark_config->n_iterations);

        // Per-iteration statistics
        for (int i = 0; i < benchmark_config->n_iterations; i++) {
            double iter_min = DBL_MAX, iter_max = 0, iter_sum = 0;
            int iter_min_rank = -1, iter_max_rank = -1;

            for (int r = 0; r < world_size; r++) {
                const double t = all_local_times[r * benchmark_config->n_iterations + i];
                if (t < iter_min) {iter_min = t; iter_min_rank = r;}
                if (t > iter_max) {iter_max = t; iter_max_rank = r;}
                iter_sum += t;
            }

            const double iter_avg = iter_sum / world_size;
            const double load_imbalance = (iter_max - iter_min) / iter_avg * 100.0;

            benchmark_unique_log(my_rank, benchmark_config->mpi_log_file,
                "Iteration %d: Global=%.6lf, Min(rank %d)=%.6lf, Max(rank %d)=%.6lf, Avg=%.6lf, Imbalance=%.2f%%\n",
                i + 1, global_times[i], iter_min_rank, iter_min, iter_max_rank, iter_max, iter_avg, load_imbalance);
        }

        free(all_local_times);
    }

    free(local_times);
    free(global_times);
}



/*void benchmark_run_c(const int my_rank, const BenchmarkConfig* benchmark_config, const Application app, void *arguments, BenchmarkResult* benchmark_result) {
    double min_time = DBL_MAX;
    double max_time = 0;
    double avg_time = 0;

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    UNIQUE_PRINT_DEBUG_INFO(my_rank, "[BENCHMARK] ***STARTING BENCHMARK***\n");

    MPI_Barrier(MPI_COMM_WORLD);
    double total_avg_time = MPI_Wtime();
    for (int i = 0; i < benchmark_config->n_iterations; i++) {
        // UNIQUE_PRINT_DEBUG_INFO(my_rank, "[BENCHMARK] Starting iteration %d/%d for benchmark '%s'\n", i + 1, benchmark_config->n_iterations, benchmark_config->name);

        const double start_time = MPI_Wtime();
        app(arguments);
        const double end_time = MPI_Wtime();
        const double local_elapsed_time = end_time - start_time;

        MPI_Barrier(MPI_COMM_WORLD);
        const double global_end_time = MPI_Wtime();
        const double global_elapsed_time = global_end_time - start_time;

        // benchmark_unique_log(my_rank, benchmark_config->log_file, "Iteration %d: %.20lf seconds\n", i + 1, elapsed_time);

        if (local_elapsed_time < min_time) {min_time = local_elapsed_time;}
        if (local_elapsed_time > max_time) {max_time = local_elapsed_time;}
        avg_time += local_elapsed_time;
    }
    total_avg_time = MPI_Wtime() - total_avg_time;
    total_avg_time /= benchmark_config->n_iterations;

    avg_time /= benchmark_config->n_iterations;

    UNIQUE_PRINT_DEBUG_INFO(my_rank, "[BENCHMARK] ***ENDING BENCHMARK***\n");

    benchmark_unique_log(my_rank, benchmark_config->mpi_log_file, "\n***\n");
    benchmark_unique_log(my_rank, benchmark_config->mpi_log_file, "Results over %d iterations, with timer resolution of %.12f seconds:\n", benchmark_config->n_iterations, MPI_Wtick());
    benchmark_unique_log(my_rank, benchmark_config->mpi_log_file, "Min Time: %.20lf seconds\n", min_time);
    benchmark_unique_log(my_rank, benchmark_config->mpi_log_file, "Max Time: %.20lf seconds\n", max_time);
    benchmark_unique_log(my_rank, benchmark_config->mpi_log_file, "Avg Time: %.20lf seconds\n", avg_time);
    benchmark_unique_log(my_rank, benchmark_config->mpi_log_file, "Amplified avg Time: %.20lf seconds\n", total_avg_time);
    benchmark_unique_log(my_rank, benchmark_config->mpi_log_file, "***\n");

    if (benchmark_result != NULL) {
        benchmark_result->min_time = min_time;
        benchmark_result->max_time = max_time;
        benchmark_result->avg_time = avg_time;
    }

}*/