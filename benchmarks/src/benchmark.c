//
// Created by Matteo Ranzi on 02/11/25.
//

#include "benchmark.h"

#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <stdarg.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
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
static void compose_benchmark_filename(const int my_rank, char *benchmark_log_filename, const char *logs_dir) {
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
            if (all_names != NULL) {
                char *recv_name = all_names + r * MPI_MAX_PROCESSOR_NAME;
                benchmark_unique_log(my_rank, file, "Rank: %d > node: %s\n", r, recv_name);
            }
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


void create_benchmark_dir_if_not_exists(const int my_rank, const char *logs_dir) {
    if (my_rank == 0) {
        if (logs_dir == NULL || logs_dir[0] == '\0') {
            PRINT_DEBUG_ERROR_R(my_rank, "[BENCHMARK] Invalid logs_dir\n");
            benchmark_abort(1);
        }

        if (mkdir(logs_dir, 0775) != 0) {
            if (errno != EEXIST) {
                PRINT_DEBUG_ERROR_R(my_rank, "[BENCHMARK] Failed to create benchmark log directory `%s`: %s\n", logs_dir, strerror(errno));
                benchmark_abort(1);
            } else {
                PRINT_DEBUG_INFO_R(my_rank, "[BENCHMARK] Benchmark log directory already exists: `%s`\n", logs_dir);
            }
        } else {
            PRINT_DEBUG_INFO_R(my_rank, "[BENCHMARK] Created benchmark log directory: `%s`\n", logs_dir);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
}


void default_benchmark_config(BenchmarkConfig* benchmark_config) {
char logs_dir_buf[MAX_LOG_FILENAME_LENGTH];
    const char *bdir = (BENCHMARK_DIR ? BENCHMARK_DIR : "");
    const char *appid = (APPLICATION_ID ? APPLICATION_ID : "default_app");
    const int len = snprintf(logs_dir_buf, sizeof(logs_dir_buf), "%s/%s", bdir, appid);
    if (len < 0 || len >= (int)sizeof(logs_dir_buf)) {
        benchmark_abort(1);
    }
    benchmark_config->logs_dir = strdup(logs_dir_buf);
    benchmark_config->repetitions = 15;
    benchmark_config->min_time = 0.5; // seconds
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

    compose_benchmark_filename(my_rank, benchmark_log_filename, benchmark_config->logs_dir);
    create_benchmark_dir_if_not_exists(my_rank, benchmark_config->logs_dir);
    open_benchmark_file(my_rank, benchmark_log_filename, benchmark_config->mpi_log_file);

    UNIQUE_PRINT_DEBUG_INFO(my_rank, "[BENCHMARK] Output log to: %s\n", benchmark_log_filename);

    write_benchmark_log_header(my_rank, benchmark_config->mpi_log_file, benchmark_config->description);

    benchmark_unique_log(my_rank, benchmark_config->mpi_log_file, "%s;repetition;iteration;global_time;time_min_rank;time_min;time_max_rank;time_max\n", benchmark_config->sweep_name);
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


static int next_iteration_count(int current, double scale) {
    double raw = current * scale;

    // Enforce growth constraints:
    // - must increase by at least 1
    // - must increase by at most 10x
    const double max_growth = current * 10.0;
    if (raw < current + 1) {
        raw = current + 1;
    } else if (raw > max_growth) {
        raw = max_growth;
    }

    return (int) raw;
}

double run_benchmark_iteration(const int my_rank, const int iterations, const Application preHook, void *preHookArgs, const Application app, void *appArgs, const Application postHook, void *postHookArgs) {
    MPI_Barrier(MPI_COMM_WORLD);
    preHook(preHookArgs);

    double elapsed_time = 0.0;
    for (int i = 0; i < iterations; i++) {
        MPI_Barrier(MPI_COMM_WORLD);
        const double start_time = MPI_Wtime();

        app(appArgs);

        MPI_Barrier(MPI_COMM_WORLD);
        const double end_time = MPI_Wtime();

        elapsed_time += end_time - start_time;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    postHook(postHookArgs);

    return elapsed_time;
}

int determine_benchmark_iterations(const int my_rank, const REPETITION_STRATEGY repetition_strategy, const BenchmarkConfig* benchmark_config, const Application preHook, void *preHookArgs, const Application app, void *appArgs, const Application postHook, void *postHookArgs) {
    if (repetition_strategy != MIN_TIME) {
        return benchmark_config->max_iterations;
    }

    int iterations = 1;
    double time = run_benchmark_iteration(my_rank, iterations, preHook, preHookArgs, app, appArgs, postHook, postHookArgs);
    UNIQUE_PRINT_DEBUG_INFO(my_rank, "[BENCHMARK] Initial single iteration time: %.6lf seconds\n", time);

    if (time >= benchmark_config->min_time) {
        UNIQUE_PRINT_DEBUG_INFO(my_rank, "[BENCHMARK] Single iteration meets min_time %.2lf seconds\n", benchmark_config->min_time);
        return 1;
    }

    while (time < benchmark_config->min_time) {
        if (my_rank == 0) {
            const double scale = benchmark_config->min_time / time * 1.1; // add 10% margin to speed up convergence
            iterations = next_iteration_count(iterations, scale);
        } else {
            MPI_Allgather(&iterations, 1, MPI_INT, &iterations, 1, MPI_INT, MPI_COMM_WORLD);
        }

        UNIQUE_PRINT_DEBUG_INFO(my_rank, "[BENCHMARK] Testing %d iterations to reach min_time %.2lf seconds\n", iterations, benchmark_config->min_time);
        time = run_benchmark_iteration(my_rank, iterations, preHook, preHookArgs, app, appArgs, postHook, postHookArgs);
        UNIQUE_PRINT_DEBUG_INFO(my_rank, "[BENCHMARK] Total global time for %d iterations: %.6lf seconds\n", iterations, time);
    }

    UNIQUE_PRINT_DEBUG_INFO(my_rank, "[BENCHMARK] Determined %d iterations to meet min_time %.2lf seconds\n", iterations, benchmark_config->min_time);

    return iterations;
}

//FIXME: if benchmarked application calls any collective MPI function, imbalance statistics may be incorrect (since benchmark cannot know which process waited more)
/**
 * Runs the benchmark for the given application and logs the execution time.
 * benchmark_run_c is a collective operation
 *
 * @param my_rank The rank of the current process.
 * @param repetition_strategy The strategy for repeating the benchmark.
 * @param benchmark_config The configuration for the benchmark.
 * @param preHook The pre-hook function to execute before the application function.
 * @param preHookArgs The arguments to pass to the pre-hook function.
 * @param app The application function to benchmark.
 * @param appArgs The arguments to pass to the application function.
 * @param postHook The post-hook function to execute after the application function.
 * @param postHookArgs The arguments to pass to the post-hook function.
 * @param benchmark_result The structure to store benchmark results.
 */
void benchmark_run(const int my_rank, const REPETITION_STRATEGY repetition_strategy, const BenchmarkConfig* benchmark_config, const Application preHook, void *preHookArgs, const Application app, void *appArgs, const Application postHook, void *postHookArgs, BenchmarkResult* benchmark_result) {
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    UNIQUE_PRINT_DEBUG_INFO(my_rank, "[BENCHMARK] -> Determining iterations for sweep: %d\n", benchmark_config->sweep_value);
    const int iterations = determine_benchmark_iterations(my_rank, repetition_strategy, benchmark_config, preHook, preHookArgs, app, appArgs, postHook, postHookArgs);

    double global_time;
    double local_time;
    double *all_local_times = NULL;

    if (my_rank == 0) {
        all_local_times = malloc(world_size * sizeof(double));
    }

    //=============================================================================
    UNIQUE_PRINT_DEBUG_INFO(my_rank, "[BENCHMARK] ***STARTING BENCHMARK (sweep: %d)***\n", benchmark_config->sweep_value);

    //fixme: currently repetitions are independent, we could gather statistics across repetitions (or perform unique time measurement across repetitions)
    // TODO
    for (int rep = 1; rep <= benchmark_config->repetitions; rep++) {
        for (int iter = 1; iter <= iterations; iter++) {
            MPI_Barrier(MPI_COMM_WORLD);
            preHook(preHookArgs);

            MPI_Barrier(MPI_COMM_WORLD);
            const double start_time = MPI_Wtime();
            app(appArgs);
            const double local_end = MPI_Wtime();

            MPI_Barrier(MPI_COMM_WORLD);
            const double global_end = MPI_Wtime();

            local_time = local_end - start_time;
            global_time = global_end - start_time;

            MPI_Barrier(MPI_COMM_WORLD);
            postHook(postHookArgs);

            MPI_Gather(&local_time, 1, MPI_DOUBLE, all_local_times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

            // Per-iteration statistics
            if (my_rank == 0) {
                double iter_min = DBL_MAX, iter_max = 0;
                int iter_min_rank = -1, iter_max_rank = -1;

                for (int r = 0; r < world_size; r++) {
                    const double t = all_local_times[r];
                    if (t <= iter_min) {iter_min = t; iter_min_rank = r;}
                    if (t >= iter_max) {iter_max = t; iter_max_rank = r;}
                }
                benchmark_unique_log(my_rank, benchmark_config->mpi_log_file,"%d;%d;%d;%.20lf;%d;%.20lf;%d;%.20lf\n",
                                     benchmark_config->sweep_value, rep, iter, global_time, iter_min_rank, iter_min, iter_max_rank, iter_max);
            }
        }
    }
    UNIQUE_PRINT_DEBUG_INFO(my_rank, "[BENCHMARK] ***ENDING BENCHMARK (sweep: %d)***\n", benchmark_config->sweep_value);
    //=============================================================================

    free(all_local_times);
}