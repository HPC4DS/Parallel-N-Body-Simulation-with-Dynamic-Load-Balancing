//
// Created by Matteo Ranzi on 19/11/25.
//

#include <iostream>
#include <mpi.h>
#include <vector>

#include "debug/print_debug.h"
#include "utils/random_utils.h"
#include "benchmark.hpp"
#include "debug/unique_print_debug.h"

#ifdef HPC_RUN
#define LOGS_DIR getenv("HPC_JOB_LOGS_DIR")
#else
#define BENCHMARK_DIR getenv("LOCAL_LOGS_DIR")
#endif

#define N_BODIES_DEFAULT 10000000
#define DIMENSIONS 3

#define BENCHMARK_LINEAR 0
#define BENCHMARK_PARALLEL 1
#define BENCHMARK_TYPE BENCHMARK_PARALLEL

// TODO SoA vs AoS benchmark (structure of arrays vs array of structures -> when to use the first and when the second?)
struct Body {
    double mass;
    double position[DIMENSIONS];
    double velocity[DIMENSIONS];

    enum Coordinate {
        X = 0,
        Y = 1,
        Z = 2
    };
};
MPI_Datatype MPI_BODY;
void register_custom_MPI_Datatypes();
void unregister_custom_MPI_Datatypes();


void generate_bodies(std::vector<Body>& bodies, int n_bodies);
void scatter_bodies(const std::vector<Body>& bodies, std::vector<Body>& local_bodies);

void compute_bounding_box_serial(const std::vector<Body>& bodies, double min_body_pos[DIMENSIONS], double max_body_pos[DIMENSIONS]);
void compute_bounding_box_parallel(const std::vector<Body>& local_bodies, double min_body_pos[DIMENSIONS], double max_body_pos[DIMENSIONS]);

void setup_bm_cnfg(BenchmarkConfig& benchmark_config, MPI_File& benchmark_log_file, const char *benchmark_description);
void linear_benchmark(int my_rank, BenchmarkConfig &benchmark_config, MPI_File &benchmark_log_file, const std::vector<Body>& bodies, double min_pos[DIMENSIONS], double max_pos[DIMENSIONS]);
void parallel_benchmark(int my_rank, BenchmarkConfig &benchmark_config, MPI_File &benchmark_log_file, const std::vector<Body>& bodies, double min_pos[DIMENSIONS], double max_pos[DIMENSIONS]);

inline void check_mpi(int my_rank, int result, const char* operation);

//*********************************************************************************
int main(int argc, char *argv[]) {
    std::vector<Body> bodies;
    double min_pos[DIMENSIONS];
    double max_pos[DIMENSIONS];
    int comm_size;
    int my_rank;
    MPI_File benchmark_log_file;
    BenchmarkConfig benchmark_config;


    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);


    if (my_rank == 0) {
        PRINT_DEBUG_INFO("Generating %d bodies...\n", N_BODIES_DEFAULT);
        generate_bodies(bodies, N_BODIES_DEFAULT);
        PRINT_DEBUG_INFO("Bodies generated.\n");
    }

    //======================================= BENCHMARK ===========================================

    if constexpr (BENCHMARK_TYPE == BENCHMARK_LINEAR) {
        linear_benchmark(my_rank, benchmark_config, benchmark_log_file, bodies, min_pos, max_pos);
    } else {
        parallel_benchmark(my_rank, benchmark_config, benchmark_log_file, bodies, min_pos, max_pos);
    }
    //=============================================================================================


    UNIQUE_PRINT_DEBUG_INFO(my_rank, "Bounding Box:\n");
    UNIQUE_PRINT_DEBUG_INFO(my_rank, "X: [%.20f, %.20f]\n", min_pos[Body::X], max_pos[Body::X]);
    UNIQUE_PRINT_DEBUG_INFO(my_rank, "Y: [%.20f, %.20f]\n", min_pos[Body::Y], max_pos[Body::Y]);
    UNIQUE_PRINT_DEBUG_INFO(my_rank, "Z: [%.20f, %.20f]\n", min_pos[Body::Z], max_pos[Body::Z]);



    MPI_Finalize();

    return 0;
}
//*********************************************************************************

void setup_bm_cnfg(BenchmarkConfig& benchmark_config, MPI_File& benchmark_log_file, const char *benchmark_description) {
    benchmark_config.logs_dir = BENCHMARK_DIR;
    benchmark_config.mpi_log_file = &benchmark_log_file;
    strcpy(benchmark_config.description, benchmark_description);
    benchmark_config.max_iterations = 5;
}

void generate_bodies(std::vector<Body>& bodies, const int n_bodies) {
    bodies.resize(n_bodies);
    for (int i = 0; i < n_bodies; ++i) {
        bodies[i].mass = random_double(1, 100);
        for (int axis = 0; axis < DIMENSIONS; ++axis) {
            bodies[i].position[axis] = random_double(1, 1000);
            bodies[i].velocity[axis] = random_double(1, 100);
        }
    }
}

void compute_bounding_box_serial(const std::vector<Body>& bodies, double min_body_pos[DIMENSIONS], double max_body_pos[DIMENSIONS]) {
    for (int axis = 0; axis < DIMENSIONS; ++axis) {
        min_body_pos[axis] = bodies[0].position[axis];
        max_body_pos[axis] = bodies[0].position[axis];
    }

    for (const auto& body : bodies) {
        for (int axis = 0; axis < DIMENSIONS; ++axis) {
            if (body.position[axis] < min_body_pos[axis]) min_body_pos[axis] = body.position[axis];
            if (body.position[axis] > max_body_pos[axis]) max_body_pos[axis] = body.position[axis];
        }
    }
}

void scatter_bodies(const std::vector<Body>& bodies, std::vector<Body>& local_bodies) {

    // Step 1: compute counts and displacements
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    const int N = static_cast<int>(bodies.size());   // only valid on rank 0

    std::vector<int> counts(world_size);
    std::vector<int> displacements(world_size);

    const int base = N / world_size;
    const int extra = N % world_size;

    if (my_rank == 0) {
        for (int r = 0; r < world_size; r++) {
            counts[r] = base + (r < extra ? 1 : 0);
            // PRINT_DEBUG_INFO("Rank %d: count = %d\n", r, counts[r]);
        }

        displacements[0] = 0;
        for (int r = 1; r < world_size; r++) {
            displacements[r] = displacements[r-1] + counts[r-1];
        }
    }

    // Step 2: communicate local count
    int local_count;

    check_mpi(my_rank, MPI_Scatter(counts.data(), 1, MPI_INT, &local_count, 1, MPI_INT, 0, MPI_COMM_WORLD), "Scatter counts");

    // Step 3: received buffer holder
    local_bodies.resize(local_count);

    // Step 4: scatterv
    check_mpi(my_rank, MPI_Scatterv(bodies.data(),counts.data(),displacements.data(),MPI_BODY,local_bodies.data(),local_count,MPI_BODY,0,MPI_COMM_WORLD), "Scatterv bodies");
}
void compute_bounding_box_parallel(const std::vector<Body>& local_bodies, double min_body_pos[DIMENSIONS], double max_body_pos[DIMENSIONS]) {
    // Step 5: local bounding box computation
    double local_min_pos[DIMENSIONS];
    double local_max_pos[DIMENSIONS];

    compute_bounding_box_serial(local_bodies, local_min_pos, local_max_pos);

    // Step 6: global reduction
    // comment out for processes profiling purposes in benchmark
    MPI_Reduce(local_min_pos, min_body_pos, DIMENSIONS, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_max_pos, max_body_pos,DIMENSIONS, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
}

void register_custom_MPI_Datatypes() {
    constexpr int block_lengths[3] = {
        1,              // mass
        DIMENSIONS,     // position[]
        DIMENSIONS      // velocity[]
    };

    MPI_Aint offsets[3];
    constexpr MPI_Datatype types[3] = {
        MPI_DOUBLE,     // mass
        MPI_DOUBLE,     // position[]
        MPI_DOUBLE      // velocity[]
    };

    offsets[0] = offsetof(Body, mass);
    offsets[1] = offsetof(Body, position);
    offsets[2] = offsetof(Body, velocity);

    MPI_Type_create_struct(3, block_lengths, offsets, types, &MPI_BODY);
    MPI_Type_commit(&MPI_BODY);
}
void unregister_custom_MPI_Datatypes() {
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Type_free(&MPI_BODY);
}

void linear_benchmark(const int my_rank, BenchmarkConfig &benchmark_config, MPI_File &benchmark_log_file, const std::vector<Body>& bodies, double min_pos[DIMENSIONS], double max_pos[DIMENSIONS]) {
    char description[256];
    snprintf(description, 256, "Linear Bounding Box Computation. N Bodies: %d", static_cast<int>(bodies.size()));
    setup_bm_cnfg(benchmark_config, benchmark_log_file, description);
    benchmark_init(my_rank, &benchmark_config);
    benchmark_run(my_rank, &benchmark_config, [&]() {
        if (my_rank == 0) {
            compute_bounding_box_serial(bodies, min_pos, max_pos);
        }
    },
    nullptr);
    benchmark_finalize(my_rank, &benchmark_config);
}

void parallel_benchmark(const int my_rank, BenchmarkConfig &benchmark_config, MPI_File &benchmark_log_file, const std::vector<Body>& bodies, double min_pos[DIMENSIONS], double max_pos[DIMENSIONS]) {
    register_custom_MPI_Datatypes();

    std::vector<Body> local_bodies;
    UNIQUE_PRINT_DEBUG_INFO(my_rank, "Scattering bodies...\n");
    scatter_bodies(bodies, local_bodies);
    UNIQUE_PRINT_DEBUG_INFO(my_rank, "Bodies scattered.\n");

    char description[256];
    snprintf(description, 256, "Parallel Bounding Box Computation + reduction (already scattered). N Bodies: %d", static_cast<int>(bodies.size()));
    setup_bm_cnfg(benchmark_config, benchmark_log_file, description);
    benchmark_init(my_rank, &benchmark_config);
    benchmark_run(my_rank, &benchmark_config, [&]() {
            compute_bounding_box_parallel(local_bodies, min_pos, max_pos);
    },
    nullptr);
    benchmark_finalize(my_rank, &benchmark_config);

    unregister_custom_MPI_Datatypes();
}

inline void check_mpi(const int my_rank, const int result, const char* operation) {
    if (result != MPI_SUCCESS) {
        char error_string[MPI_MAX_ERROR_STRING];
        int length;
        MPI_Error_string(result, error_string, &length);
        UNIQUE_PRINT_DEBUG_ERROR(my_rank, "MPI Error - %s: %s", operation, error_string);
        MPI_Abort(MPI_COMM_WORLD, result);
    }
}