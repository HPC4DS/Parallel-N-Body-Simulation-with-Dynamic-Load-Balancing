//
// Created by Matteo Ranzi on 19/11/25.
//

#include <iostream>
#include <mpi.h>
#include <ostream>
#include <vector>

#include "debug/print_debug.h"
#include "utils/random_utils.h"
#include "benchmark.hpp"
#include "debug/unique_print_debug.h"

#ifdef HPC_RUN
#define LOGS_DIR getenv("HPC_JOB_LOGS_DIR")
#else
#define LOGS_DIR getenv("LOCAL_LOGS_DIR")
#endif

#define N_BODIES 80000000
#define DIMENSIONS 3

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


void generate_bodies(std::vector<Body>& bodies, int n_bodies);
void linear_compute_bounding_box(const std::vector<Body>& bodies, double min_body_pos[DIMENSIONS], double max_body_pos[DIMENSIONS]);
void parallelMPI_compute_bounding_box(const std::vector<Body>& bodies, double min_body_pos[DIMENSIONS], double max_body_pos[DIMENSIONS]);

void setup_bm(BenchmarkConfig* benchmark_config, MPI_File* benchmark_log_file);

//*********************************************************************************
int main(int argc, char *argv[]) {
    std::vector<Body> bodies;
    double min_pos[DIMENSIONS];
    double max_pos[DIMENSIONS];

    int comm_size;
    int my_rank;

    MPI_File benchmark_log_file;


    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    BenchmarkConfig benchmark_config;
    setup_bm(&benchmark_config, &benchmark_log_file);
    benchmark_init(my_rank, &benchmark_config);

    if (my_rank == 0) {
        PRINT_DEBUG_INFO("Generating %d bodies...\n", N_BODIES);
        generate_bodies(bodies, N_BODIES);
        PRINT_DEBUG_INFO("Bodies generated.\n");
    }

    {
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
    // compute_bounding_box(bodies, min_pos, max_pos);
    //=============================== BENCHMARK ===================================
    benchmark_run(my_rank, &benchmark_config, [&]() {
        /*if (my_rank == 0) {
            linear_compute_bounding_box(bodies, min_pos, max_pos);
        }*/
            parallelMPI_compute_bounding_box(bodies, min_pos, max_pos);
    },
    nullptr);
    //=============================================================================
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Type_free(&MPI_BODY);

    benchmark_finalize(my_rank, &benchmark_config);

    UNIQUE_PRINT_DEBUG_INFO(my_rank, "(Linear) Bounding Box:\n");
    UNIQUE_PRINT_DEBUG_INFO(my_rank, "X: [%.20f, %.20f]\n", min_pos[Body::X], max_pos[Body::X]);
    UNIQUE_PRINT_DEBUG_INFO(my_rank, "Y: [%.20f, %.20f]\n", min_pos[Body::Y], max_pos[Body::Y]);
    UNIQUE_PRINT_DEBUG_INFO(my_rank, "Z: [%.20f, %.20f]\n", min_pos[Body::Z], max_pos[Body::Z]);

    MPI_Finalize();


    return 0;
}
//*********************************************************************************

void setup_bm(BenchmarkConfig* benchmark_config, MPI_File* benchmark_log_file) {
    benchmark_config->logs_dir = LOGS_DIR;
    benchmark_config->mpi_log_file = benchmark_log_file;
    strcpy(benchmark_config->description, "Linear Bounding Box Computation. N Bodies: ");
    strcat(benchmark_config->description, std::to_string(N_BODIES).c_str());
    benchmark_config->n_iterations = 5;
}

void generate_bodies(std::vector<Body>& bodies, const int n_bodies) {
    bodies.resize(n_bodies);
    for (int i = 0; i < n_bodies; ++i) {
        bodies[i].mass = random_real(1, 100);
        for (int axis = 0; axis < DIMENSIONS; ++axis) {
            bodies[i].position[axis] = random_real(1, 1000);
            bodies[i].velocity[axis] = random_real(1, 100);
        }
    }
}

void linear_compute_bounding_box(const std::vector<Body>& bodies, double min_body_pos[DIMENSIONS], double max_body_pos[DIMENSIONS]) {
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

void parallelMPI_compute_bounding_box(const std::vector<Body>& bodies, double min_body_pos[DIMENSIONS], double max_body_pos[DIMENSIONS]) {

    // Step 1: compute counts and displacements
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    const int N = static_cast<int>(bodies.size());   // only valid on rank 0

    std::vector<int> counts(world_size);
    std::vector<int> displacements(world_size);

    const int base = N / world_size;
    const int extra = N % world_size;

    if (world_rank == 0) {
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

    if (world_rank == 0) {
        local_count = counts[0];
        for (int r = 1; r < world_size; r++)
            MPI_Send(&counts[r], 1, MPI_INT, r, 0, MPI_COMM_WORLD);
    } else {
        MPI_Recv(&local_count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // Step 3: received buffer holder
    std::vector<Body> local_bodies(local_count);

    // Step 4: scatterv
    MPI_Scatterv(
    bodies.data(),
    counts.data(),
    displacements.data(),
    MPI_BODY,
    local_bodies.data(),
    local_count,
    MPI_BODY,
    0,
    MPI_COMM_WORLD
    );


    // Step 5: local bounding box computation
    double local_min_pos[DIMENSIONS];
    double local_max_pos[DIMENSIONS];

    linear_compute_bounding_box(local_bodies, local_min_pos, local_max_pos);

    // Step 6: global reduction
    // comment out for processes profiling purposes in benchmark
/*
    {
        MPI_Reduce(
        local_min_pos,
        min_body_pos,
        DIMENSIONS,
        MPI_DOUBLE,
        MPI_MIN,
        0,
        MPI_COMM_WORLD
        );
    MPI_Reduce(
        local_max_pos,
        max_body_pos,
        DIMENSIONS,
        MPI_DOUBLE,
        MPI_MAX,
        0,
        MPI_COMM_WORLD
        );
    }
*/
}
