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

#define N_BODIES 1000000
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

void generate_bodies(std::vector<Body>& bodies, int n_bodies);
void linear_compute_bounding_box(const std::vector<Body>& bodies, double min_body_pos[DIMENSIONS], double max_body_pos[DIMENSIONS]);
void parallelMPI_compute_bounding_box(int rank, int comm_size, const std::vector<Body>& bodies, double min_body_pos[DIMENSIONS], double max_body_pos[DIMENSIONS]);

void benchmark_setup(int rank, BenchmarkConfig* benchmark_config, MPI_File* benchmark_log_file);

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
    benchmark_setup(my_rank, &benchmark_config, &benchmark_log_file);
    benchmark_init(my_rank, &benchmark_config);

    if (my_rank == 0) {
        PRINT_DEBUG_INFO("Generating %d bodies...\n", N_BODIES);
        generate_bodies(bodies, N_BODIES);
        PRINT_DEBUG_INFO("Bodies generated.\n");
    }

    // compute_bounding_box(bodies, min_pos, max_pos);
    //=============================================================================
    benchmark_run(my_rank, &benchmark_config, [&]() {
        if (my_rank == 0) {
            linear_compute_bounding_box(bodies, min_pos, max_pos);
        }
    },
    nullptr);
    //=============================================================================


    benchmark_finalize(my_rank, &benchmark_config);

    UNIQUE_PRINT_DEBUG_INFO(my_rank, "(Linear) Bounding Box:\n");
    UNIQUE_PRINT_DEBUG_INFO(my_rank, "X: [%.20f, %.20f]\n", min_pos[Body::X], max_pos[Body::X]);
    UNIQUE_PRINT_DEBUG_INFO(my_rank, "Y: [%.20f, %.20f]\n", min_pos[Body::Y], max_pos[Body::Y]);
    UNIQUE_PRINT_DEBUG_INFO(my_rank, "Z: [%.20f, %.20f]\n", min_pos[Body::Z], max_pos[Body::Z]);

    MPI_Finalize();


    return 0;
}
//*********************************************************************************

void benchmark_setup(int rank, BenchmarkConfig* benchmark_config, MPI_File* benchmark_log_file) {
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

void parallelMPI_compute_bounding_box(int rank, int comm_size, const std::vector<Body>& bodies, double min_body_pos[DIMENSIONS], double max_body_pos[DIMENSIONS]) {

}
