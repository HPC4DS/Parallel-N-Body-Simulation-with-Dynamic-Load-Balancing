//
// Created by Matteo Ranzi on 19/11/25.
//

#include <mpi.h>
#include <vector>

#include "debug/print_debug.h"
#include "utils/random_utils.h"
#include "benchmark.h"

#ifdef HPC_RUN
#define LOGS_DIR getenv("HPC_JOB_LOGS_DIR")
#else
#define LOGS_DIR getenv("LOCAL_LOGS_DIR")
#endif

#define N_BODIES 10000000
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
void compute_bounding_box(const std::vector<Body>& bodies, double min_body_pos[DIMENSIONS], double max_body_pos[DIMENSIONS]);

//=============================================================================
int main(int argc, char *argv[]) {
    std::vector<Body> bodies;
    double min_pos[DIMENSIONS];
    double max_pos[DIMENSIONS];

    generate_bodies(bodies, N_BODIES);
    compute_bounding_box(bodies, min_pos, max_pos);

    PRINT_DEBUG_INFO("Bounding Box:\n");
    PRINT_DEBUG_INFO("X: [%.20f, %.20f]\n", min_pos[Body::X], max_pos[Body::X]);
    PRINT_DEBUG_INFO("Y: [%.20f, %.20f]\n", min_pos[Body::Y], max_pos[Body::Y]);
    PRINT_DEBUG_INFO("Z: [%.20f, %.20f]\n", min_pos[Body::Z], max_pos[Body::Z]);


    return 0;
}
//=============================================================================


void generate_bodies(std::vector<Body>& bodies, const int n_bodies) {
    bodies.reserve(n_bodies);
    for (int i = 0; i < n_bodies; ++i) {
        const double mass = random_real(1, 100);
        const double position[DIMENSIONS] = {random_real(1, 1000), random_real(1, 1000), random_real(1, 1000)};
        const double velocity[DIMENSIONS] = {random_real(1, 100), random_real(1, 100), random_real(1, 100)};
        bodies.emplace_back(Body{mass, {position[0], position[1], position[2]}, {velocity[0], velocity[1], velocity[2]}});
    }
}

void compute_bounding_box(const std::vector<Body>& bodies, double min_body_pos[DIMENSIONS], double max_body_pos[DIMENSIONS]) {
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