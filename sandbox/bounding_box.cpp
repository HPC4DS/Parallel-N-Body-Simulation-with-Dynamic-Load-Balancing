//
// Created by Matteo Ranzi on 19/11/25.
//

#include <mpi.h>
#include <vector>

#include "debug/print_debug.h"
#include "utils/random_utils.h"
#include "benchmark.h"

#define N_BODIES 1000

// TODO SoA vs AoS benchmark (structure of arrays vs array of structures -> when to use the first and when the second?)
struct Body {
    double mass;
    double position[3];
    double velocity[3];
};

void generate_bodies(std::vector<Body>& bodies, const int n_bodies) {
    bodies.reserve(n_bodies);
    for (int i = 0; i < n_bodies; ++i) {
        const double mass = random_real(1, 100);
        const double position[3] = {random_real(1, 1000), random_real(1, 1000), random_real(1, 1000)};
        const double velocity[3] = {random_real(1, 100), random_real(1, 100), random_real(1, 100)};
        bodies.push_back(Body{mass, {position[0], position[1], position[2]}, {velocity[0], velocity[1], velocity[2]}});
    }
}

int main(int argc, char *argv[]) {
    std::vector<Body> bodies(N_BODIES);

    generate_bodies(bodies, N_BODIES);



    for (const auto i : bodies) {
        PRINT_DEBUG_INFO("value: %.35f\n", i.position[0]˛);
    }

    return 0;
}