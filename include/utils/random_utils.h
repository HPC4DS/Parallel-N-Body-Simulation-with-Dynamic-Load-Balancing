//
// Created by Matteo Ranzi on 19/11/25.
//

#ifndef N_BODY_SIMULATION_RANDOM_UTILS_H
#define N_BODY_SIMULATION_RANDOM_UTILS_H

#include <random>

#define DOUBLE_RANDOM_MIN (-32.768)
#define DOUBLE_RANDOM_MAX (32.768)

inline double random_real(const double min = DOUBLE_RANDOM_MIN, const double max = DOUBLE_RANDOM_MAX) {
    static std::random_device rd;
    static std::default_random_engine generator(rd());
    static std::uniform_real_distribution<double> distribution;

    return distribution(generator, std::uniform_real_distribution<double>::param_type(min, max));
}

//TODO check if it works correctly
inline int random_int(const int start, const int end) {
    static std::random_device rd;
    static std::mt19937 mt(rd());
    static std::uniform_int_distribution<int> distribution;

    return distribution(mt, std::uniform_int_distribution<int>::param_type(start, end));
}

#endif //N_BODY_SIMULATION_RANDOM_UTILS_H