//
// Created by Matteo Ranzi on 19/11/25.
//

#ifndef N_BODY_SIMULATION_RANDOM_UTILS_H
#define N_BODY_SIMULATION_RANDOM_UTILS_H

#include <vector>
#include <cstdint>
#include <random>
#include <chrono>
#include <type_traits>
#include <algorithm>

#define DOUBLE_RANDOM_MIN (-32.768)
#define DOUBLE_RANDOM_MAX (32.768)

/**
 * @brief Generates a random double value within the specified range.
 *
 * @param min The minimum value of the range (default: DOUBLE_RANDOM_MIN).
 * @param max The maximum value of the range (default: DOUBLE_RANDOM_MAX).
 * @return A random double value in the range [min, max].
 */
inline double random_double(const double min = DOUBLE_RANDOM_MIN, const double max = DOUBLE_RANDOM_MAX) {
    static std::random_device rd;
    static std::default_random_engine generator(rd());
    static std::uniform_real_distribution<double> distribution;

    return distribution(generator, std::uniform_real_distribution<double>::param_type(min, max));
}

/**
 * @brief Generates a random integer value within the specified range.
 *
 * @param start The minimum value of the range (inclusive).
 * @param end The maximum value of the range (inclusive).
 * @return A random integer value in the range [start, end].
 *
 * @note This function uses a static random device and generator for efficiency.
 * @warning Ensure the range [start, end] is valid. Behavior is undefined otherwise.
 */
inline int random_int(const int start, const int end) {
    static std::random_device rd;
    static std::mt19937 mt(rd());
    static std::uniform_int_distribution<int> distribution;

    return distribution(mt, std::uniform_int_distribution<int>::param_type(start, end));
}



/**
 * @brief Fills a vector with random unsigned integer values within the specified range.
 *
 * @tparam T The type of the elements in the vector. Must be an unsigned integer type.
 * @param vec The vector to fill with random values.
 * @param n The number of random values to generate.
 * @param min_value The minimum value of the range (inclusive).
 * @param max_value The maximum value of the range (inclusive).
 *
 * @note If n is 0, the vector will be cleared.
 * @note Uses a custom inline splitmix64 generator for fast, high-quality random numbers.
 * @warning T must be an unsigned integer type. Compilation will fail otherwise.
 */
template<typename T>
void random_fill(std::vector<T>& vec, const size_t n, const uint64_t min_value, const uint64_t max_value) {
    static_assert(std::is_unsigned_v<T>, "T must be an unsigned integer type");
    if (n == 0) { vec.clear(); return; }
    vec.resize(n);

    // Seed once (combine random_device and steady clock)
    std::random_device rd;
    uint64_t seed = (static_cast<uint64_t>(rd()) << 32) ^ rd() ^
                    static_cast<uint64_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    uint64_t state = seed;

    // Inline splitmix64 generator (fast, high-quality for non-cryptographic use)
    auto next64 = [&state]() -> uint64_t {
        state += 0x9E3779B97F4A7C15ULL;
        uint64_t z = state;
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    };

    const uint64_t range = max_value - min_value;
    if (range == 0) {
        std::fill(vec.begin(), vec.end(), static_cast<T>(min_value));
        return;
    }

    // If range == UINT64_MAX mapping via multiply would use (range+1)==0, so use raw output
    if (range == std::numeric_limits<uint64_t>::max()) {
        for (size_t i = 0; i < n; ++i) {
            vec[i] = static_cast<T>(next64() + min_value);
        }
    } else {
        const unsigned __int128 mult = static_cast<unsigned __int128>(range) + 1;
        for (size_t i = 0; i < n; ++i) {
            const uint64_t r = next64();
            // uniform mapping to [0, range]: (r * (range+1)) >> 64
            const auto scaled = static_cast<uint64_t>((static_cast<unsigned __int128>(r) * mult) >> 64);
            vec[i] = static_cast<T>(min_value + scaled);
        }
    }
}


#endif //N_BODY_SIMULATION_RANDOM_UTILS_H