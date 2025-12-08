//
// Created by Matteo Ranzi on 08/12/25.
//

#ifndef N_BODY_SIMULATION_LINEAR_ALGORITHMS_HPP
#define N_BODY_SIMULATION_LINEAR_ALGORITHMS_HPP

#include <iostream>
#include <vector>
#include <cstdint>
#include <random>
#include <type_traits>
#include <algorithm>


/**
 * @brief Performs a radix sort on a range of unsigned integers.
 *
 * @tparam BITS The number of bits to process per pass (default: 10).
 * @tparam RandomAccessIterator The type of the iterator, which must point to unsigned integers.
 * @param first The beginning of the range to sort.
 * @param last The end of the range to sort.
 * @param reverse If true, sorts the range in descending order; otherwise, sorts in ascending order.
 *
 * @note This function requires the value type of the iterator to be an unsigned integer.
 * @note The algorithm uses a stable, in-place sorting approach with a time complexity of O(n * log(max_value)).
 * @warning Ensure the range [first, last) is valid and the value type is unsigned. Behavior is undefined otherwise.
 */
template<const int BITS = 10, typename RandomAccessIterator>
inline void radix_sort(RandomAccessIterator first, RandomAccessIterator last, const bool reverse = false)
{
    using value_type = typename std::iterator_traits<RandomAccessIterator>::value_type;
    static_assert(std::is_unsigned_v<value_type>, "radix_sort requires unsigned value_type (use uint64_t).");

    if (first == last) return;
    const size_t n = std::distance(first, last);
    if (n <= 1) return;

    constexpr int RADIX = 1 << BITS;
    constexpr value_type MASK = (RADIX - 1);
    const int max_bits = [&]() {
        auto it = max_element(first, last);
        value_type max_value = (it == last) ? 0 : *it;
        if (max_value == 0) return 0;
        return 64 - __builtin_clzll(max_value); // number of significant bits
    }();

    int passes = (max_bits + BITS - 1) / BITS;
    if (passes == 0) passes = 1;

    std::vector<value_type> src;
    src.assign(first, last);           // assume cost ok for caller (caller used vector)
    std::vector<value_type> dst(n);

    std::vector<size_t> count(RADIX);
    std::vector<size_t> offset(RADIX);

    int shift = 0;
    for (int pass = 0; pass < passes; ++pass, shift += BITS) {
        std::fill(count.begin(), count.end(), 0);

        // counting phase
        for (size_t i = 0; i < n; ++i) {
            auto idx = static_cast<uint32_t>((src[i] >> shift) & MASK);
            if (reverse) idx ^= static_cast<uint32_t>(MASK);
            ++count[idx];
        }

        // prefix-sum -> offsets
        size_t sum = 0;
        for (int i = 0; i < RADIX; ++i) {
            offset[i] = sum;
            sum += count[i];
        }

        // scatter
        for (size_t i = 0; i < n; ++i) {
            auto idx = static_cast<uint32_t>((src[i] >> shift) & MASK);
            if (reverse) idx ^= static_cast<uint32_t>(MASK);
            dst[offset[idx]++] = src[i];
        }

        // swap src/dst for next pass
        src.swap(dst);
    }

    // ensure sorted data is written back to original range
    std::copy(src.begin(), src.end(), first);
}


#endif //N_BODY_SIMULATION_LINEAR_ALGORITHMS_HPP