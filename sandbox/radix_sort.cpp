//
// Created by Matteo Ranzi on 06/12/25.
//

#include <iostream>
#include <vector>
#include <cstdint>
#include <random>
#include <chrono>
#include <type_traits>
#include <algorithm>

// https://codeforces.com/blog/entry/122438
/*
 * BITS = 4 is faster than BITS = 8
Short explanation — why BITS = 4 was faster here:
The implementation uses std::vector per bucket (P[PART]) and push_back. With BITS = 8 you have 256 buckets; push_back across many buckets causes pointer-chasing, lots of small growths/allocations/fragmentation and very poor streaming writes.
BITS = 4 uses only 16 buckets, so fewer vectors, less scattered writes and better effective locality — that overcame the extra passes.
Fix: use a single auxiliary contiguous buffer + a counting (histogram) phase per pass (prefix-sum + scatter). That removes per-bucket allocations, makes writes streaming, and restores the expected benefit of larger BITS.
Brief description of the replacement radix_sort:
Assumes contiguous random-access range (e.g. std::vector<T>).
Works for unsigned integer element types (use uint64_t in your case).
For each pass: build counts, compute offsets (prefix sums), scatter into dst, then swap src/dst. At the end, ensure result is in the original range.


radix implementation becomes memory-bandwidth and memory-pressure bound as N grows — extra passes, extra buffer allocations and writes make it scale worse than std::sort once data no longer fits comfortably in cache/RAM.
Why this happens (concise):
Pass count is fixed by value width: passes = ceil((log2(max)+1) / BITS). With BITS = 4 and max = 1 << 63 you do ~16 full passes over the array. Each pass reads and writes the whole array → huge memory traffic.
Memory traffic grows linearly with N. Example (64-bit elements, 16 passes): per-pass read+write ≈ 16·N bytes, total ≈ 16·passes·N ≈ 256·N bytes. For N=500M: ~128 GB of memory traffic; for N=1B: ~256 GB. That can saturate system memory bandwidth.
Extra memory footprint: your bucket std::vector arrays hold the whole array each pass (additional ~N·8 bytes). Doubling N may exceed available RAM and force paging or heavy OS memory pressure, killing performance.
Allocator costs and reallocation: many push_back calls cause reallocations and heap fragmentation; even if capacity persists, managing many large vectors increases overhead.
Cache effects: larger N moves working set out of caches/DRAM locality; std::sort (in-place, comparison-based) uses less extra memory and benefits from locality, so it can win once memory bandwidth / cache misses dominate.
Threading/CPU: std::sort's CPU-bound work can be faster relative to the memory-bound radix at large N.
How to confirm:
Monitor RAM / swap while running (Activity Monitor, top, vm_stat on macOS) and watch page faults.
Profile with Instruments (Allocations, Time Profiler, VM Tracker) or a memory-bandwidth benchmark to see saturation.
Measure per-pass time or add counters to see time spent on bucket push_back vs copying back.
Practical fixes / improvements:
Use a single auxiliary buffer and counting histogram (two-phase stable distribution) instead of many std::vector buckets. That removes per-element allocations and reduces write amplification. (Classic LSD radix implementation: count[], prefix-sum, scatter into temp[], swap buffers.)
Increase BITS (e.g. 8) to reduce passes (fewer passes → less total memory traffic). Tradeoff: larger PART and count array; prefer single-buffer approach when increasing BITS.
Pre-reserve bucket capacities if you keep the vector-buckets approach: reserve(N / PART) for each bucket once, to avoid repeated reallocations.
Avoid extra copies: write directly into an auxiliary contiguous array and swap pointers each pass instead of push_back into many small vectors.
Ensure you have sufficient RAM; if not, reduce N or use external sorting / disk-based approaches.
Consider parallel radix that chunks data and uses per-thread histograms to reduce contention and improve throughput (but still memory-bandwidth bound).
Small concrete suggestion: implement two-phase counting + scatter with one temp buffer and BITS = 8. That typically reduces memory traffic and avoids many tiny allocations, and will scale much better to 1B elements.
 */
template<const bool negative = true, const int BITS = 4, typename RandomAccessIterator>
inline void codeforces_radix_sort(RandomAccessIterator _first, RandomAccessIterator _last, const bool reverse = false)
{
    typedef typename std::iterator_traits<RandomAccessIterator>::value_type value_type;
    if(negative)
    {
        if(!reverse) {
            RandomAccessIterator _middle = partition(_first, _last, [](value_type a) -> bool {return a<0;});
            for(auto it = _first; it != _middle; ++it)*it *= -1;
            codeforces_radix_sort<false>(_first, _middle, true);
            for(auto it = _first; it != _middle; ++it)*it *= -1;
            codeforces_radix_sort<false>(_middle, _last);
        } else {
            RandomAccessIterator _middle = partition(_first, _last, [](value_type a) -> bool {return a>=0;});
            codeforces_radix_sort<false>(_first, _middle, true);
            for(auto it = _middle; it != _last; ++it)*it *= -1;
            codeforces_radix_sort<false>(_middle, _last);
            for(auto it = _middle; it != _last; ++it)*it *= -1;
        }
        return;
    }
    constexpr int PART = (1<<BITS);
    constexpr int FULL = (1<<BITS)-1;
    int shift = 0;
    const auto MAX = log2(*max_element(_first, _last));
    std::vector<value_type> P[PART];

    while(MAX >= shift)
    {
        for(auto it = _first; it != _last; ++it)P[(((*it)>>shift)&FULL)^(reverse?FULL:0)].push_back(*it);
        auto it = _first;
        for(int i = 0; i < PART; i++)
        {
            for(auto u: P[i]) *it++ = u;
            P[i].clear();
        }
        shift += BITS;
    }
}


// only works with max vector size: 2^32
template<const int BITS = 8, typename RandomAccessIterator>
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

int main() {
    std::vector<uint64_t> numbers;

    std::cout << "Filling vector with random numbers..." << std::endl;
    random_fill(numbers, 500000000, 0LU, /*9209641988550369538*/ 1LU << 63);



    std::cout << "Sorting vector elements..." << std::endl;
    const auto start = std::chrono::high_resolution_clock::now();
    codeforces_radix_sort(numbers.begin(), numbers.end());
    radix_sort<11>(numbers.begin(), numbers.end());
    //std::sort(numbers.begin(), numbers.end());
    const auto elapsed = std::chrono::high_resolution_clock::now() - start;

    const unsigned long microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
    std::cout << "Sorting completed in: " << microseconds << " us ["
    << static_cast<double>(microseconds) / 1000000.0 << " s]" << std::endl;


    // Verify sorting
    for (size_t i = 1; i < numbers.size(); ++i) {
        if (numbers[i - 1] > numbers[i]) {
            std::cerr << "Sorting failed at index " << i << ": " << numbers[i - 1] << " > " << numbers[i] << std::endl;
            return 1;
        }
    }

    std::cout << "\nSorting successful!" << std::endl;

    return 0;
}

/*
 * // TODO test better BITS value on HPC remote machine
 * // TODO test RAM/cache limits performances and how the algorithm scale with larger N on single machine vs distributed (multi node computing)
 *
 * === NOTES FROM LOCAL BENCHMARKS ===
 *
 * N_ELEMENTS: 500 million
 * MIN value: 0LU
 * MAX value: 1LU << 63 (9,223,372,036,854,775,808)
 *
 * _____________ std::sort: ~13.43 s _____________
 *
 *          codeforces_radix_sort<4>: ~10.26 s
 * [BAD]    codeforces_radix_sort<8>: ~17.52 s
 *
 *          new_radix_sort<8>: ~7.39 s
 *          radix_sort<10>: ~6.74 s
 * [BEST]   radix_sort<11>: ~5.93 s
 *
 * --->         SPEEDUP codeforces_radix_sort<4>:   ~1.31x (-23.6%)
 * --->         SPEEDUP radix_sort<8>:              ~1.82x (-44.97%)
 * --->         SPEEDUP radix_sort<10>:             ~1.99x (-49.8%)
 * ---> [BEST]  SPEEDUP radix_sort<11>:             ~2.26x (-55.8%)
 */