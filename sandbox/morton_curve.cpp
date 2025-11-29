//
// Created by Matteo Ranzi on 22/11/25.
//

#include <iostream>
#include <bitset>
#include <algorithm>


// #if defined(__i386__) || defined(__x86_64__)
/* Header file for SIMD intrinsics (MMX, SSE, AVX, ...)*/
// #include <immintrin.h>
/* example of functions: _pdep_u64() */
// #endif

#include "debug/print_debug.h"
#include "utils/random_utils.h"

void visualize_morton_steps();
void visualize_morton_positions();

struct Particle {
    double x, y, z;
    double mass;
    uint64_t morton_code;
};

struct Node {
    double mass;
    double com[3]; // center of mass
    double size;
    uint32_t start;
    uint32_t count;
    uint32_t child[8];
};

inline uint32_t morton3D_u32(const uint32_t x, const uint32_t y, const uint32_t z) {
    auto splitBy2 = [](uint32_t v) {
        v = (v | (v << 16)) & 0x030000FF;
        v = (v | (v << 8))  & 0x0300F00F;
        v = (v | (v << 4))  & 0x030C30C3;
        v = (v | (v << 2))  & 0x09249249;
        return v;
    };

    return (splitBy2(z) << 2) | (splitBy2(y) << 1) | splitBy2(x);
}

/////////////////////////////////////////////////////////
// Morton / Z-order encoder (3D) - 21 bits per coord
// Produces a 63-bit Morton code in uint64_t
//
// The helper splitBy2() "spreads" the bits of a 21-bit value
// so that between each original bit there are two zero bits:
//   b20  -> b20__ (positions spaced by 3)
//   b19  -> b19__ ...
//
// Note: masks are chosen to isolate and move chunks safely.
// TODO: improve performance using AVX instructions
/////////////////////////////////////////////////////////
inline uint64_t morton3D_u64(const uint32_t x, const uint32_t y, const uint32_t z) {
    auto splitBy2 = [](const uint32_t _v) {
        uint64_t v = _v & 0x1FFFFF; // we only look at the first 21 bits
        v = (v | (v << 32)) & 0x1f00000000ffffULL;
        v = (v | (v << 16)) & 0x1f0000ff0000ffULL;
        v = (v | (v << 8))  & 0x100f00f00f00f00fULL;
        v = (v | (v << 4))  & 0x10c30c30c30c30c3ULL;
        v = (v | (v << 2))  & 0x1249249249249249ULL;
        return v;
    };

    return (splitBy2(z) << 2) | (splitBy2(y) << 1) | splitBy2(x);
}

std::vector<Particle> randomParticles(const size_t nParticles) {
    std::vector<Particle> particles(nParticles);
    for (int i = 0; i < nParticles; i++) {
        particles[i] = {random_double(0.0, 1.0), random_double(0.0, 1.0), random_double(0.0, 1.0), 1.0};
    }

    return std::move(particles);
}

std::vector<Particle> generate_gaussian_clusters(size_t N, int numClusters) {
    std::vector<Particle> out;
    out.reserve(N);

    std::mt19937_64 rng(12345); // Fixed seed for reproducibility
    std::uniform_real_distribution<double> uni(0.0, 1.0);

    // Random cluster centers
    struct C { double x,y,z; };
    std::vector<C> centers(numClusters);
    for(int i=0;i<numClusters;i++)
        centers[i] = { uni(rng), uni(rng), uni(rng) };

    // Extremely small spread (the key!)
    double sigma = 1e-5;  // cluster radius ~1e-6
    std::normal_distribution<double> gauss(0.0, sigma);

    for(size_t i=0;i<N;i++) {
        int k = i % numClusters;
        double x = centers[k].x + gauss(rng);
        double y = centers[k].y + gauss(rng);
        double z = centers[k].z + gauss(rng);

        // Keep inside [0,1]
        x = std::clamp(x,0.0,1.0);
        y = std::clamp(y,0.0,1.0);
        z = std::clamp(z,0.0,1.0);

        out.push_back({x,y,z,1.0});
    }
    return out;
}

void sortParticlesByMortonCode(std::vector<Particle>& particles) {
    // struct MortonParticle {uint32_t code; Particle particle;};
    // std::vector<MortonParticle> morton_particles;
    // morton_particles.reserve(particles.size());

    int index = 0;
    for (auto& particle : particles) {
        // assuming positions are normalized [0,1)
        const auto x = static_cast<uint32_t>(particle.x * ((1<<21) -1)); // expand number from [0,1) to [0, 2^21 -1)
        const auto y = static_cast<uint32_t>(particle.y * ((1<<21) -1));
        const auto z = static_cast<uint32_t>(particle.z * ((1<<21) -1));
        particle.morton_code = morton3D_u64(x, y, z);
        // morton_particles.push_back({morton3D_u32(x, y, z), particle});
        if (index++ % 10000 == 0) {
            // std::cout << "Processed particles: " << index << std::endl;
        }
    }

    std::cout << "Sorting particles..." << std::endl;
    std::sort(particles.begin(), particles.end(),
    [](const Particle& a, const Particle& b) {return a.morton_code < b.morton_code;});
    std::cout << "Sorting COMPLETED!" << std::endl;
}

int findMaxCommonPrefixPosition(const std::vector<Particle>& particles, const int firstIndex, const int lastIndex) {
    // Find the split position that maximizes the common prefix
    // This determines how to partition the range [first, last)

    if (firstIndex == lastIndex) {
        return firstIndex;
    }

    // Find the longest common prefix between first and last codes
    const uint64_t firstCode = particles[firstIndex].morton_code;
    const uint64_t lastCode = particles[lastIndex].morton_code;

    if (firstCode == lastCode) {
        // All particles have the same Morton code in this range
        return (firstIndex + lastIndex) / 2;
    }

    // Find the position of highest differing bit
    const int commonPrefixLength = __builtin_clzll(firstCode ^ lastCode) / 3; // divide by 3 for 3D Morton code

    // Binary search to find the split point where prefix changes
    int splitIndex = firstIndex; // initial guess
    int splitPrefixLength = commonPrefixLength;

    for (int i = firstIndex + 1; i < lastIndex; i++) {
        const int currentPrefixLength = __builtin_clzll(firstCode ^ particles[i].morton_code) / 3;
        if (currentPrefixLength > splitPrefixLength) {
            splitIndex = i;
            splitPrefixLength = currentPrefixLength;
        }
    }

    return splitIndex;
}


int main() {
    constexpr size_t nParticles = 10;
    std::vector<Particle> particles = randomParticles(nParticles);

    sortParticlesByMortonCode(particles);

    std::cout << "\n=====================================\n" << std::endl;
    for (const auto& particle : particles) {
        std::bitset<64> mc_bits(particle.morton_code);
        std::cout << mc_bits << std::endl;
    }
    std::cout << "\n=====================================\n" << std::endl;

    int firstIndex = 0;
    int lastIndex = nParticles - 1;
    std::cout << "Split [ " << firstIndex << ", " << lastIndex << "] at index: "
              << findMaxCommonPrefixPosition(particles, firstIndex, lastIndex) << std::endl;
}

int old_main () {
    constexpr int nParticles = 10000000;
    // std::vector<Particle> particles = generate_gaussian_clusters(nParticles, 25);
    std::vector<Particle> particles(nParticles);
    // particles.reserve(nParticles);
    for (int i = 0; i < nParticles; i++) {
        particles.push_back({random_double(0.0, 1.0), random_double(0.0, 1.0), random_double(0.0, 1.0), 1.0});
        if (i % 10000 == 0) {
            // std::cout << "Generated particle: " << i << std::endl;
        }
    }

    // std::cout << "***\nBEFORE morton sorting\n***" << std::endl;
    // for (auto particle: particles) {
    //     std::cout << "Particle: (" << particle.x << ", " << particle.y << ", " << particle.z << ")\n";
    // }

    sortParticlesByMortonCode(particles);

    // std::cout << "\n***\nAFTER morton sorting\n***" << std::endl;
    // for (auto particle: particles) {
    //     std::cout << "Particle: (" << particle.x << ", " << particle.y << ", " << particle.z << ")\n";
    // }

    std::cout << "\n=====================================\n" << std::endl;


    bool flag = true;
    int parent_level_count[64] = {0};
    for (int i = 0; i < particles.size(); i++) {
        long left_prefix = -1;
        long right_prefix = -1;

        if (i > 0)                      left_prefix = __builtin_clzll(particles[i-1].morton_code ^ particles[i].morton_code);
        if (i < particles.size() - 1)   right_prefix = __builtin_clzll(particles[i].morton_code ^ particles[i+1].morton_code);

        const long parent_prefix = std::max(left_prefix, right_prefix);
        parent_level_count[parent_prefix]++;
        if (flag && parent_prefix == 23) {
            flag = false;
            std::bitset<64> mb1, mb2;

            if (left_prefix > right_prefix) {
                mb1 = particles[i-1].morton_code;
                mb2 = particles[i].morton_code;
                PRINT_DEBUG_INFO("M1\nx = %.30f\ny = %.30f\nz = %0.30f\n", particles[i-1].x, particles[i-1].y, particles[i-1].z);
                PRINT_DEBUG_INFO("M2\nx = %.30f\ny = %.30f\nz = %0.30f\n", particles[i].x, particles[i].y, particles[i].z);
            } else {
                mb1 = particles[i].morton_code;
                mb2 = particles[i+1].morton_code;
                PRINT_DEBUG_INFO("M1\nx = %.30f\ny = %.30f\nz = %0.30f\n", particles[i].x, particles[i].y, particles[i].z);
                PRINT_DEBUG_INFO("M2\nx = %.30f\ny = %.30f\nz = %0.30f\n", particles[i+1].x, particles[i+1].y, particles[i+1].z);
            }

            std:: cout << "\n" << mb1 << "\n" << mb2 << "\n";

        }
        std::bitset<64>mc_bits(particles[i].morton_code);
        // std::cout << "Particle " << i << ": morton_code = " << mc_bits << ", parent_prefix = " << parent_prefix << std::endl;

        if (i % 10000 == 0) {
            // std::cout << "Parented particle: " << i << std::endl;
        }
    }

    int totalCount = 0;
    for (int i = 0; i < 64; i++) {
        std::cout << "Level " << i << ": " << parent_level_count[i] << " particles" << std::endl;
        totalCount += parent_level_count[i];
    }
    std::cout << "Total particles counted: " << totalCount << std::endl;

    return 0;
}

void visualize_morton_steps() {
    std::bitset<32> mask_bits;

    uint32_t v = 0b1100101;
    std::bitset<32>v_bits(v);

    //==========================================
    std::cout << "    v_bits: " << v_bits << std::endl;
    std::cout << "v_bits<<16: " << (v_bits<<16) << std::endl;
    v = (v | (v << 16)) & 0b00000011000000000000000011111111; //0x030000FF
    v_bits = v;
    mask_bits = 0x030000FF;
    std::cout << " mask_bits: " << mask_bits << std::endl;
    //==========================================

    //==========================================
    std::cout << "\n    v_bits: " << v_bits << std::endl;
    std::cout << "v_bits<<8 : " << (v_bits<<8) << std::endl;
    v = (v | (v << 8))  & 0b00000011000000001111000000001111; // 0x0300F00F
    v_bits = v;
    mask_bits = 0x0300F00F;
    std::cout << " mask_bits: " << mask_bits << std::endl;
    //==========================================

    //==========================================
    std::cout << "\n    v_bits: " << v_bits << std::endl;
    std::cout << "v_bits<<4 : " << (v_bits<<4) << std::endl;
    v = (v | (v << 4))  & 0b00000011000011000011000011000011; //0x030C30C3
    v_bits = v;
    mask_bits = 0x030C30C3;
    std::cout << " mask_bits: " << mask_bits << std::endl;
    //==========================================

    //==========================================
    std::cout << "\n    v_bits: " << v_bits << std::endl;
    std::cout << "v_bits<<2 : " << (v_bits<<2) << std::endl;
    v = (v | (v << 2))  & 0b00001001001001001001001001001001; //0x09249249
    v_bits = v;
    mask_bits = 0x09249249;
    std::cout << " mask_bits: " << mask_bits << std::endl;
    //==========================================

    std::cout << "\n    v_bits: " << v_bits << std::endl;
}

void visualize_morton_positions() {
    constexpr int max = 0b1111111111;
    constexpr int min = 0b1111111100;

    for (int _z = max; _z >= min; _z--) {
        for (int y = max; y >= min; y--) {
            for (int _x = max; _x >= min; _x--) {
                const uint32_t morton_code = morton3D_u32(_x, y, _z);
                const std::bitset<32>mc_bits(morton_code);
                std::cout << "morton3D(" << _x << ", " << y << ", " << _z << ") = " << morton_code << "\t|\t " << mc_bits << std::endl;
            }
        }
    }
}