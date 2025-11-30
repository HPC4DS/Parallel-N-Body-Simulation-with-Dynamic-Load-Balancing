//
// Created by Matteo Ranzi on 29/11/25.
//

#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <stack>
#include <algorithm>

#include "debug/print_debug.h"
#include "utils/random_utils.h"

#define SPACE_DIMENSIONS 3
#define MAX_PARTICLES 10
#define NO_CHILDREN (-1)


// TODO set with proper values
#define SIMULATION_ITERATIONS 1 //(25*60*1) // 25 FPS for 1 minutes
#define THETA 0.5 // Opening angle threshold
#define SOFTENING_EPSILON 1e-3 // Softening factor to avoid singularities
#define GRAVITATIONAL_CONSTANT 1.0 // Arbitrary units
#define DELTA_TIME 0.01 // Arbitrary units

struct Particle {
    double position[SPACE_DIMENSIONS]; // center of mass
    double velocity[SPACE_DIMENSIONS];
    double mass;
    uint64_t morton_code;

    enum Coordinate {
        X = 0,
        Y = 1,
        Z = 2
    };
};

struct BarnesHutNode {
    double mass;
    double COM[SPACE_DIMENSIONS]; // center of mass
    double boundingBoxMin[SPACE_DIMENSIONS];
    double boundingBoxMax[SPACE_DIMENSIONS];
    size_t childrenIndex[2]; // Indices of child nodes (NO_CHILDREN if leaf)
    const Particle* p_particle;
    double squaredDiagonalLength; // Length of the bounding box diagonal

    // Fixme these may be redundant
    size_t particleCount; // Number of particles (1 for leaves, sum for internal nodes)

    enum Children {
        LEFT = 0,
        RIGHT = 1
    };
};

std::vector<Particle> randomParticles(const size_t nParticles) {
    std::vector<Particle> particles(nParticles);
    for (int i = 0; i < nParticles; i++) {
        particles[i] = {random_double(0.0, 1.0), random_double(0.0, 1.0), random_double(0.0, 1.0),
                        //random_double(-0.01, 0.01), random_double(-0.01, 0.01), random_double(-0.01, 0.01),
                        0,0,0,
                        random_double(1.0, 10.0)};
    }

    return particles;
}
inline uint64_t morton3D_u64(const uint32_t position[SPACE_DIMENSIONS]) {
    auto splitBy2 = [](const uint32_t _v) {
        uint64_t v = _v & 0x1FFFFF; // we only look at the first 21 bits
        v = (v | (v << 32)) & 0x1f00000000ffffULL;
        v = (v | (v << 16)) & 0x1f0000ff0000ffULL;
        v = (v | (v << 8))  & 0x100f00f00f00f00fULL;
        v = (v | (v << 4))  & 0x10c30c30c30c30c3ULL;
        v = (v | (v << 2))  & 0x1249249249249249ULL;
        return v;
    };

    uint64_t morton_code = 0;
    for (int i = 0; i < SPACE_DIMENSIONS; i++) {
        morton_code |= splitBy2(position[i]) << i;
    }
    return morton_code;
}
void computeAndSortParticlesByMortonCode(std::vector<Particle>& particles) {
    for (auto& particle : particles) {
        // assuming positions are normalized [0,1)
        uint32_t discretized_position[SPACE_DIMENSIONS];
        for (int i = 0; i < SPACE_DIMENSIONS; i++) {
            discretized_position[i] = static_cast<uint32_t>(particle.position[i] * ((1<<21) -1));
        }
        particle.morton_code = morton3D_u64(discretized_position);
    }

    PRINT_DEBUG_INFO("Sorting particles...\n");
    std::sort(particles.begin(), particles.end(), //TODO use radix sort
    [](const Particle& a, const Particle& b) {return a.morton_code < b.morton_code;});
    PRINT_DEBUG_INFO("Sorting COMPLETED.\n");

}
size_t findMaxCommonPrefixPosition(const std::vector<Particle>& particles, const size_t firstIndex, const size_t lastIndex) {
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
        return (firstIndex + lastIndex) >> 1; // return the midpoint
    }

    // Find the position of highest differing bit
    const int commonPrefixLength = __builtin_clzll(firstCode ^ lastCode);

    // Use binary search to find where the next bit differs.
    // Specifically, we are looking for the highest object that
    // shares more than commonPrefix bits with the first one.
    size_t splitIndex = firstIndex;
    size_t step = lastIndex - firstIndex;

    do {
        step = (step + 1) >> 1;
        const size_t newSplitIndex = splitIndex + step;
        if (newSplitIndex < lastIndex) {
            const int newSplitPrefixLength = __builtin_clzll(firstCode ^ particles[newSplitIndex].morton_code);
            if (newSplitPrefixLength > commonPrefixLength) {
                splitIndex = newSplitIndex;
            }
        }

    } while (step > 1);

    return splitIndex;
}

inline double computeBoundingBoxSquaredDiagonal(const BarnesHutNode& node) {
    double diagSquared = 0.0;
    for (int d = 0; d < SPACE_DIMENSIONS; d++) {
        const double side = node.boundingBoxMax[d] - node.boundingBoxMin[d];
        diagSquared += side * side;
    }
    return diagSquared;
}
void mergeBoundingBoxes(BarnesHutNode& parentNode, const BarnesHutNode& leftChildNode, const BarnesHutNode& rightChildNode) {
    for (int d = 0; d < SPACE_DIMENSIONS; d++) {
        parentNode.boundingBoxMin[d] = std::min(leftChildNode.boundingBoxMin[d], rightChildNode.boundingBoxMin[d]);
        parentNode.boundingBoxMax[d] = std::max(leftChildNode.boundingBoxMax[d], rightChildNode.boundingBoxMax[d]);
    }
    parentNode.squaredDiagonalLength = computeBoundingBoxSquaredDiagonal(parentNode);
}
void aggregateMassAndCOM(BarnesHutNode& parentNode, const BarnesHutNode& leftChildNode, const BarnesHutNode& rightChildNode) {
    parentNode.mass = leftChildNode.mass + rightChildNode.mass;

    // Weighted average for center of mass
    for (int d = 0; d < SPACE_DIMENSIONS; d++) {
        parentNode.COM[d] = (leftChildNode.mass * leftChildNode.COM[d] + rightChildNode.mass * rightChildNode.COM[d]) / parentNode.mass;
    }
}


void createLeafNode(BarnesHutNode& node, const Particle& particle) {
    node.particleCount = 1;
    node.mass = particle.mass;
    node.childrenIndex[BarnesHutNode::LEFT] = NO_CHILDREN;
    node.childrenIndex[BarnesHutNode::RIGHT] = NO_CHILDREN;
    for (int d = 0; d < SPACE_DIMENSIONS; d++) {
        node.COM[d] = particle.position[d];
        node.boundingBoxMin[d] = particle.position[d];
        node.boundingBoxMax[d] = particle.position[d];
    }
    node.squaredDiagonalLength = 0.0; // Single particle, no size
    node.p_particle = &particle;
}

void createInnerNode(std::vector<BarnesHutNode>& barnesHutTree, const size_t parentIndex, const size_t leftChildIndex, const size_t rightChildIndex) {
    barnesHutTree[parentIndex].childrenIndex[BarnesHutNode::LEFT] = leftChildIndex;
    barnesHutTree[parentIndex].childrenIndex[BarnesHutNode::RIGHT] = rightChildIndex;

    // Aggregate data from children
    barnesHutTree[parentIndex].p_particle = nullptr; // not a leaf
    barnesHutTree[parentIndex].particleCount = barnesHutTree[leftChildIndex].particleCount + barnesHutTree[rightChildIndex].particleCount;
    aggregateMassAndCOM(barnesHutTree[parentIndex], barnesHutTree[leftChildIndex], barnesHutTree[rightChildIndex]);
    mergeBoundingBoxes(barnesHutTree[parentIndex], barnesHutTree[leftChildIndex], barnesHutTree[rightChildIndex]);
}

size_t buildBarnesHutTreeRecursive(std::vector<BarnesHutNode>& barnesHutTree, const std::vector<Particle>& particles, const size_t firstIndex, const size_t lastIndex) {
    // Create a new node
    const size_t currentIndex = barnesHutTree.size();
    barnesHutTree.emplace_back();

    // Base case: leaf node
    if (firstIndex == lastIndex) {
        // PRINT_DEBUG_INFO("Leaf node %lu for particle %lu\n", currentIndex, firstIndex);
        createLeafNode(barnesHutTree[currentIndex], particles[firstIndex]);
        return currentIndex;
    }

    // Internal node: find split position
    const size_t splitIndex = findMaxCommonPrefixPosition(particles, firstIndex, lastIndex);

    // Build inner nodes recursively
    const size_t leftChildIndex = buildBarnesHutTreeRecursive(barnesHutTree, particles, firstIndex, splitIndex);
    const size_t rightChildIndex = buildBarnesHutTreeRecursive(barnesHutTree, particles, splitIndex + 1, lastIndex);
    // PRINT_DEBUG_INFO("Node %lu split at %lu\t LEFT(%lu, %lu)\tRIGHT(%lu, %lu)\t|\tLeftNode: %lu\tRightNode: %lu\n", currentIndex, splitIndex, firstIndex, splitIndex, splitIndex + 1, lastIndex, leftChildIndex, rightChildIndex);

    createInnerNode(barnesHutTree, currentIndex, leftChildIndex, rightChildIndex);

    return currentIndex;
}

void rebuildBarnesHutTree(std::vector<BarnesHutNode>& barnesHutTree, std::vector<Particle>& particles) {
    computeAndSortParticlesByMortonCode(particles);

    barnesHutTree.clear();
    buildBarnesHutTreeRecursive(barnesHutTree, particles, 0, particles.size() - 1);
}

/*
*he resize() method (and passing argument to constructor is equivalent to that) will insert or delete appropriate number of elements to the vector to make it given size (it has optional second argument to specify their value). It will affect the size(), iteration will go over all those elements, push_back will insert after them and you can directly access them using the operator[].
The reserve() method only allocates memory, but leaves it uninitialized. It only affects capacity(), but size() will be unchanged. There is no value for the objects, because nothing is added to the vector. If you then insert the elements, no reallocation will happen, because it was done in advance, but that's the only effect.
*/
std::vector<BarnesHutNode> buildBarnesHutTree(std::vector<Particle>& particles) {
    computeAndSortParticlesByMortonCode(particles);

    std::vector<BarnesHutNode> barnesHutTree;
    // Preallocate memory (with reserve() the size remains 0, but capacity is set)
    barnesHutTree.reserve(2 * particles.size()); // Maximum number of nodes in a binary tree (typically octree uses fewer)

    // Build tree recursively and flatten into array
    buildBarnesHutTreeRecursive(barnesHutTree, particles, 0, particles.size() - 1);

    return barnesHutTree;
}


inline double computeSquaredDistance(const Particle& particle, const BarnesHutNode& node, double distance_along_axis[SPACE_DIMENSIONS]) {
    double distSquared = 0.0;
    for (int d = 0; d < SPACE_DIMENSIONS; d++) {
        distance_along_axis[d] = node.COM[d] - particle.position[d];
        distSquared += distance_along_axis[d] * distance_along_axis[d];
    }
    return distSquared;
}


std::array<double, SPACE_DIMENSIONS> computeAcceleration(const std::vector<BarnesHutNode>& barnesHutTree, const Particle& particle, const size_t nodeIndex = 0) {
    std::array<double, SPACE_DIMENSIONS> acceleration{0.0, 0.0, 0.0};
    std::stack<size_t> indexStack;
    indexStack.push(nodeIndex);

    constexpr double theta2 = THETA * THETA;
    constexpr double eps2 = SOFTENING_EPSILON * SOFTENING_EPSILON;

    while (!indexStack.empty()) {
        const size_t currentIndex = indexStack.top();
        indexStack.pop();
        const BarnesHutNode& currentNode = barnesHutTree[currentIndex];

        double delta[SPACE_DIMENSIONS];
        const double r2 = computeSquaredDistance(particle, currentNode, delta);

        const bool isLeafOtherParticle = (currentNode.p_particle != nullptr && currentNode.p_particle != &particle);
        const bool open = ((r2 > 0.0) && (currentNode.squaredDiagonalLength < theta2 * r2));

        if (open || isLeafOtherParticle) {
            if (r2 > 0.0) {
                // Acceleration magnitude: G * M / (r^2 + epsilon^2)^(3/2)
                const double r2_soft = r2 + eps2;
                const double inv_r3 = 1.0 / (r2_soft * std::sqrt(r2_soft));

                const double accMagnitude = GRAVITATIONAL_CONSTANT * currentNode.mass * inv_r3;
                for (int d = 0; d < SPACE_DIMENSIONS; d++) {
                    // Accumulate acceleration vector components
                    acceleration[d] += accMagnitude * delta[d];
                }
            }
        } else {
            // Traverse children
            for (auto childIndex : currentNode.childrenIndex) {
                if (childIndex != NO_CHILDREN) {
                    indexStack.push(childIndex);
                }
            }
        }
    }

    return acceleration;
}


void initLeapfrogVelocityKick(const std::vector<BarnesHutNode>& barnesHutTree, std::vector<Particle>& particles) {
    // Single half-kick to get velocities to t = dt/2
    constexpr double half_dt = 0.5 * DELTA_TIME;

    for (auto& particle : particles) {
        const auto acceleration = computeAcceleration(barnesHutTree, particle);

        for (int d = 0; d < SPACE_DIMENSIONS; d++) {
            particle.velocity[d] += acceleration[d] * half_dt;
        }
    }

}
void leapfrogVelocityKick(const std::vector<BarnesHutNode>& barnesHutTree, std::vector<Particle>& particles) {
    for (auto& particle : particles) {
        //Fixme: for parallelization, compute force for individual particle vs compute in batches storing value in Particle.ax, Particle.ay, Particle.az ?
        const auto acceleration = computeAcceleration(barnesHutTree, particle);

        for (int d = 0; d < SPACE_DIMENSIONS; d++) {
            particle.velocity[d] += acceleration[d] * DELTA_TIME;
        }
    }
}

void leapfrogPositionDrift(std::vector<Particle>& particles) {
    for (auto& particle : particles) {
        for (int d = 0; d < SPACE_DIMENSIONS; d++) {
            particle.position[d] += particle.velocity[d] * DELTA_TIME;
        }
    }

    //TODO optimize (all over the source code) by removing for loop (?) and using SIMD operations (improve readability creating x, y , z, vx, vy, vz variables?)
    /*
     * What's the difference in performance in SIMD instruction grouping the operations above in a loop for all particles?
     * E.g.:
            #pragma omp parallel for simd
            for i from 0 to N-1:
                particles[i].x += dt * particles[i].vx
                particles[i].y += dt * particles[i].vy
                particles[i].z += dt * particles[i].vz
     */
}



void printBarnesHutTree(const std::vector<BarnesHutNode>& tree) {
    std::cout << "\n======================================== BARNES HUT TREE ============================================" << std::endl;
    for (size_t i = 0; i < tree.size(); i++) {
        const BarnesHutNode& node = tree[i];
        std::cout << "[" << std::setw(  static_cast<int>(std::log10(tree.size()))+1) << i
                          << "] mass: " << std::defaultfloat << std::setw(static_cast<int>(std::log10(MAX_PARTICLES))+1) << node.mass
                          << "\t|\tCOM: (" << std::fixed << std::setprecision(5) << node.COM[0];
                for (int d = 1; d < SPACE_DIMENSIONS; d++) {
                    std::cout << ", " << std::fixed << std::setprecision(5) << node.COM[d];
                }
                std::cout << ") \t|\tparticleCount:" << std::setw(3) << node.particleCount
                          << "\t|\tchildrenIndex:[LEFT:" << std::setw(static_cast<int>(std::log10(MAX_PARTICLES * 2)) +2) << (node.childrenIndex[BarnesHutNode::LEFT] == NO_CHILDREN ? "N/A" : std::to_string(node.childrenIndex[BarnesHutNode::LEFT]))
                            << ", RIGHT:" << std::setw(static_cast<int>(std::log10(MAX_PARTICLES * 2)) +2) << (node.childrenIndex[BarnesHutNode::RIGHT] == NO_CHILDREN ? "N/A" : std::to_string(node.childrenIndex[BarnesHutNode::RIGHT])) << "]\n";
                          // << ", particleIndex=" << std::setw(5) << (node.particleIndex  == NOT_A_LEAF ? "N/A" : std::to_string(node.particleIndex)) << "\n";
    }
    std::cout << "=====================================================================================================\n" << std::endl;
}

//fixme: complete simulation loop with energy conservation checks
//fixme: currently simulation is based on fixed number of iterations, independently from time step (is it ok?)
void runSimulation(std::vector<Particle>& particles, const unsigned long iterations) {
    auto barnesHutTree = buildBarnesHutTree(particles);

    initLeapfrogVelocityKick(barnesHutTree, particles);
    // From now on, positions are at integer times: t = 0, dt, 2dt, ...
    // And velocities are at half-integer times: t = dt/2, 3dt/2, 5dt/2, ...

    //TODO compute initial energy to validate conservation through the simulation

    // ===================== MAIN SIMULATION LOOP =====================
    double time = 0.0;
    for (unsigned long i = 0; i < iterations; i++) {
        // Drift: update positions x^(n) -> x^(n+1)
        leapfrogPositionDrift(particles);

        rebuildBarnesHutTree(barnesHutTree, particles);

        // Kick: update velocities v^(n+1/2) -> v^(n+3/2)
        leapfrogVelocityKick(barnesHutTree, particles);

        time += DELTA_TIME;
    }
}

/*void old_runSimulation(std::vector<Particle>& particles, const unsigned long iterations) {
    std::cout << "\n=====================================" << std::endl;
    std::cout << "STARTING SIMULATION" << std::endl;
    std::cout << "=====================================\n" << std::endl;

    computeAndSortParticlesByMortonCode(particles);
    auto barnesHutTree = buildBarnesHutTree(particles);

    initializeLeapfrogIntegration(barnesHutTree, particles);

    for (unsigned long i = 1; i <= iterations; i++) {
        computeAndSortParticlesByMortonCode(particles);
        barnesHutTree = buildBarnesHutTree(particles);

        int counter = 1;
        int percentage = 0;
        for (auto& particle : particles) {
            // Compute force on particle using Barnes-Hut tree
            //leapfrogStep(particle, DELTA_TIME);

            //TODO move code into function to update progress bar of the current iteration
            constexpr int percentageResolution = 100;
            if (counter / static_cast<double>(MAX_PARTICLES) * 100.0 >= percentage * percentageResolution) {
                std::cout << std::setw(static_cast<int>(log10(iterations)) +1) << i << "/" << iterations << ") Processed " << std::setw(3) << percentage * percentageResolution << "% (" << std::setw(static_cast<int>(std::log10(MAX_PARTICLES))+1) << counter << "/" << MAX_PARTICLES << ")" << std::endl;
                percentage++;
            }
            counter++;
        }
        //TODO move code into function to update progress bar of the entire simulation
        std::cout << "#" << i << " iteration finished: " << static_cast<double>(i) / static_cast<double>(iterations) * 100.0 << "%" << std::endl;
        std::cout << "----------------------------------------" << std::endl;

        printBarnesHutTree(barnesHutTree);
    }
    std::cout << "\n=====================================" << std::endl;
    std::cout << "END OF SIMULATION" << std::endl;
    std::cout << "=====================================\n" << std::endl;
}*/

void printParticles(const std::vector<Particle>& particles) {
    std::cout << "\n================ PARTICLES =====================" << std::endl;
    for (size_t i = 0; i < particles.size(); i++) {
        std::bitset<64> mc_bits(particles[i].morton_code);
        std::cout << "[" << std::setw(3) << i << "] " << mc_bits << std::endl;
    }
    std::cout << "=====================================\n" << std::endl;
}

int main() {
    constexpr size_t nParticles = MAX_PARTICLES;

    std::vector<Particle> particles = randomParticles(nParticles);
    runSimulation(particles, SIMULATION_ITERATIONS);


    //printParticles(particles);
    return 0;
}

//==============================================

/*
 * EFFICIENCY SUGGESTIONS
* TODO Use radix sort for Morton codes to avoid O(n log n) sort cost each step.
* TODO Cache and reuse the tree when possible; consider incremental updates or rebuild only every k steps if accuracy allows.
* TODO Precompute and reuse cellSize or store node size/half-diagonal to avoid repeated bounding-box scans.
* fixme Parallelize per-particle acceleration with OpenMP; ensure tree is read-only during traversal.
* Avoid heap allocations in hot loops (use fixed-size arrays, reserve stacks if needed).
* TODO Consider a true octree (8 children) for tighter opening criterion in 3D; if staying binary BVH, refine the size metric (use bounding box diagonal).
* TODO Use SoA layout (separate arrays for positions/velocities/masses) for better SIMD and cache locality.
* TODO Replace repeated computeDistance with inlined arithmetic using precomputed COM to reduce function overhead.
 */