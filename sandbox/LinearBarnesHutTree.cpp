//
// Created by Matteo Ranzi on 29/11/25.
//

#include <iostream>
#include <iomanip>
#include <vector>

#include "debug/print_debug.h"
#include "utils/random_utils.h"

#define SPACE_DIMENSIONS 3
#define MAX_PARTICLES 1000000
#define NO_CHILDREN (-1)


// TODO set with proper values
#define SIMULATION_ITERATIONS 1
#define THETA 0.5
#define SOFTENING_EPSILON 0.01
#define GRAVITATIONAL_CONSTANT 1.0 // Arbitrary units

struct Particle {
    double position[SPACE_DIMENSIONS]; // center of mass
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

    // Fixme these may be redundant
    bool isLeaf; // Index into particles array (NOT_A_LEAF if internal node)
    size_t particleCount; // Number of particles (1 for leaves, sum for internal nodes)

    enum Children {
        LEFT = 0,
        RIGHT = 1
    };
};

std::vector<Particle> randomParticles(const size_t nParticles) {
    std::vector<Particle> particles(nParticles);
    for (int i = 0; i < nParticles; i++) {
        particles[i] = {random_double(0.0, 1.0), random_double(0.0, 1.0), random_double(0.0, 1.0), 1.0};
    }

    return std::move(particles);
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
void sortParticlesByMortonCode(std::vector<Particle>& particles) {
    for (auto& particle : particles) {
        // assuming positions are normalized [0,1)
        uint32_t discretized_position[SPACE_DIMENSIONS];
        for (int i = 0; i < SPACE_DIMENSIONS; i++) {
            discretized_position[i] = static_cast<uint32_t>(particle.position[i] * ((1<<21) -1));
        }
        particle.morton_code = morton3D_u64(discretized_position);
    }

    PRINT_DEBUG_INFO("Sorting particles...\n");
    std::sort(particles.begin(), particles.end(),
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

void mergeBoundingBoxes(BarnesHutNode& parentNode, const BarnesHutNode& leftChildNode, const BarnesHutNode& rightChildNode) {
    for (int d = 0; d < SPACE_DIMENSIONS; d++) {
        parentNode.boundingBoxMin[d] = std::min(leftChildNode.boundingBoxMin[d], rightChildNode.boundingBoxMin[d]);
        parentNode.boundingBoxMax[d] = std::max(leftChildNode.boundingBoxMax[d], rightChildNode.boundingBoxMax[d]);
    }
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
    node.isLeaf = true;
}

void createInnerNode(std::vector<BarnesHutNode>& barnesHutTree, const size_t parentIndex, const size_t leftChildIndex, const size_t rightChildIndex) {
    barnesHutTree[parentIndex].childrenIndex[BarnesHutNode::LEFT] = leftChildIndex;
    barnesHutTree[parentIndex].childrenIndex[BarnesHutNode::RIGHT] = rightChildIndex;

    // Aggregate data from children
    barnesHutTree[parentIndex].isLeaf = false;
    barnesHutTree[parentIndex].particleCount = barnesHutTree[leftChildIndex].particleCount + barnesHutTree[rightChildIndex].particleCount;
    aggregateMassAndCOM(barnesHutTree[parentIndex], barnesHutTree[leftChildIndex], barnesHutTree[rightChildIndex]);
    mergeBoundingBoxes(barnesHutTree[parentIndex], barnesHutTree[leftChildIndex], barnesHutTree[rightChildIndex]);
}

size_t buildTreeRecursive(std::vector<BarnesHutNode>& barnesHutTree, const std::vector<Particle>& particles, const size_t firstIndex, const size_t lastIndex) {
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
    const size_t leftChildIndex = buildTreeRecursive(barnesHutTree, particles, firstIndex, splitIndex);
    const size_t rightChildIndex = buildTreeRecursive(barnesHutTree, particles, splitIndex + 1, lastIndex);
    // PRINT_DEBUG_INFO("Node %lu split at %lu\t LEFT(%lu, %lu)\tRIGHT(%lu, %lu)\t|\tLeftNode: %lu\tRightNode: %lu\n", currentIndex, splitIndex, firstIndex, splitIndex, splitIndex + 1, lastIndex, leftChildIndex, rightChildIndex);

    createInnerNode(barnesHutTree, currentIndex, leftChildIndex, rightChildIndex);

    return currentIndex;
}
/*
*he resize() method (and passing argument to constructor is equivalent to that) will insert or delete appropriate number of elements to the vector to make it given size (it has optional second argument to specify their value). It will affect the size(), iteration will go over all those elements, push_back will insert after them and you can directly access them using the operator[].
The reserve() method only allocates memory, but leaves it uninitialized. It only affects capacity(), but size() will be unchanged. There is no value for the objects, because nothing is added to the vector. If you then insert the elements, no reallocation will happen, because it was done in advance, but that's the only effect.
*/
std::vector<BarnesHutNode> buildLinearTree(const std::vector<Particle>& particles) {
    std::vector<BarnesHutNode> barnesHutTree;
    // Preallocate memory (with reserve() the size remains 0, but capacity is set)
    barnesHutTree.reserve(2 * particles.size()); // Maximum number of nodes in a binary tree (typically octree uses fewer)

    // Build tree recursively and flatten into array
    buildTreeRecursive(barnesHutTree, particles, 0, particles.size() - 1);

    return barnesHutTree;
}


inline double computeDistance(const Particle& particle, const BarnesHutNode& node, double distance_along_axis[SPACE_DIMENSIONS]) {
    double distSquared = 0.0;
    for (int d = 0; d < SPACE_DIMENSIONS; d++) {
        distance_along_axis[d] = node.COM[d] - particle.position[d];
        distSquared += distance_along_axis[d] * distance_along_axis[d];
    }
    return std::sqrt(distSquared);
}

inline double cellSize(const BarnesHutNode& node) {
    double size = 0.0;
    for (int d = 0; d < SPACE_DIMENSIONS; d++) {
        size = std::max(size, node.boundingBoxMax[d] - node.boundingBoxMin[d]);
    }
    return size;
}

std::vector<double> computeForce(const std::vector<BarnesHutNode>& barnesHutTree, const Particle& particle, const size_t nodeIndex = 0) {
    std::vector<double> forces(SPACE_DIMENSIONS, 0.0);
    std::stack<size_t> indexStack;
    indexStack.push(nodeIndex);

    while (!indexStack.empty()) {
        const size_t currentIndex = indexStack.top();
        indexStack.pop();
        const BarnesHutNode& currentNode = barnesHutTree[currentIndex];

        double distance_along_axis[SPACE_DIMENSIONS];
        const double distance = computeDistance(particle, currentNode, distance_along_axis);

        if (cellSize(currentNode) / distance < THETA || currentNode.isLeaf) {
            // if (currentNode.isLeaf) {
            //     PRINT_DEBUG_INFO("Calculating force from leaf %lu to particle\n", currentIndex);
            // } else {
            //     PRINT_DEBUG_INFO("SINGLE BODY APPROXIMATION: calculating force from inner node %lu to particle\n", currentIndex);
            // }
            // Treat this node as a single body, or it's a leaf node
            // Avoid division by zero (self-interaction) fixme: this solution doesn't properly prevent self-interaction
            if (distance > 0) {
                const double forceMagnitude = GRAVITATIONAL_CONSTANT * (currentNode.mass * particle.mass) / (distance * distance + SOFTENING_EPSILON * SOFTENING_EPSILON);
                for (int d = 0; d < SPACE_DIMENSIONS; d++) {
                    forces[d] += forceMagnitude * distance_along_axis[d] / distance;
                }
            }
        } else {
            // Traverse children
            for (auto childIndex : currentNode.childrenIndex) {
                indexStack.push(childIndex);
            }
        }
    }

    return forces;
}

void updateParticle(Particle& particle, const std::vector<double>& forces, const double timeStep) {

}

void printBarnesHutTree(const std::vector<BarnesHutNode>& tree) {
    std::cout << "\n======================================== BARNES HUT TREE ============================================" << std::endl;
    for (size_t i = 0; i < tree.size(); i++) {
        const BarnesHutNode& node = tree[i];
        std::cout << "[" << std::setw(3) << i
                          << "] mass=" << std::defaultfloat << std::setw(log(MAX_PARTICLES)+1) << node.mass
                          << ",\t|\tCOM=(" << std::fixed << std::setprecision(5) << node.COM[0];
                for (int d = 1; d < SPACE_DIMENSIONS; d++) {
                    std::cout << ", " << std::fixed << std::setprecision(5) << node.COM[d];
                }
                std::cout << ")\t|\tparticleCount=" << std::setw(3) << node.particleCount
                          << "\t|\tchildrenIndex=[LEFT:" << std::setw(log(MAX_PARTICLES*2) +1) << (node.childrenIndex[BarnesHutNode::LEFT] == NO_CHILDREN ? "N/A" : std::to_string(node.childrenIndex[BarnesHutNode::LEFT]))
                            << ", RIGHT:" << std::setw(5) << (node.childrenIndex[BarnesHutNode::RIGHT] == NO_CHILDREN ? "N/A" : std::to_string(node.childrenIndex[BarnesHutNode::RIGHT])) << "]\n";
                          // << ", particleIndex=" << std::setw(5) << (node.particleIndex  == NOT_A_LEAF ? "N/A" : std::to_string(node.particleIndex)) << "\n";
    }
    std::cout << "=====================================================================================================\n" << std::endl;
}


void runSimulation(std::vector<Particle>& particles, const unsigned long iterations) {
    for (unsigned long i = 0; i < iterations; i++) {
        sortParticlesByMortonCode(particles);
        const auto barnesHutTree = buildLinearTree(particles);

        int counter = 1;
        for (auto& particle : particles) {
            // Compute force on particle using Barnes-Hut tree
            updateParticle(particle, computeForce(barnesHutTree, particle), 0.01 /*TODO set a proper time step*/);

            if (counter % 5000 == 0) {
                std::cout << "Processed " << counter << " particles..." << std::endl;
            }
            counter++;
        }

        // printBarnesHutTree(barnesHutTree);
    }
}

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