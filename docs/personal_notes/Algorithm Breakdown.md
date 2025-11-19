# Algorithm Breakdown
Each step involves building the octree and calculating forces. The parallel octree construction and force calculation phases are outlined below:

## 1. Bounding Box Computation
- **Goal**: Determine the minimum spatial extent containing all particles.
- **Linear Method**: Iterate through all particle positions to find the minimum and maximum coordinates along each axis.
- **Parallel Method**: Scatter particle positions across multiple MPI processes, each computing local minima and maxima, then use MPI reduction operations to find global minima and maxima along each axis.
- **Input**: Array of particle positions.
- **Output**: Minimum and maximum corner coordinates (e.g., `min_x`, `min_y`, `min_z`, `max_x`, `max_y`, `max_z`).

# References
- [N-Bodies Algorithm Breakdown - Marc Viva](https://github.com/MarcVivas/N-body/blob/main/BarnesHutExplained.md#algorithm-breakdown)
