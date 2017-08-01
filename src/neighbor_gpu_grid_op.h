#ifndef NEIGHBOR_GPU_GRID_OP_H_
#define NEIGHBOR_GPU_GRID_OP_H_

#include <utility>
#include <vector>

#include "gpu_grid.cuh"

namespace bdm {

/// A class that sets up an uniform grid on GPU to perform operations that require
/// knowledge about neighboring simulation objects
class NeighborGpuGridOp {
 public:
  NeighborGpuGridOp() {}
  virtual ~NeighborGpuGridOp() {}

  template <typename TContainer>
  void Compute(TContainer* cells) const {
    float* d_positions = cells->GetGpuPositions();
    // calculate grid hash
    calculate_hash(d_grid_particle_hash, d_grid_particle_index, 
      d_positions, cells->size());

    // for (size_t i = 0; i < (cells->size()); i++) {
    //   std::cout << "[" << d_positions[i*4] << "," << d_positions[i*4+1] << ", " << d_positions[i*4+2] << ", " << d_positions[i*4+3] << "]" << std::endl;
    // }

    // sort particles based on hash

    // reorder particle arrays into sorted order and find start 
    // and end of each cell (i.e. box)

  }
 private:
  uint *d_grid_particle_hash;   // grid hash value for each particle
  uint *d_grid_particle_index;  // particle index for each particle
};

}  // namespace bdm

#endif  // NEIGHBOR_GPU_GRID_OP_H_
