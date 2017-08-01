#ifndef GPU_GRID_H_
#define GPU_GRID_H_

#include "gpu_grid.cuh"
#include "gpu_kernels.cuh"
#include "helper_cuda.h"
#include <stdio.h>

extern "C"
{

namespace bdm {

// void copyArrayToDevice(void *device, const void *host, int offset, int size)
// {
//   checkCudaErrors(cudaMemcpy((char *) device + offset, host, size, cudaMemcpyHostToDevice));
// }

// void ParticleSystem::setArray(const float *data, int start, int count) {
//   copyArrayToDevice(m_dVel, data, start*4*sizeof(float), count*4*sizeof(float));
// }

//Round a / b to nearest higher integer value
uint iDivUp(uint a, uint b) {
  return (a % b != 0) ? (a / b + 1) : (a / b);
}

void computeGridSize(uint n, uint blockSize, uint &numBlocks, uint &numThreads) {
  numThreads = min(blockSize, n);
  numBlocks = iDivUp(n, numThreads);
}

void calculate_hash(uint  *gridParticleHash,
                  uint  *gridParticleIndex,
                  float *pos,
                  int    numParticles) {
  uint numThreads, numBlocks;
  computeGridSize(numParticles, 64, numBlocks, numThreads);

  // hello_world<<<5,5>>>();

  printf("Launching calculate_hash_d<<<%d, %d>>>...\n", numBlocks, numThreads);

  // execute the kernel
  calculate_hash_d<<< numBlocks, numThreads >>>(gridParticleHash,
                                         gridParticleIndex,
                                         (float4 *) pos,
                                         numParticles);

  // check if kernel invocation generated an error
  getLastCudaError("Kernel execution failed");
  cudaDeviceReset();
}
}  // namespace bdm
}

#endif  // GPU_GRID_H_
