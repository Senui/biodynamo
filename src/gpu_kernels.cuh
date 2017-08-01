#ifndef GPU_KERNELS_H_
#define GPU_KERNELS_H_

#define USE_TEX 0

#if USE_TEX
#define FETCH(t, i) tex1Dfetch(t##Tex, i)
#else
#define FETCH(t, i) t[i]
#endif

#include "vector_types.h"
#include <stdio.h>
typedef unsigned int uint;

struct SimParams
{
    float3 colliderPos;
    float  colliderRadius;

    float3 gravity;
    float globalDamping;
    float particleRadius;

    uint3 gridSize;
    uint numCells;
    float3 worldOrigin;
    float3 cellSize;

    uint numBodies;
    uint maxParticlesPerCell;

    float spring;
    float damping;
    float shear;
    float attraction;
    float boundaryDamping;
};

__constant__ SimParams params;

namespace bdm {

// calculate position in uniform grid
__device__ int3 calcGridPos(float3 p) {
  int3 gridPos;
  gridPos.x = floor((p.x - params.worldOrigin.x) / params.cellSize.x);
  gridPos.y = floor((p.y - params.worldOrigin.y) / params.cellSize.y);
  gridPos.z = floor((p.z - params.worldOrigin.z) / params.cellSize.z);
  // printf("GridPos = %d,%d,%d\n", gridPos.x, gridPos.y, gridPos.z);
  return gridPos;
}

// calculate address in grid from position (clamping to edges)
__device__ uint calcGridHash(int3 gridPos) {
  gridPos.x = gridPos.x & (params.gridSize.x-1);  // wrap grid, assumes size is power of 2
  gridPos.y = gridPos.y & (params.gridSize.y-1);
  gridPos.z = gridPos.z & (params.gridSize.z-1);
  return __umul24(__umul24(gridPos.z, params.gridSize.y), params.gridSize.x) + __umul24(gridPos.y, params.gridSize.x) + gridPos.x;
}

__global__
void hello_world() {
  if (threadIdx.x == 0) {
    printf("Hello from block %d\n", blockIdx.x);
  }
}

// calculate grid hash value for each particle
__global__
void calculate_hash_d(uint   *gridParticleHash,  // output
               uint   *gridParticleIndex, // output
               float4 *pos,               // input: positions
               uint    numParticles) {
  uint index = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;

  // pos is not on the device yet! So this kernel will segfault atm

  // if (index == 0) {
  //   for (size_t i = 0; i < numParticles; i++) {
  //     printf("[%f, %f, %f]\n", pos[i*4], pos[i*4+1], pos[i*4+2], pos[i*4+3]);
  //   }
  // }

  // if (index >= numParticles) return;

  // needs to be volatile?
  // volatile float4 p = pos[index];
  // float4 p = pos[index];

  // get address in grid
  // int3 gridPos = calcGridPos(make_float3(p.x, p.y, p.z));

  // uint hash = calcGridHash(gridPos);

  // // store grid hash and particle index
  // gridParticleHash[index] = hash;
  // gridParticleIndex[index] = index;
}

}  // namespace bdm

#endif  // GPU_KERNELS_H_
