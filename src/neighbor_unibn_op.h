#ifndef NEIGHBOR_UNIBN_OP_H_
#define NEIGHBOR_UNIBN_OP_H_

#include <fstream>
#include <cmath>
#include <utility>
#include <vector>
#include <chrono>

#include "../third_party/unibn_octree.h"
#include "inline_vector.h"

using std::ofstream;

namespace bdm {

class NeighborUnibnOp {
 public:
  NeighborUnibnOp() {}
  explicit NeighborUnibnOp(double distance) : distance_(distance) {}
  ~NeighborUnibnOp() {}

  template <typename TContainer>
  void Compute(TContainer* cells, const char* filename = nullptr) const {
    ofstream outfile;
    if (filename != nullptr) {
      outfile.open(filename, std::ofstream::out | std::ofstream::app);
    }

    unibn::OctreeParams params;
    params.bucketSize = 16;

    std::vector<unibn::Octant*> octant_mapping(cells->GetAllPositions().size());
    for( auto& element : octant_mapping) { element = nullptr; }

    unibn::Octree<std::array<double, 3> > octree;

    std::chrono::steady_clock::time_point begin_build = std::chrono::steady_clock::now();
    octree.initialize(cells->GetAllPositions(), &octant_mapping, params);
    std::chrono::steady_clock::time_point end_build = std::chrono::steady_clock::now();

    size_t counter = 0;
    for (auto element : octant_mapping) {
      if (element == nullptr) { counter++; }
    }
    std::cout << "uninitized octant mappings " << counter << std::endl;

    if (print_terminal == 1) {
      std::cout << "==================[UniBn]====================" << std::endl;
      std::cout << "Octree build time    = " << std::chrono::duration_cast<std::chrono::milliseconds>(end_build - begin_build).count() << "ms\n";
    }
    if (filename != nullptr) {
      outfile << cells->size() << ",";
      outfile << std::chrono::duration_cast<std::chrono::microseconds>(end_build - begin_build).count() << ",";
    }

// calc neighbors
std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
int avg_num_neighbors = 0;
// holds the indices of found neighbors
std::vector<uint32_t> found_neighbors;
std::vector<float> distances;
double search_radius = sqrt(distance_);
#pragma omp parallel for firstprivate(found_neighbors, search_radius, distances)
    for (size_t i = 0; i < cells->size(); i++) {
      // fixme make param
      // according to roman 50 - 100 micron
      const auto& query = (*cells)[i].GetPosition();
      found_neighbors.clear();

      // calculate neighbors
      // octree.radiusNeighborsCached<unibn::L2Distance<std::array<double, 3> > >(octant_mapping[i], i, query, search_radius, found_neighbors, distances);
      octree.radiusNeighborsRelative<unibn::L2Distance<std::array<double, 3> > >(octant_mapping[i], i, query, distance_, found_neighbors);
      // octree.radiusNeighbors<unibn::L2Distance<std::array<double, 3> > >(query, distance_, found_neighbors);
      avg_num_neighbors += found_neighbors.size();

      // // transform result
      // InlineVector<int, 8> neighbors;
      // neighbors.reserve(found_neighbors.size() - 1);
      // for (size_t j = 0; j < found_neighbors.size(); j++) {
      //   if (found_neighbors[j] != i) {
      //     neighbors.push_back(found_neighbors[j]);
      //   }
      // }
      // (*cells)[i].SetNeighbors(neighbors);
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    if (print_terminal == 1) {
      std::cout << "Neighbor search time = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms\n";
      std::cout << "# of neighbors found = " << (avg_num_neighbors/(cells->size())) << std::endl;
    }
    if (filename != nullptr) {
      outfile << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << ",";
      outfile.close();
    }
  }

 private:
  double distance_ = 3000;
};

}  // namespace bdm

#endif  // NEIGHBOR_UNIBN_OP_H_
