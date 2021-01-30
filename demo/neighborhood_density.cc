#include <omp.h>
#include <cmath>
#include <functional>
#include <iostream>
#include <random>
#include <sstream>

#include "backend.h"
#include "cell.h"
#include "displacement_op.h"
#include "dividing_cell_op.h"
#include "exporter.h"
#include "neighbor_nanoflann_op.h"
#include "neighbor_op.h"
#include "param.h"
#include "resource_manager.h"
#include "scheduler.h"
#include "timing.h"
#include "timing_aggregator.h"
// #include <ittnotify.h>

using bdm::Cell;
using bdm::Exporter;
using bdm::Param;
using bdm::Scalar;
using bdm::Soa;
using bdm::Timing;
using bdm::TimingAggregator;

void execute(size_t num_cells, size_t iterations, double min, double max) {
  // Will be used to obtain a seed for the random number engine
  std::random_device rd;
  std::mt19937 gen(rd());  // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(min, max);
  auto cells = Cell<>::NewEmptySoa();
  cells.reserve(num_cells);
  for (size_t i = 0; i < num_cells; i++) {
    std::array<double, 3> position;
    position[0] = dis(gen);
    position[1] = dis(gen);
    position[2] = dis(gen);
    Cell<Scalar> cell(position);
    cell.SetDiameter(10);
    cell.SetAdherence(0.4);
    cell.SetMass(1.0);
    cell.UpdateVolume();
    cells.push_back(cell);
  }

  bdm::NeighborNanoflannOp neighborhood(20);
  bdm::DividingCellOp biology;
  bdm::DisplacementOp displacement;

  Timing timer;
  auto start = timer.timestamp();

  for (size_t i = 0; i < iterations; i++) {
    neighborhood.Compute(&cells);
    // size_t total_neighbors = cells.CountNeighbors();
    // std::cout << "Neighbors per agent = " << total_neighbors / num_cells
    //           << std::endl;
    biology.Compute(&cells);
    displacement.Compute(&cells);
  }

  auto stop = timer.timestamp();
  std::cout << (stop - start) << std::endl;
}

int main(int args, char** argv) {
  double min_bound = 0;
  double max_bound;
  if (args != 4 ||
      (std::string(argv[1]) == "help" || std::string(argv[1]) == "--help")) {
    // clang-format off
    std::cout << "SYNOPSIS\n"
              << "  ./neighborhood_density <num_cells> "
              << "<iterations> <max_bound> \n"
              << std::endl;
    // clang-format on
    return 1;
  }
  if (Param::kSimulationMaximalDisplacement > 1e-9) {
    std::cout << "ERROR - kSimulationMaximalDisplacement must be set to 0 for "
                 "this benchmark!"
              << std::endl;
    return 1;
  }
  size_t num_cells;
  size_t iterations;
  std::istringstream(std::string(argv[1])) >> num_cells;
  std::istringstream(std::string(argv[2])) >> iterations;
  std::istringstream(std::string(argv[3])) >> max_bound;

  execute(num_cells, iterations, min_bound, max_bound);
  return 0;
}
