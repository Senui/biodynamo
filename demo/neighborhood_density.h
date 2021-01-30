#ifndef DEMO_NEIGHBORHOOD_DENSITY_H_
#define DEMO_NEIGHBORHOOD_DENSITY_H_

#include "biodynamo.h"

#include "timing.h"

namespace bdm {

template <typename Backend>
struct CompileTimeParam : public DefaultCompileTimeParam<Backend> {
  using BiologyModules = Variant<GrowDivide>;
};

inline int Simulate(int argc, const char** argv) {
  InitializeBiodynamo(argc, argv);

  size_t num_cells;
  size_t iterations;
  double max_bound;

  if (argc != 4) {
    std::cout
        << "Usage: ./neighborhood_density <population> <iterations> <max_bound>"
        << std::endl;
    return 1;
  }

  if (Param::simulation_max_displacement_ > 1e-9) {
    std::cout
        << "ERROR - Param::simulation_max_displacement_ must be set to 0 for "
           "this benchmark! It's currently set to "
        << Param::simulation_max_displacement_ << std::endl;
    return 1;
  }

  std::istringstream(std::string(argv[1])) >> num_cells;
  std::istringstream(std::string(argv[2])) >> iterations;
  std::istringstream(std::string(argv[3])) >> max_bound;

  auto construct = [](const std::array<float, 3>& position) {
    Cell cell(position);
    cell.SetDiameter(10);
    cell.SetAdherence(0.4);
    cell.SetMass(1.0);
    cell.AddBiologyModule(GrowDivide());
    return cell;
  };
  ModelInitializer::CreateCellsRandom(0, max_bound, num_cells, construct);

  Scheduler<> scheduler;
  Timing t;
  auto start = t.Timestamp();
  scheduler.Simulate(iterations);
  auto stop = t.Timestamp();
  std::cout << stop - start << std::endl;
  return 0;
}

}  // namespace bdm

#endif  // DEMO_NEIGHBORHOOD_DENSITY_H_
