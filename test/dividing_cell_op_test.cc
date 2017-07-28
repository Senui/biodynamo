#include "dividing_cell_op.h"
#include "cell.h"
#include "gtest/gtest.h"
#include "test_util.h"
#include "transactional_vector.h"

namespace bdm {
namespace dividing_cell_op_test_internal {

template <typename T>
void RunTest(T* cells) {
  cells->push_back(Cell<>(41.0));
  cells->push_back(Cell<>(19.0));

  float volume_mother = (*cells)[0].GetVolume();

  DividingCellOp op;
  op.Compute(cells);

  EXPECT_EQ(3u, cells->size());
  EXPECT_FLOAT_EQ(19.005288996600001, (*cells)[1].GetDiameter());
  EXPECT_FLOAT_EQ(3594.3640018287319, (*cells)[1].GetVolume());

  // cell got divided so it must be smaller than before
  // more detailed division test can be found in `cell_test.h`
  EXPECT_GT(41, (*cells)[0].GetDiameter());
  EXPECT_GT(41, (*cells)[2].GetDiameter());
  // volume of two daughter cells must be equal to volume of the mother
  EXPECT_FLOAT_EQ(volume_mother, (*cells)[0].GetVolume() + (*cells)[2].GetVolume());
}

TEST(DividingCellOpTest, ComputeAos) {
  TransactionalVector<Cell<Scalar>> cells;
  RunTest(&cells);
}

TEST(DividingCellOpTest, ComputeSoa) {
  auto cells = Cell<>::NewEmptySoa();
  RunTest(&cells);
}

}  // namespace dividing_cell_op_test_internal
}  // namespace bdm
