#include "matrix.h"
#include <gtest/gtest.h>
#include "test_util.h"

namespace bdm {

TEST(MatrixTest, Add) {
  array<float, 3> a = {0.5, 0.7, 1.2};
  array<float, 3> b = {0.6, 1.5, 2.1};
  auto result = Matrix::Add(a, b);

  EXPECT_NEAR(1.1, result[0], abs_error<float>::value);
  EXPECT_NEAR(2.2, result[1], abs_error<float>::value);
  EXPECT_NEAR(3.3, result[2], abs_error<float>::value);
}

TEST(MatrixTest, Subtract) {
  array<float, 3> a = {0.6, 1.5, 2.1};
  array<float, 3> b = {0.5, 0.7, 0.8};
  auto result = Matrix::Subtract(a, b);

  EXPECT_NEAR(0.1, result[0], abs_error<float>::value);
  EXPECT_NEAR(0.8, result[1], abs_error<float>::value);
  EXPECT_NEAR(1.3, result[2], abs_error<float>::value);
}

TEST(MatrixTest, Dot) {
  array<float, 3> a = {0.5, 0.7, 0.8};
  array<float, 3> b = {0.6, 1.5, 2.1};
  float result = Matrix::Dot(a, b);

  EXPECT_NEAR(3.03, result, abs_error<float>::value);
}

}  // namespace bdm
