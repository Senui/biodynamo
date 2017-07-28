#include "math_util.h"
#include "backend.h"
#include "gtest/gtest.h"
#include "test_util.h"

namespace bdm {

TEST(MathUtilTest, Norm) {
  std::array<float, 3> vector = {1.1, 2.2, 3.3};
  auto result = Math::Norm(vector);

  EXPECT_FLOAT_EQ(16.94, result);
}

TEST(MathUtilTest, NormZero) {
  std::array<float, 3> vector = {0, 0, 0};
  auto result = Math::Norm(vector);

  EXPECT_NEAR(1, result, abs_error<float>::value);
}

TEST(MathUtilTest, NormalizeZero) {
  std::array<float, 3> vector = {0, 0, 0};
  auto result = Math::Normalize(vector);

  EXPECT_NEAR(0, result[0], abs_error<float>::value);
  EXPECT_NEAR(0, result[1], abs_error<float>::value);
  EXPECT_NEAR(0, result[2], abs_error<float>::value);
}

TEST(MathUtilTest, Normalize) {
  std::array<float, 3> vector = {1.1, 2.2, 3.3};
  auto result = Math::Normalize(vector);

  EXPECT_NEAR(0.0649350649350649351, result[0], abs_error<float>::value);
  EXPECT_NEAR(0.1298701298701298701, result[1], abs_error<float>::value);
  EXPECT_NEAR(0.1948051948051948052, result[2], abs_error<float>::value);
}

}  // namespace bdm
