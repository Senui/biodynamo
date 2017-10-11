#include "math_util.h"
#include "backend.h"
#include "gtest/gtest.h"
#include "unit/test_util.h"

namespace bdm {

TEST(MathUtilTest, Norm) {
  std::array<double, 3> vector = {1.1, 2.2, 3.3};
  auto result = Math::Norm(vector);

  EXPECT_NEAR(16.94, result, abs_error<double>::value);
}

TEST(MathUtilTest, NormZero) {
  std::array<double, 3> vector = {0, 0, 0};
  auto result = Math::Norm(vector);

  EXPECT_NEAR(1, result, abs_error<double>::value);
}

TEST(MathUtilTest, NormalizeZero) {
  std::array<double, 3> vector = {0, 0, 0};
  auto result = Math::Normalize(vector);

  EXPECT_NEAR(0, result[0], abs_error<double>::value);
  EXPECT_NEAR(0, result[1], abs_error<double>::value);
  EXPECT_NEAR(0, result[2], abs_error<double>::value);
}

TEST(MathUtilTest, Normalize) {
  std::array<double, 3> vector = {1.1, 2.2, 3.3};
  auto result = Math::Normalize(vector);

  EXPECT_NEAR(0.0649350649350649351, result[0], abs_error<double>::value);
  EXPECT_NEAR(0.1298701298701298701, result[1], abs_error<double>::value);
  EXPECT_NEAR(0.1948051948051948052, result[2], abs_error<double>::value);
}

TEST(MathUtilTest, CrossProduct) {
  std::array<double, 3> a = {1.1, 2.2, 3.3};
  std::array<double, 3> b = {5.8, 7.3, 11.87};

  auto&& result = Math::CrossProduct(a, b);
  EXPECT_ARR_NEAR(result, {2.024, 6.083, -4.73});
}

TEST(MathUtilTest, RotAroundAxis) {
  std::array<double, 3> axis = {1.0,1.0,0.0};
  std::array<double, 3> vector = {4,5,6};
  double theta = Math::kPi;

  auto&& result = Math::RotAroundAxis(vector, theta, axis);
  EXPECT_ARR_NEAR(result, {0.5,-0.5, -6});
}

TEST(MathUtilTest, Perp3) {
  std::array<double, 3> vector = {4,5,6};
  double random = 1.1234;

  auto&& result = Math::Perp3(vector, random);
  EXPECT_ARR_NEAR(result, {0.086162044419162448, -0.057217207991650247, -0.0097603562863997576});
}

TEST(MathUtilTest, AngleRadian) {
  std::array<double, 3> a = {1,2,3};
  std::array<double, 3> b = {9,8,7};

  double result = Math::AngleRadian(a, b);
  EXPECT_NEAR(1.5538588453980886, result, abs_error<double>::value);
}

}  // namespace bdm
