#ifndef MATH_UTIL_H_
#define MATH_UTIL_H_

#include <array>

namespace bdm {

struct Math {
  /// value of pi
  static constexpr float kPi = 3.141592653589793238462643383279502884;

  template <typename T>
  static T Norm(const std::array<T, 3>& a) {
    T norm = a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
    if (norm == 0.0) {
      norm = 1.0;
    }
    return norm;
  }

  template <typename T>
  static std::array<T, 3> Normalize(const std::array<T, 3>& a) {
    auto norm = Norm(a);
    std::array<T, 3> ret;
    ret[0] = a[0] / norm;
    ret[1] = a[1] / norm;
    ret[2] = a[2] / norm;
    return ret;
  }
};

}  // namespace bdm

#endif  // MATH_UTIL_H_
