#ifndef SPATIAL_ORGANIZATION_BOUND_H_
#define SPATIAL_ORGANIZATION_BOUND_H_

#include <cmath>
#include <utility>
#include "spatial_organization/point.h"

namespace bdm {
namespace spatial_organization {

using std::pair;
using std::make_pair;

/// Class 'Bound' represents rectangular Bound of the node of the tree
class Bound {
 public:
  /// far left Bottom, Near Right Top
  Point far_left_bottom_point_, near_right_top_point_;

  Bound() {
    far_left_bottom_point_ = Point(0, 0, 0);
    near_right_top_point_ = Point(1, 1, 1);
  }

  Bound(float x1, float y1, float z1, float x2, float y2, float z2) {
    far_left_bottom_point_ = Point(x1, y1, z1);
    near_right_top_point_ = Point(x2, y2, z2);
  }

  Bound(Point p1, Point p2)
      : far_left_bottom_point_(p1), near_right_top_point_(p2) {}

  /// Check if 'x' is between 'a' and 'b' on the line
  bool IsBetween(float x, float a, float b) const {
    float min = fmin(a, b);
    float max = fmax(a, b);
    return (x >= min && x <= max);
  }

  /// Calculate distance between two line segments on the line.
  /// This method assumed that segments are not overlaped.
  float DistanceBetweenSegments(float x, float y, float a, float b) const {
    float min_xy = fmin(x, y);
    float max_xy = fmax(x, y);
    float min_ab = fmin(a, b);
    float max_ab = fmax(a, b);

    if (min_xy >= max_ab) {
      return min_xy - max_ab;
    }
    return min_ab - max_xy;
  }

  /// Calculate squared distance between two boundaries in 3-d space.
  float SquaredDistance(Bound const &b) const {
    bool is_overlap_x;
    bool is_overlap_y;
    bool is_overlap_z;

    float bx[2][2] = {{Far(), Near()}, {b.Far(), b.Near()}};
    float by[2][2] = {{Left(), Right()}, {b.Left(), b.Right()}};
    float bz[2][2] = {{Bottom(), Top()}, {b.Bottom(), b.Top()}};

    // check axis_ that have overlaped projections
    is_overlap_x = (IsBetween(bx[0][0], bx[1][0], bx[1][1])) ||
                   (IsBetween(bx[0][1], bx[1][0], bx[1][1])) ||
                   (IsBetween(bx[1][0], bx[0][0], bx[0][1])) ||
                   (IsBetween(bx[1][1], bx[0][0], bx[0][1]));
    is_overlap_y = (IsBetween(by[0][0], by[1][0], by[1][1])) ||
                   (IsBetween(by[0][1], by[1][0], by[1][1])) ||
                   (IsBetween(by[1][0], by[0][0], by[0][1])) ||
                   (IsBetween(by[1][1], by[0][0], by[0][1]));
    is_overlap_z = (IsBetween(bz[0][0], bz[1][0], bz[1][1])) ||
                   (IsBetween(bz[0][1], bz[1][0], bz[1][1])) ||
                   (IsBetween(bz[1][0], bz[0][0], bz[0][1])) ||
                   (IsBetween(bz[1][1], bz[0][0], bz[0][1]));

    float dx = 0, dy = 0, dz = 0;

    // calculate distance only if there is no overlaping
    if (!is_overlap_x) {
      dx = DistanceBetweenSegments(bx[0][0], bx[0][1], bx[1][0], bx[1][1]);
    }

    if (!is_overlap_y) {
      dy = DistanceBetweenSegments(by[0][0], by[0][1], by[1][0], by[1][1]);
    }

    if (!is_overlap_z) {
      dz = DistanceBetweenSegments(bz[0][0], bz[0][1], bz[1][0], bz[1][1]);
    }
    return dx * dx + dy * dy + dz * dz;
  }

  Point Center() const {
    return (far_left_bottom_point_ + near_right_top_point_) * 0.5;
  }

  float Near() const { return near_right_top_point_.x_; }

  float Far() const { return far_left_bottom_point_.x_; }

  float Left() const { return far_left_bottom_point_.y_; }

  float Right() const { return near_right_top_point_.y_; }

  float Top() const { return near_right_top_point_.z_; }

  float Bottom() const { return far_left_bottom_point_.z_; }

  float Length() const {
    return near_right_top_point_.x_ - far_left_bottom_point_.x_;
  }

  float Width() const {
    return near_right_top_point_.y_ - far_left_bottom_point_.y_;
  }

  float Height() const {
    return near_right_top_point_.z_ - far_left_bottom_point_.z_;
  }

  float HalfSurfaceArea() const {
    return Width() * Length() + Height() * Length() + Width() * Height();
  }
};
}  // namespace spatial_organization
}  // namespace bdm
#endif  // SPATIAL_ORGANIZATION_BOUND_H_
