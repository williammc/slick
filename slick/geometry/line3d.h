// Copyright 2014 The Slick Authors. All rights reserved.
#pragma once
#include <Eigen/Dense>
#include "slick/datatypes.h"
#include "slick/slick_api.h"

namespace slick {
struct SLICK_API Line3d {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Line3d() {}

  Line3d(Eigen::Matrix<SlickScalar, 3, 1> pt1, Eigen::Matrix<SlickScalar, 3, 1> pt2, bool linevec = false)
    : point1(pt1), point2(pt2) {
    if (linevec) {
      is_bounded = false;
      point1 = pt1;
      point2 = pt1 + pt2.normalized();
    } else {
      is_bounded = true;
    }
    line_vector = (point2 - point1).normalized();
  }

  const Eigen::Matrix<SlickScalar, 3, 1>& linevec() const {
    return line_vector;
  }

  const SlickScalar length() const {
    return line_vector.norm();
  }

  const bool isbounded() const {
    return is_bounded;
  }

  void set_linevec(Eigen::Matrix<SlickScalar, 3, 1> linevec) {
    line_vector = linevec.normalized();
  }


  bool contains_point(const Eigen::Matrix<SlickScalar, 3, 1>& pt) const {
    const SlickScalar d = perpendicular_distance(pt);
    if (d > epsilon )
      return false;
    if (!is_bounded) {
      return true;  // line is unbounded
    } else {
      // projected length
      SlickScalar length_ = (pt - point1).transpose() * line_vector;
      if (length_ > 0.0)
        return (length_ < length());
      else
        return false;
    }
  }

  // Perpendicular distance from a point.
  SlickScalar perpendicular_distance(const Eigen::Matrix<SlickScalar, 3, 1>& point) const {
    const Eigen::Matrix<SlickScalar, 3, 1> p = project_pt(point);
    return (p-point).norm();
  }

  // projects a 3D point onto the line
  Eigen::Matrix<SlickScalar, 3, 1> project_pt(const Eigen::Matrix<SlickScalar, 3, 1> &point) const {
    SlickScalar s = (point - point1).transpose() * line_vector;
    return point1 + line_vector * s;
  }

  static const SlickScalar epsilon;  // Define line precision.

  Eigen::Matrix<SlickScalar, 3, 1> point1;
  Eigen::Matrix<SlickScalar, 3, 1> point2;
  Eigen::Matrix<SlickScalar, 3, 1> line_vector;

  bool is_bounded;
};
}  // namespace slick
