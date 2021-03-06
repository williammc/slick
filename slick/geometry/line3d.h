// Copyright 2014, 2015 The Slick Authors. All rights reserved.
#pragma once
#include <Eigen/Geometry>
#include "slick/datatypes.h"

namespace slick {

// 3D line segment
template <typename T> struct Line3DBase {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  using PointType = Eigen::Matrix<T, 3, 1>;

  Line3DBase() = default;
  Line3DBase(const PointType &pt1, const PointType &pt2)
      : pt1_(pt1), pt2_(pt2) {}

  PointType project(const PointType &pt) const {
    return project_point(pt, pt1_, pt2_);
  }

  T perpendicular_distance(const PointType &pt) const {
    return std::sqrt(perpendicular_squared_distance(pt));
  }

  T perpendicular_squared_distance(const PointType &pt) const {
    return perpen_sqdist_to_line(pt, pt1_, pt2_);
  }

  PointType &point1() { return pt1_; }
  const PointType &point1() const { return pt1_; }

  PointType &point2() { return pt2_; }
  const PointType &point2() const { return pt2_; }

  T length() const { return (pt2_ - pt1_).norm(); }
  PointType line_vector() const { return pt2_ - pt1_; }

  // useful functionalities ====================================================
  // @ref: http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
  static T perpen_sqdist_to_line(const PointType &p, const PointType &pt1,
                                      const PointType &pt2) {
    const auto v12 = (pt2 - pt1);
    const auto t = v12.cross(pt1 - p);
    return t.squaredNorm() / v12.squaredNorm();
  }

  static PointType project_point(const PointType &pt, const PointType &pt1,
                                 const PointType &pt2) {
    const auto line_vec = (pt2 - pt1).normalized();
    T s = (pt - pt1).dot(line_vec);
    return pt1 + line_vec * s;
  }

  /// check if input point @p is inside two end points of given line-segment
  static bool is_point_in(const PointType &p, const PointType &p1,
                          const PointType &p2) {
    const PointType proj_p = project_point(p, p1, p2);
    return (proj_p - p1).dot(proj_p - p2) <= 0;
  }

protected:
  PointType pt1_, pt2_;
};

using Line3d = Line3DBase<double>;

} // namespace slick
