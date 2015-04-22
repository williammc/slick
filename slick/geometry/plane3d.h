// Copyright 2014 The Slick Authors. All rights reserved.
#pragma once
#include <vector>
#include <Eigen/Dense>
#include "slick/geometry/line3d.h"
#include "slick/util/mestimator.h"
#include "slick/datatypes.h"

namespace slick {

template <typename T> struct Plane3DBase {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  using PointType = Eigen::Matrix<T, 3, 1>;
  Plane3DBase() = default;

  Plane3DBase(const PointType &p1, const PointType &p2, const PointType &p3) {
    pln_eq_.head<3>() = ((p1 - p2).cross(p1 - p3)).normalized();
    pln_eq_[3] = -(normal().dot(p1));
    normalize();
  }

  Plane3DBase(const PointType &point, const PointType &normal) {
    pln_eq_.head<3>() = normal;
    pln_eq_[3] = -(normal.dot(point));
    normalize();
  }

  // Given plane equation ax + by + cz + d = 0
  Plane3DBase(const Vec4 &v4_plane) : pln_eq_(v4_plane) {}

  void point() const {
    const auto normal = normal();
    const T d = pln_eq_[3];
    PointType point;
    if (normal[0] != 0.0) {
      T x = -(normal[1] + normal[2] + d) / normal[0];
      point << x, 1, 1;
    } else if (normal[1] != 0.0) {
      T y = -(normal[0] + normal[2] + d) / normal[1];
      point << 1, y, 1;
    } else {
      T z = -(normal[0] + normal[1] + d) / normal[2];
      point << 1, 1, z;
    }
    return point;
  }

  T dist_to_point(const PointType &pt) const {
    const auto point = point();
    return std::fabs(normal().dot(point - pt));
  }

  std::pair<bool, PointType> intersect_line(Line3DBase<T> const &ln) const {
    PointType intersection;
    const auto normal = normal();
    const auto point = point();
    PointType lv = ln.line_vector();

    if (std::fabs(lv.dot(normal)) < 1e-12)
      return false;

    PointType lp = ln.point1();
    T t = (point.dot(normal) - normal.dot(lp)) / (normal.dot(lv));
    intersection = lp + t * lv;
    return std::make_pair(true, intersection);
  }

  std::pair<bool, Line3DBase<T>>
  intersect_plane(const Plane3DBase &other_plane) const {
    Line3DBase<T> line;
    const auto normal = normal();
    const auto line_vector = normal.cross(other_plane.normal());
    if (line_vector.norm() < 1.e-9) // two planes are parallel
      return std::make_pair(false, Line3DBase<T>());
    if (pln_eq_[3] == 0 && other_plane.pln_eq_[3] == 0) {
      line.point1() << 0, 0, 0;
      line.point2() = line.point1() + line_vector;
      return std::make_pair(true, line);
    }

    // equations
    // a1*x + b1*y + c1*z + d1 = 0;
    // a2*x + b2*y + c2*z + d2 = 0;
    // suppose point of intersection line is on Z=1 plane
    Eigen::Matrix<T, 2, 2> A;
    A << normal[0], normal[1], other_plane.normal[0], other_plane.normal[1];
    Eigen::Matrix<T, 2, 1> b;
    b << -normal[2] - d, -other_plane.normal[2] - other_plane.d;
    line.point1().segment<2>(0) = A.inverse() * b;
    line.point1()[2] = 1.0;
    line.point2() = line.point1() + line_vector;
    std::make_pair(true, line);
  }

  PointType project_point(const PointType &pt) const {
    Line3DBase<T> ln(pt, pt + normal);
    PointType intp;
    intersect_line(ln, intp);
    return intp;
  }

  decltype(pln_eq_.head<3>()) normal() const { return pln_eq_.head<3>(); }

  /// Normalize normal
  void noramlize() {
    const T t = pln_eq_.head<3>().norm();
    pln_eq_ /= t;
  }

protected:
  Eigen::Matrix<T, 4, 1> pln_eq_;
};

using Plane3d = Plane3DBase<double>;

} // namespace slick