// Copyright 2014 The Slick Authors. All rights reserved.
#pragma once
#include <memory>
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

  PointType point() const {
    const auto n = normal();
    const T d = pln_eq_[3];
    PointType point;
    if (n[0] != 0.0) {
      T x = -(n[1] + n[2] + d) / n[0];
      point << x, 1, 1;
    } else if (n[1] != 0.0) {
      T y = -(n[0] + n[2] + d) / n[1];
      point << 1, y, 1;
    } else {
      T z = -(n[0] + n[1] + d) / n[2];
      point << 1, 1, z;
    }
    return point;
  }

  /// @ref: http://www.9math.com/book/projection-point-plane
  PointType project(const PointType &pt) const {
    const T t = pln_eq_.head<3>().dot(pt) + pln_eq_[3];
    const T s = pln_eq_.head<3>().squaredNorm();
    const T t1 = t / s;
    const T x = pt[0] - pln_eq_[0]*t1;
    const T y = pt[1] - pln_eq_[1]*t1;
    const T z = pt[2] - pln_eq_[2]*t1;
    return PointType(x, y, z);
  }

  T perpendicular_distance(const PointType &pt) const {
    return std::sqrt(perpendicular_sqdistance(pt));
  }

  T perpendicular_sqdistance(const PointType &pt) const {
    const T t = pln_eq_.head<3>().dot(pt) + pln_eq_[3];
    return t * t / normal().squaredNorm();
  }

  std::unique_ptr<PointType> intersect(Line3DBase<T> const &ln) const {
    PointType intersection;
    const auto n = normal();
    const auto p = point();
    PointType lv = ln.line_vector();

    if (std::fabs(lv.dot(n)) < 1e-12)
      return std::unique_ptr<PointType>();

    PointType lp = ln.point1();
    T t = (p.dot(n) - n.dot(lp)) / (n.dot(lv));
    intersection = lp + t * lv;
    return std::unique_ptr<PointType>(new PointType(intersection));
  }

  /// @ref: http://mathworld.wolfram.com/Plane-PlaneIntersection.html
  std::unique_ptr<Line3DBase<T>> intersect(const Plane3DBase &other_plane) const {
    Line3DBase<T> line;
    const auto n = normal();
    const auto line_vector = n.cross(other_plane.normal());
    if (line_vector.norm() < 1.e-9) // two planes are parallel
      return std::unique_ptr<Line3DBase<T>>();
    if (pln_eq_[3] == 0 && other_plane.pln_eq_[3] == 0) {
      line.point1() << 0, 0, 0;
      line.point2() = line.point1() + line_vector;
      return std::unique_ptr<Line3DBase<T>>(new Line3DBase<T>(line));
    }

    auto solve_onaxis = [&](int axis) {
      // equations
      // a1*x + b1*y + c1*z + d1 = 0;
      // a2*x + b2*y + c2*z + d2 = 0;
      // suppose point of intersection line is on axis=1 plane
      int i1 = (axis + 1) % 3;
      int i2 = (axis + 2) % 3;
      Eigen::Matrix<T, 2, 2> A;
      A << n[i1], n[i2], other_plane.normal()[i1], other_plane.normal()[i2];
      Eigen::Matrix<T, 2, 1> b;
      b << -n[axis] - pln_eq_[3], -other_plane.normal()[axis] - other_plane.pln_eq_[3];
      const auto v = A.inverse() * b;
      line.point1()[axis] = 1.0;
      line.point1()[i1] = v[0];
      line.point1()[i2] = v[1];
      line.point2() = line.point1() + line_vector;
    };

    if (std::fabs(line_vector.dot(Eigen::Matrix<T, 3, 1>(0.0, 0.0, 1.0))) > 1.e-7) {
      solve_onaxis(2);
    } else if (std::fabs(line_vector.dot(Eigen::Matrix<T, 3, 1>(0.0, 1.0, 0.0))) > 1.e-7) {
      solve_onaxis(1);
    } else {
      solve_onaxis(0);
    }
    return std::unique_ptr<Line3DBase<T>>(new Line3DBase<T>(line));
  }

  Eigen::Matrix<T, 3, 1> normal() const { return pln_eq_.head<3>(); }

  Eigen::Matrix<T, 4, 1> equation() const { return pln_eq_; }

  /// Normalize normal
  void normalize() {
    const T t = pln_eq_.head<3>().norm();
    pln_eq_ /= t;
  }

protected:
  Eigen::Matrix<T, 4, 1> pln_eq_;
};

using Plane3d = Plane3DBase<double>;

} // namespace slick