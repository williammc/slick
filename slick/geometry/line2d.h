// Copyright 2014 The Slick Authors. All rights reserved.
#pragma once
#include <vector>
#include "slick/datatypes.h"

namespace slick {

/// Perpendicular distance from point to line in 2D space
/// @ref: http://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
template <typename Scalar>
inline Scalar perpen_dist_to_line(const Eigen::Matrix<Scalar, 2, 1>& p,
                                  const Eigen::Matrix<Scalar, 3, 1>& leq) {
  const Scalar t = leq[0] * p[0] + leq[1] * p[1] + leq[2];
  return std::fabs(t) / std::sqrt(leq[0] * leq[0] + leq[1] * leq[1]);
}

// Calculate signed normal distance between line and point
template <typename Scalar>
inline Scalar perpen_dist_to_line_signed(
    const Eigen::Matrix < Scalar, 2, 1> & pt,
    const Eigen::Matrix<Scalar, 3, 1>& leq) {
  const Scalar t = leq[0] * p[0] + leq[1] * p[1] + leq[2];
  return t / std::sqrt(leq[0] * leq[0] + leq[1] * leq[1]);
}

/// Perpen. distance from point @p to line (@line_p1, @line_p2)
/// @ref: http://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
template <typename Scalar>
inline Scalar perpen_dist_to_line(const Eigen::Matrix<Scalar, 2, 1>& p,
                                  const Eigen::Matrix<Scalar, 2, 1>& line_p1,
                                  const Eigen::Matrix<Scalar, 2, 1>& line_p2) {
  const Scalar dx = line_p2[0] - line_p1[0];
  const Scalar dy = line_p2[1] - line_p1[1];
  return std::fabs(dy * p[0] - dx * p[1] - line_p1[0] * line_p2[1] +
                   line_p1[1] * line_p2[0]) /
         std::sqrt(dx * dx + dy * dy);
}

/// compute 2D line equation given 2 2D points on the line
template <typename Scalar>
inline Eigen::Matrix<Scalar, 3, 1> line_equation(
    const Eigen::Matrix<Scalar, 2, 1>& p1,
    const Eigen::Matrix<Scalar, 2, 1>& p2) {
  if (p1[0] != p2[0]) {
    const Scalar m = (p2[1] - p1[1]) / (p2[0] - p1[0]);
    const Scalar b = p1[1] - m * p1[0];
    return Eigen::Matrix<Scalar, 3, 1>(m, -1, b);
  } else {
    const Scalar m = (p2[0] - p1[0]) / (p2[1] - p1[1]);
    const Scalar b = p1[0] - m * p1[1];
    return Eigen::Matrix<Scalar, 3, 1>(-1, m, b);
  }
}

/// get a point on given line (@line_eq)
/// each input @index has unique output point
template <typename Scalar>
inline Eigen::Matrix<Scalar, 2, 1> get_a_point(
    const Eigen::Matrix<Scalar, 3, 1>& line_eq, unsigned index) {
  index = index + 1;                                 // zero index cause problem
  Eigen::Matrix<Scalar, 2, 1> apoint(index, index);  // online point
  if (line_eq[1] != 0) {
    apoint[1] = (-line_eq[0] * Scalar(index) - line_eq[2]) / line_eq[1];
  } else {
    apoint[0] = (-line_eq[1] * Scalar(index) - line_eq[2]) / line_eq[0];
  }
  return apoint;
}

/// @return projected point of input point (@p) onto given line (@line_eq)
template <typename Scalar>
inline Eigen::Matrix<Scalar, 2, 1> project_point(
    const Eigen::Matrix<Scalar, 2, 1>& p,
    const Eigen::Matrix<Scalar, 3, 1>& line_eq) {
  const Eigen::Matrix<Scalar, 2, 1> n = line_eq.segment<2>(0).normalized();
  const Eigen::Matrix<Scalar, 2, 1> line_vec(-n[1], n[0]);
  const Eigen::Matrix<Scalar, 2, 1> apoint = get_a_point(line_eq, 1);
  const Scalar s = (p - apoint).transpose() * line_vec;

  return apoint + line_vec * s;
}

/// @return projected point of input point (@p) onto given line (@line_points)
/// @line_points is (x1, y1, x2, y2)
template <typename Scalar>
inline Eigen::Matrix<Scalar, 2, 1> project_point(
    const Eigen::Matrix<Scalar, 2, 1>& p, const Eigen::Matrix<Scalar, 2, 1>& p1,
    const Eigen::Matrix<Scalar, 2, 1>& p2) {
  const Eigen::Matrix<Scalar, 2, 1> line_vec = (p1 - p2).normalized();
  const Scalar s = (p - p1).transpose() * line_vec;

  return p1 + line_vec * s;
}

/// check if input point @p is inside two end points of given line-segment
template <typename Scalar>
inline bool is_point_in(const Eigen::Matrix<Scalar, 2, 1>& p,
                        const Eigen::Matrix<Scalar, 2, 1>& p1,
                        const Eigen::Matrix<Scalar, 2, 1>& p2) {
  const Eigen::Matrix<Scalar, 2, 1> proj_p = project_point(p, p1, p2);
  return (proj_p - p1).dot(proj_p - p2) <= 0;
}

/// 2D line segment (with 2 end points) ========================================
struct Line2d {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  using Point2d = Eigen::Matrix<double, 2, 1>;
  Line2d() = default;
  Line2d(const Point2d& pt1, const Point2d& pt2) : point1(pt1), point2(pt2) {}

  Point2d project(const Point2d& pt) const {
    return project_point(pt, point1, point2);
  }

  double perpendicular_distance(const Point2d& pt) const {
    return perpen_dist_to_line(pt, point1, point2);
  }

  double squared_perpendicular_distance(const Point2d& p) const {
    const double dx = point2[0] - point1[0];
    const double dy = point2[1] - point1[1];
    const double t = (dy * p[0] - dx * p[1] - point1[0] * point2[1] +
                     point1[1] * point2[0]);
    return t*t / (dx * dx + dy * dy);
  }

  double length() const {
    return (point1 - point2).norm();
  }
  Point2d point1, point2;
};
}  // namespace slick
