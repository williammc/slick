// Copyright 2014 The Slick Authors. All rights reserved.
#pragma once
#include <vector>
#include "slick/datatypes.h"

namespace slick {
// 2D line segment
template <typename Scalar> struct Line2DBase {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  using PointType = Eigen::Matrix<Scalar, 2, 1>;
  Line2DBase() = default;
  Line2DBase(const PointType &pt1, const PointType &pt2)
      : pt1_(pt1), pt2_(pt2) {}

  PointType project(const PointType &pt) const {
    return project_point(pt, pt1_, pt2_);
  }

  Scalar perpendicular_distance(const PointType &pt) const {
    return perpen_dist_to_line(pt, pt1_, pt2_);
  }

  Scalar perpendicular_sqdistance(const PointType &pt) const {
    return perpen_sqdist_to_line(pt, pt1_, pt2_);
  }

  Scalar squared_perpendicular_distance(const PointType &p) const {
    const Scalar dx = pt2_[0] - pt1_[0];
    const Scalar dy = pt2_[1] - pt1_[1];
    const Scalar t =
        (dy * p[0] - dx * p[1] - pt1_[0] * pt2_[1] + pt1_[1] * pt2_[0]);
    return t * t / (dx * dx + dy * dy);
  }

  PointType &point1() { return pt1_; }
  const PointType &point1() const { return pt1_; }

  PointType &point2() { return pt2_; }
  const PointType &point2() const { return pt2_; }

  Scalar length() { return (pt2_ - pt1_).norm(); }
  PointType line_vector() { return pt2_ - pt1_; }

  // Useful functions ==========================================================
  /// Perpendicular distance from point to line in 2D space
  /// @ref: http://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
  static Scalar perpen_dist_to_line(const PointType &p,
                                    const Eigen::Matrix<Scalar, 3, 1> &leq) {
    const Scalar t = leq[0] * p[0] + leq[1] * p[1] + leq[2];
    return std::fabs(t) / std::sqrt(leq[0] * leq[0] + leq[1] * leq[1]);
  }

  // Calculate signed normal distance between line and point
  static Scalar
  perpen_dist_to_line_signed(const Eigen::Matrix<Scalar, 2, 1> &pt,
                             const Eigen::Matrix<Scalar, 3, 1> &leq) {
    const Scalar t = leq[0] * p[0] + leq[1] * p[1] + leq[2];
    return t / std::sqrt(leq[0] * leq[0] + leq[1] * leq[1]);
  }

  /// Perpen. distance from point @p to line (@line_p1, @line_p2)
  /// @ref: http://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
  static Scalar perpen_dist_to_line(const PointType &p,
                                    const PointType &line_p1,
                                    const PointType &line_p2) {
    const Scalar dx = line_p2[0] - line_p1[0];
    const Scalar dy = line_p2[1] - line_p1[1];
    return std::fabs(dy * p[0] - dx * p[1] + line_p1[1] * line_p2[0] -
                     line_p1[0] * line_p2[1]) /
           std::sqrt(dx * dx + dy * dy);
  }

  /// Perpen. distance from point @p to line (@line_p1, @line_p2)
  /// @ref: http://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
  static Scalar perpen_sqdist_to_line(const PointType &p,
                                      const PointType &line_p1,
                                      const PointType &line_p2) {
    const Scalar dx = line_p2[0] - line_p1[0];
    const Scalar dy = line_p2[1] - line_p1[1];
    const Scalar t = dy * p[0] - dx * p[1] + line_p1[1] * line_p2[0] -
                     line_p1[0] * line_p2[1];
    return t * t / (dx * dx + dy * dy);
  }

  /// compute 2D line equation given 2 2D points on the line
  static Eigen::Matrix<Scalar, 3, 1>
  line_equation(const Eigen::Matrix<Scalar, 2, 1> &p1,
                const Eigen::Matrix<Scalar, 2, 1> &p2) {
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
  static Eigen::Matrix<Scalar, 2, 1>
  get_a_point(const Eigen::Matrix<Scalar, 3, 1> &line_eq, unsigned index) {
    index = index + 1; // zero index cause problem
    Eigen::Matrix<Scalar, 2, 1> apoint(index, index); // online point
    if (line_eq[1] != 0) {
      apoint[1] = (-line_eq[0] * Scalar(index) - line_eq[2]) / line_eq[1];
    } else {
      apoint[0] = (-line_eq[1] * Scalar(index) - line_eq[2]) / line_eq[0];
    }
    return apoint;
  }

  /// @return projected point of input point (@p) onto given line (@line_eq)
  static PointType project_point(const PointType &p,
                                 const Eigen::Matrix<Scalar, 3, 1> &line_eq) {
    const PointType n = line_eq.segment<2>(0).normalized();
    const PointType line_vec(-n[1], n[0]);
    const PointType apoint = get_a_point(line_eq, 1);
    const Scalar s = (p - apoint).transpose() * line_vec;

    return apoint + line_vec * s;
  }

  /// @return projected point of input point (@p) onto given line (@line_points)
  /// @line_points is (x1, y1, x2, y2)
  static PointType project_point(const PointType &p, const PointType &p1,
                                 const PointType &p2) {
    const PointType line_vec = (p1 - p2).normalized();
    const Scalar s = (p - p1).dot(line_vec);

    return p1 + line_vec * s;
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

using Line2d = Line2DBase<double>;
} // namespace slick
