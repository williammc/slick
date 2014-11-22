// Copyright 2012, 2013, 2014 The Look3D Authors. All rights reserved.
#ifndef LOOK3D_GEOMETRY_LINE2D_H_
#define LOOK3D_GEOMETRY_LINE2D_H_
#include <vector>
#include <Eigen/Dense>
#include <opencv2/core/core.hpp>
#include "math/common.h"
#include "geometry/look3d_geometry_api.h"

namespace look3d {

// representation of 2D line/line segment
struct LOOK3D_GEOMETRY_API Line2d {
  // Default Constructor.
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Line2d() {}

  // Constructor with 2 points
  Line2d(const Eigen::Vector2d& pt1, const Eigen::Vector2d& pt2, bool linevec = false);

  // Constructor for common line-equation (a,b,c)
  Line2d(DefaultScalarType a, DefaultScalarType b, DefaultScalarType c);

  // best-fitting line through point cloud (using SVD)
  explicit Line2d(const std::vector<Eigen::Vector2d> &pts);

  // value consider zero (e.x: identical points) if less than line precision.
  bool is_zero(DefaultScalarType d) const{
    if (d < epsilon)
      return true;
    return false;
  }

  // Shifts the line by given vector
  Line2d shift_line_by_vector(const Eigen::Vector2d &vector) const {
    return Line2d(point1 + vector, point2 + vector);
  }

  Line2d scale_line(DefaultScalarType scale) {
    return Line2d(point1 * scale, point2 * scale);
  }

  bool operator == (const Line2d &line) {
    if ((point1 == line.point1) && (point2 == line.point2))
      return true;
    else
      return false;
  }

  // Return length of the line segment
  DefaultScalarType length() const {
    return ((point1 - point2).norm());
  }

  // Intersection of two lines - return true if intersection exists
  bool IntersectLines(const Line2d& line, Eigen::Vector2d& intersection) const;


  // check if two line segments intersect (does not return intersection point)
  bool CheckIntersectLines(const Line2d& line) const;


  // Returns orthogonal line through point pt
  inline Line2d perpendicular_line(const Eigen::Vector2d &pt) const {
    Line2d back(pt, pt + normal);
    return back;
  }

  // Calculate absolute normal distance between line and point
  DefaultScalarType perpendicular_distance(const Eigen::Vector2d& pt) const {
    return (fabs((a*pt.x()+b*pt.y()+c) / sqrt(a*a+b*b)));
  }

  // Calculate signed normal distance between line and point
  DefaultScalarType perpendicular_distance_signed(const Eigen::Vector2d& pt) const {
    return ((a*pt.x()+b*pt.y()+c) / sqrt(a*a+b*b));
  }

  // Calculates squared normal distance between line and point (speedup)
  DefaultScalarType squared_perpendicular_distance(const Eigen::Vector2d& pt) const {
    DefaultScalarType a = a*pt.x()+b*pt.y()+c;
    return((a*a)/ (a*a+b*b));
  }

  // Project a 2D point onto the line
  Eigen::Vector2d project_pt(const Eigen::Vector2d &point) const {
    DefaultScalarType s = (point - point1).transpose() * line_vector;
    return point1 + line_vector * s;
  }

  // returns the quadrant (1..4) of the line vector (pt2 - pt1)
  int GetQuadrant() const;

  // returns the enclosed angle (in radians!!!)
  DefaultScalarType enclosed_angle(const Line2d &line) const {
    return std::acos(std::fabs(line_vector.transpose() * line.line_vector));
  }

  // Check if a point lies on this line subject to line precision (epsilon).
  bool contains_point(const Eigen::Vector2d& pt) const {
    const DefaultScalarType d = perpendicular_distance(pt);
    if (d > epsilon )
      return false;
    if (!is_bounded) {
      return true;  // line is unbounded
    } else {
      // projected length
      DefaultScalarType length_ = (pt - point1).transpose() * line_vector;
      if (length_ > 0.0)
        return (length_ < length());
      else
        return false;
    }
  }

  // Check if this line align with another line
  bool is_aligned(const Line2d& another_line,
                  const DefaultScalarType angle_threshold,
                  const DefaultScalarType dist_threshold) const {
    const DefaultScalarType cos_angle = line_vector.dot(another_line.line_vector);
    if (std::fabs(cos_angle) > std::cos(M_PI/2-angle_threshold)) {
      const DefaultScalarType len1 = perpendicular_distance(another_line.point1);
      if (len1 < dist_threshold) {
        const DefaultScalarType len2 = perpendicular_distance(another_line.point2);
        if (len2 < dist_threshold)
          return true;
      }
    }
    return false;
  }

  // compare function for sorting (descending order - longest lines first)
  static bool compare_length(const Line2d &l1, const Line2d &l2) {
    return l1.length() > l2.length();
  }

  // Calculate a, b, and c out of two points on the line
  void CalcParams(void);
  void CalcPointsFromLineEquation(int offset = 100);

  void calc_line_vector() {
    line_vector = point2 - point1;
    line_vector.normalize();
    normal = Eigen::Vector2d(line_vector.y(), -line_vector.x());
  }

  static const DefaultScalarType epsilon;  // Define line precision.
  Eigen::Vector2d point1;
  Eigen::Vector2d point2;
  Eigen::Vector2d normal;
  bool is_bounded;
  DefaultScalarType a, b, c;
  Eigen::Vector2d line_vector;
  DefaultScalarType score;
};
}  // namespace look3d
#endif  // LOOK3D_GEOMETRY_LINE2D_H_
