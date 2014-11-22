// Copyright 2014 The Slick Authors. All rights reserved.
#pragma once

#include <vector>
#include <Eigen/Dense>
#include "slick/geometry/line3d.h"
#include "slick/datatypes.h"
#include "slick/slick_api.h"

namespace slick {

struct SLICK_API Plane3d {
  Plane3d() {
    point = Eigen::Matrix<DefaultScalarType, 3, 1>::Zero();
    normal = Eigen::Matrix<DefaultScalarType, 3, 1>::Zero();
    d = 0;
  }

  Plane3d(const Eigen::Matrix<DefaultScalarType, 3, 1> &p1,
          const Eigen::Matrix<DefaultScalarType, 3, 1> &p2,
          const Eigen::Matrix<DefaultScalarType, 3, 1> &p3) {
    point = p1;
    normal = ((p1 - p2).cross(p1 - p3)).normalized();
    d = - (normal.dot(point));
  }

  Plane3d(const Eigen::Matrix<DefaultScalarType, 3, 1> &point,
          const Eigen::Matrix<DefaultScalarType, 3, 1> &normal)
    : point(point), normal(normal.normalized()) {
    d = - (normal.dot(point));
  }

  // Given plane equation ax + by + cz + d = 0
  Plane3d(const Eigen::Matrix<DefaultScalarType, 4, 1> &v4_plane) {
    SetEquation(v4_plane);
  }

  void SetEquation(const Eigen::Matrix<DefaultScalarType, 4, 1> &v4_plane) {
    normal = v4_plane.segment<3>(0).normalized();
    d = v4_plane[3]/v4_plane.segment<3>(0).norm();
    if (d < 0) {  // normalize show that normal point toward origin
      d = -d;
      normal = - normal;
    }
    if (normal[0] != 0.0) {
      DefaultScalarType x = -(normal[1] + normal[2] + d)/normal[0];
      point << x, 1, 1;
    } else if(normal[1] != 0.0) {
      DefaultScalarType y = -(normal[0] + normal[2] + d)/normal[1];
      point << 1, y, 1;
    } else {
      DefaultScalarType z = -(normal[0] + normal[1] + d)/normal[2];
      point << 1, 1, z;
    }
  }

  static std::pair<std::vector<unsigned int>, DefaultScalarType> fit_plane(
      const std::vector<Eigen::Matrix<DefaultScalarType, 3, 1> >& points, Plane3d &plane);

  // fit dominant plane containing: ZERO, 1, or 2 anchor points
  // @return (inlier points, average squared error)
  static std::pair<std::vector<unsigned int>, DefaultScalarType> fit_plane_advanced(
    const std::vector<Eigen::Matrix<DefaultScalarType, 3, 1> >& points, Plane3d& plane,
    Eigen::Matrix<DefaultScalarType, 3, 1>* anchor_point1 = NULL, Eigen::Matrix<DefaultScalarType, 3, 1>* anchor_point2 = NULL);

  DefaultScalarType dist_to_point(const Eigen::Matrix<DefaultScalarType, 3, 1> &pt) const {
    return fabs(normal.dot(point - pt));
  }

  bool intersect_line(slick::Line3d const&ln,
                         Eigen::Matrix<DefaultScalarType, 3, 1> &intersection) const {
    Eigen::Matrix<DefaultScalarType, 3, 1> lv = ln.line_vector;

    if (fabs(lv.dot(normal)) < 1e-12)
      return false;

    Eigen::Matrix<DefaultScalarType, 3, 1> lp = ln.point1;
    DefaultScalarType t =  (point.dot(normal) - normal.dot(lp)) / (normal.dot(lv));
    intersection = lp + t * lv;
    return true;
  }

  bool intersect_plane(const Plane3d& other_plane, Line3d& line) const {
    line.line_vector = normal.cross(other_plane.normal);
    if (line.line_vector.norm() < 1.e-9)  // two planes are parallel
      return false;
    if (d == 0 && other_plane.d == 0) {
      line.point1 << 0, 0, 0;
      line.point2 = line.point1 + line.line_vector;
      return true;
    }

    // equations
    // a1*x + b1*y + c1*z + d1 = 0;
    // a2*x + b2*y + c2*z + d2 = 0;
    // suppose point of intersection line is on Z=1 plane
    Eigen::Matrix2d A;
    A << normal[0], normal[1],
        other_plane.normal[0], other_plane.normal[1];
    Eigen::Vector2d b;
    b << -normal[2] - d,
        -other_plane.normal[2] - other_plane.d;
    line.point1.segment<2>(0) = A.inverse()*b;
    line.point1[2] = 1.0;
    line.point2 = line.point1 + line.line_vector;
    return true;
  }

  Eigen::Matrix<DefaultScalarType, 3, 1> project_point(const Eigen::Matrix<DefaultScalarType, 3, 1> &pt) const {
    Line3d ln(pt, pt + normal);
    Eigen::Matrix<DefaultScalarType, 3, 1> intp;
    intersect_line(ln, intp);
    return intp;
  }

  bool is_same(const Plane3d& another_plane,
               DefaultScalarType angle_threshold = 5, DefaultScalarType dist_threshold = 0.01) const {
    if (is_parallel(another_plane, angle_threshold)) {
      Eigen::Matrix<DefaultScalarType, 3, 1> p = project_point(another_plane.point);
      if(((p-another_plane.point).norm() < dist_threshold) || dist_threshold < 0)
        return true;
    }
    return false;
  }

  bool is_parallel(const Plane3d& another_plane,
                   DefaultScalarType angle_threshold = 5) const{
    if (std::cos(std::fabs(normal.dot(another_plane.normal))) >
        std::cos(M_PI*angle_threshold/180))  {
        return true;
    }
    return false;
  }

  static const DefaultScalarType epsilon;  // plane precision
  Eigen::Matrix<DefaultScalarType, 3, 1> point;
  Eigen::Matrix<DefaultScalarType, 3, 1> normal;
  DefaultScalarType d;
};
}  // namespace slick
