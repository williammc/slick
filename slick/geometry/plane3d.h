// Copyright 2014 The Slick Authors. All rights reserved.
#pragma once
#include <vector>
#include <Eigen/Dense>
#include "slick/geometry/line3d.h"
#include "slick/util/mestimator.h"
#include "slick/datatypes.h"

namespace slick {

struct Plane3d {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Plane3d() {
    point = Vec3::Zero();
    normal = Vec3::Zero();
    d = 0;
  }

  Plane3d(const Vec3 &p1,
          const Vec3 &p2,
          const Vec3 &p3) {
    point = p1;
    normal = ((p1 - p2).cross(p1 - p3)).normalized();
    d = - (normal.dot(point));
  }

  Plane3d(const Vec3 &point,
          const Vec3 &normal)
    : point(point), normal(normal.normalized()) {
    d = - (normal.dot(point));
  }

  // Given plane equation ax + by + cz + d = 0
  Plane3d(const Vec4& v4_plane) {
    SetEquation(v4_plane);
  }

  void SetEquation(const Vec4 &v4_plane) {
    normal = v4_plane.segment<3>(0).normalized();
    d = v4_plane[3]/v4_plane.segment<3>(0).norm();
    if (d < 0) {  // normalize show that normal point toward origin
      d = -d;
      normal = - normal;
    }
    if (normal[0] != 0.0) {
      SlickScalar x = -(normal[1] + normal[2] + d)/normal[0];
      point << x, 1, 1;
    } else if(normal[1] != 0.0) {
      SlickScalar y = -(normal[0] + normal[2] + d)/normal[1];
      point << 1, y, 1;
    } else {
      SlickScalar z = -(normal[0] + normal[1] + d)/normal[2];
      point << 1, 1, z;
    }
  }

  static std::pair<std::vector<unsigned int>, SlickScalar> FitPlane(
      const std::vector<Vec3 >& points, Plane3d &plane);

  // fit dominant plane containing: ZERO, 1, or 2 anchor points
  // @return (inlier points, average squared error)
  static std::pair<std::vector<unsigned int>, SlickScalar> FitPlaneAdvance(
    const std::vector<Vec3 >& points, Plane3d& plane,
    Vec3* anchor_point1 = nullptr, Vec3* anchor_point2 = nullptr);

  SlickScalar dist_to_point(const Vec3 &pt) const {
    return fabs(normal.dot(point - pt));
  }

  bool intersect_line(slick::Line3d const&ln,
                         Vec3 &intersection) const {
    Vec3 lv = ln.line_vector;

    if (fabs(lv.dot(normal)) < 1e-12)
      return false;

    Vec3 lp = ln.point1;
    SlickScalar t =  (point.dot(normal) - normal.dot(lp)) / (normal.dot(lv));
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

  Vec3 project_point(const Vec3 &pt) const {
    Line3d ln(pt, pt + normal);
    Vec3 intp;
    intersect_line(ln, intp);
    return intp;
  }

  bool is_same(const Plane3d& another_plane,
               SlickScalar angle_threshold = 5, SlickScalar dist_threshold = 0.01) const {
    if (is_parallel(another_plane, angle_threshold)) {
      Vec3 p = project_point(another_plane.point);
      if(((p-another_plane.point).norm() < dist_threshold) || dist_threshold < 0)
        return true;
    }
    return false;
  }

  bool is_parallel(const Plane3d& another_plane,
                   SlickScalar angle_threshold = 5) const{
    if (std::cos(std::fabs(normal.dot(another_plane.normal))) >
        std::cos(M_PI*angle_threshold/180))  {
        return true;
    }
    return false;
  }

  static const SlickScalar epsilon;  // plane precision
  Vec3 point;
  Vec3 normal;
  SlickScalar d;
};


inline std::pair<std::vector<unsigned int>, SlickScalar> Plane3d::FitPlane(
    const std::vector<Vec3 >& points, Plane3d& plane) {
  Eigen::MatrixXd v4_points(points.size(), 4);
  for (std::vector<Vec3 >::size_type i = 0; i < points.size(); i++) {
    v4_points.row(i) = Vec4(points[i][0], points[i][1], points[i][2], 1);
  }
  Eigen::Matrix4d A = v4_points.transpose()*v4_points;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> sym(A);
  SlickScalar smallest = 1.e+9;
  int sm_id = 0;
  for (int i = 0; i < 3; i++) {
    if (smallest > sym.eigenvalues()[i]) {
      smallest = sym.eigenvalues()[i];
      sm_id = i;
    }
  }
  Vec4 v4_plane = sym.eigenvectors().col(sm_id);

  // calculates error
  std::vector<SlickScalar> v_errsq;
  Tukey<SlickScalar> est;
  for (unsigned int i = 0; i < points.size(); i++) {
    const SlickScalar error = v4_points.row(i).dot(v4_plane);
    v_errsq.push_back(error*error);
  }
  est.compute_sigma_squared(v_errsq);

  std::vector<unsigned int> inlier_indices;
  SlickScalar sum_error = 0.0;
  for (unsigned int i = 0; i < points.size(); i++) {
    const SlickScalar error = v4_points.row(i).dot(v4_plane);
    const SlickScalar w = est.weight(error*error);
    sum_error += error*error*w;
    if (w > 0.7)
      inlier_indices.push_back(i);
  }
  plane = Plane3d(v4_plane);
  return std::make_pair(inlier_indices, sum_error/inlier_indices.size());
}

inline std::pair<std::vector<unsigned int>, SlickScalar> 
Plane3d::FitPlaneAdvance(const std::vector<Vec3 >& points, 
                            Plane3d& plane,
                            Vec3* anchor_point1,
                            Vec3* anchor_point2) {
  using namespace Eigen;
  int npoints = points.size();
  Vec3 v3_bestmean, v3_bestnormal;
  SlickScalar best_distsq = 1.e+9;
  int nransac = std::min(100, int(points.size()));
  Vec3 point2, point3;
  for (int i = 0; i < nransac; ++i) {
    int id1 = rand() % npoints;
    if (anchor_point1 == nullptr) {
      int id2 = id1;
      while (id2 == id1)
        id2 = rand()%npoints;
      point2 = points[id2];
      int id3 = id2;
      while (id3 == id2 || id3 == id1)
        id3 = rand() % npoints;
      point3 = points[id3];
    } else {
      point2 = *anchor_point1;
      if (anchor_point2 == nullptr) {
        int id3=id1;
        while (id3 == id1)
          id3 = rand()%npoints;
        point3 = points[id3];
      } else {
        point3 = *anchor_point2;
      }
    }

    Vec3 v3_mean = 0.33333333 *(points[id1]+point2+point3);
    Vec3 v3_ca = point3 - points[id1];
    Vec3 v3_ba = point2 - points[id1];
    Vec3 v3_normal = v3_ca.cross(v3_ba).normalized();

    if (v3_normal.norm() == 0)
      continue;
    v3_normal.normalized();

    std::vector<SlickScalar> v_errsq;
    Tukey<SlickScalar> est;
    for (int i = 0; i < npoints; i++) {
      Vector3d v3_diff = (points[i] - v3_mean).normalized();
      SlickScalar distsq = v3_diff.norm();
      if (distsq == 0.0)
        continue;
      v_errsq.push_back(fabs(v3_diff.dot(v3_normal)));
    }
    est.compute_sigma_squared(v_errsq);

    SlickScalar sum_error = 0.0;
    for (int i = 0; i < npoints; i++) {
      Vector3d v3_diff = (points[i] - v3_mean).normalized();
      SlickScalar distsq = v3_diff.norm();
      if (distsq == 0.0)
        continue;
      SlickScalar normdist = fabs(v3_diff.dot(v3_normal));

      if (normdist > 0.05)
        normdist = 0.05;
      sum_error += normdist*est.weight(normdist);
    }
    if (sum_error < best_distsq) {
      best_distsq = sum_error;
      v3_bestmean = v3_mean;
      v3_bestnormal = v3_normal;
    }
  }
  std::vector<unsigned> inlier_indices;
  for (int i = 0; i < npoints; i++) {
    Vector3d v3_diff = points[i] - v3_bestmean;
    SlickScalar distsq = v3_diff.norm();
    if (distsq == 0.0)
      continue;
    SlickScalar normdist = fabs(v3_diff.dot(v3_bestnormal));
    if (normdist < 0.05)
      inlier_indices.push_back(i);
  }

  // With these inliers, calculate mean and cov
  Vec3 plane_point = Vector3d::Zero();
  for (unsigned int i = 0; i < inlier_indices.size(); i++)
    plane_point += points[inlier_indices[i]];
  plane_point *= (1.0 / inlier_indices.size());

  Matrix3d m3_cov = Matrix3d::Zero();
  for (unsigned int i = 0; i < inlier_indices.size(); ++i) {
    Vector3d v3_diff = points[inlier_indices[i]] - plane_point;
    m3_cov += v3_diff * v3_diff.transpose();
  };

  // find the principal component with the minimal variance: this is the plane normal
  SelfAdjointEigenSolver<Matrix3d> sym(m3_cov);
  SlickScalar smallest = 1.e+9;
  int sm_id = 0;
  for (int i = 0; i < 3; i++) {
    if (smallest > sym.eigenvalues()[i]) {
      smallest = sym.eigenvalues()[i];
      sm_id = i;
    }
  }
  Vec3 plane_normal = sym.eigenvectors().col(sm_id);
  plane = Plane3d(plane_point, plane_normal);

  // calculates error
  std::vector<SlickScalar> v_errsq;
  Tukey<SlickScalar> est;
  for (unsigned int i = 0; i < inlier_indices.size(); i++) {
    int id2 = i;
    while (id2 == i)
      id2 = rand() % npoints;
    point2 = points[id2];
    int id3 = id2;
    while (id3 == id2 || id3 == i)
      id3 = rand() % npoints;
    point3 = points[id3];
    Vec3 v3_mean = 0.33333333 *(points[inlier_indices[i]]+point2+point3);
    Vec3 v3_diff = (points[inlier_indices[i]] - v3_mean).normalized();
    SlickScalar distsq = v3_diff.norm();
    if (distsq == 0.0)
      continue;
    v_errsq.push_back(fabs(v3_diff.dot(plane_normal)));
  }
  est.compute_sigma_squared(v_errsq);

  SlickScalar sum_error = 0.0;
  for (unsigned int i = 0; i < inlier_indices.size(); i++) {
    int id2 = i;
    while (id2 == i)
      id2 = rand() % npoints;
    point2 = points[id2];
    int id3 = id2;
    while (id3 == id2 || id3 == i)
      id3 = rand() % npoints;
    point3 = points[id3];
    Vec3 v3_mean = 0.33333333*(points[inlier_indices[i]]+point2+point3);
    Vec3 v3_diff = (points[inlier_indices[i]] - v3_mean).normalized();
    SlickScalar distsq = v3_diff.norm();
    if (distsq == 0.0)
      continue;
    SlickScalar normdist = fabs(v3_diff.dot(plane_normal));

    if (normdist > 0.05)
      normdist = 0.05;
    sum_error += normdist*est.weight(normdist);
  }

  return std::make_pair(inlier_indices, sum_error/inlier_indices.size());
}
}  // namespace slick