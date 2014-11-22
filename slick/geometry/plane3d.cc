// Copyright 2014 The Slick Authors. All rights reserved.
#include "slick/geometry/plane3d.h"
#include "slick/util/mestimator.h"

namespace slick {
const DefaultScalarType Plane3d::epsilon = 1.e-30;

std::pair<std::vector<unsigned int>, DefaultScalarType> Plane3d::fit_plane(
    const std::vector<Eigen::Matrix<DefaultScalarType, 3, 1> >& points, Plane3d& plane) {
  Eigen::MatrixXd v4_points(points.size(), 4);
  for (std::vector<Eigen::Matrix<DefaultScalarType, 3, 1> >::size_type i = 0; i < points.size(); i++) {
    v4_points.row(i) = Eigen::Matrix<DefaultScalarType, 4, 1>(points[i][0], points[i][1], points[i][2], 1);
  }
  Eigen::Matrix4d A = v4_points.transpose()*v4_points;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> sym(A);
  DefaultScalarType smallest = 1.e+9;
  int sm_id = 0;
  for (int i = 0; i < 3; i++) {
    if (smallest > sym.eigenvalues()[i]) {
      smallest = sym.eigenvalues()[i];
      sm_id = i;
    }
  }
  Eigen::Matrix<DefaultScalarType, 4, 1> v4_plane = sym.eigenvectors().col(sm_id);

  // calculates error
  std::vector<DefaultScalarType> v_errsq;
  Tukey<DefaultScalarType> est;
  for (unsigned int i = 0; i < points.size(); i++) {
    const DefaultScalarType error = v4_points.row(i).dot(v4_plane);
    v_errsq.push_back(error*error);
  }
  est.compute_sigma_squared(v_errsq);

  std::vector<unsigned int> inlier_indices;
  DefaultScalarType sum_error = 0.0;
  for (unsigned int i = 0; i < points.size(); i++) {
    const DefaultScalarType error = v4_points.row(i).dot(v4_plane);
    const DefaultScalarType w = est.weight(error*error);
    sum_error += error*error*w;
    if (w > 0.7)
      inlier_indices.push_back(i);
  }
  plane = Plane3d(v4_plane);
  return std::make_pair(inlier_indices, sum_error/inlier_indices.size());
}

std::pair<std::vector<unsigned int>, DefaultScalarType> Plane3d::fit_plane_advanced(
    const std::vector<Eigen::Matrix<DefaultScalarType, 3, 1> >& points, Plane3d& plane,
    Eigen::Matrix<DefaultScalarType, 3, 1>* anchor_point1, Eigen::Matrix<DefaultScalarType, 3, 1>* anchor_point2) {
  using namespace Eigen;
  int npoints = points.size();
  Eigen::Matrix<DefaultScalarType, 3, 1> v3_bestmean, v3_bestnormal;
  DefaultScalarType best_distsq=1.e+9;
  int nransac = std::min(100, int(points.size()));
  Eigen::Matrix<DefaultScalarType, 3, 1> point2, point3;
  for (int i = 0; i < nransac; ++i) {
    int id1 = rand() % npoints;
    if (anchor_point1 == NULL) {
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
      if (anchor_point2 == NULL) {
        int id3=id1;
        while (id3 == id1)
          id3 = rand()%npoints;
        point3 = points[id3];
      } else {
        point3 = *anchor_point2;
      }
    }

    Eigen::Matrix<DefaultScalarType, 3, 1> v3_mean = 0.33333333 *(points[id1]+point2+point3);
    Eigen::Matrix<DefaultScalarType, 3, 1> v3_ca = point3 - points[id1];
    Eigen::Matrix<DefaultScalarType, 3, 1> v3_ba = point2 - points[id1];
    Eigen::Matrix<DefaultScalarType, 3, 1> v3_normal = v3_ca.cross(v3_ba).normalized();

    if (v3_normal.norm() == 0)
      continue;
    v3_normal.normalized();

    std::vector<DefaultScalarType> v_errsq;
    Tukey<DefaultScalarType> est;
    for (int i = 0; i < npoints; i++) {
      Vector3d v3_diff = (points[i] - v3_mean).normalized();
      DefaultScalarType distsq = v3_diff.norm();
      if (distsq == 0.0)
        continue;
      v_errsq.push_back(fabs(v3_diff.dot(v3_normal)));
    }
    est.compute_sigma_squared(v_errsq);

    DefaultScalarType sum_error = 0.0;
    for (int i = 0; i < npoints; i++) {
      Vector3d v3_diff = (points[i] - v3_mean).normalized();
      DefaultScalarType distsq = v3_diff.norm();
      if (distsq == 0.0)
        continue;
      DefaultScalarType normdist = fabs(v3_diff.dot(v3_normal));

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
    DefaultScalarType distsq = v3_diff.norm();
    if (distsq == 0.0)
      continue;
    DefaultScalarType normdist = fabs(v3_diff.dot(v3_bestnormal));
    if (normdist < 0.05)
      inlier_indices.push_back(i);
  }

  // With these inliers, calculate mean and cov
  Eigen::Matrix<DefaultScalarType, 3, 1> plane_point = Vector3d::Zero();
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
  DefaultScalarType smallest = 1.e+9;
  int sm_id = 0;
  for (int i = 0; i < 3; i++) {
    if (smallest > sym.eigenvalues()[i]) {
      smallest = sym.eigenvalues()[i];
      sm_id = i;
    }
  }
  Eigen::Matrix<DefaultScalarType, 3, 1> plane_normal = sym.eigenvectors().col(sm_id);
  plane = Plane3d(plane_point, plane_normal);

  // calculates error
  std::vector<DefaultScalarType> v_errsq;
  Tukey<DefaultScalarType> est;
  for (unsigned int i = 0; i < inlier_indices.size(); i++) {
    int id2 = i;
    while (id2 == i)
      id2 = rand() % npoints;
    point2 = points[id2];
    int id3 = id2;
    while (id3 == id2 || id3 == i)
      id3 = rand() % npoints;
    point3 = points[id3];
    Eigen::Matrix<DefaultScalarType, 3, 1> v3_mean = 0.33333333 *(points[inlier_indices[i]]+point2+point3);
    Eigen::Matrix<DefaultScalarType, 3, 1> v3_diff = (points[inlier_indices[i]] - v3_mean).normalized();
    DefaultScalarType distsq = v3_diff.norm();
    if (distsq == 0.0)
      continue;
    v_errsq.push_back(fabs(v3_diff.dot(plane_normal)));
  }
  est.compute_sigma_squared(v_errsq);

  DefaultScalarType sum_error = 0.0;
  for (unsigned int i = 0; i < inlier_indices.size(); i++) {
    int id2 = i;
    while (id2 == i)
      id2 = rand() % npoints;
    point2 = points[id2];
    int id3 = id2;
    while (id3 == id2 || id3 == i)
      id3 = rand() % npoints;
    point3 = points[id3];
    Eigen::Matrix<DefaultScalarType, 3, 1> v3_mean = 0.33333333*(points[inlier_indices[i]]+point2+point3);
    Eigen::Matrix<DefaultScalarType, 3, 1> v3_diff = (points[inlier_indices[i]] - v3_mean).normalized();
    DefaultScalarType distsq = v3_diff.norm();
    if (distsq == 0.0)
      continue;
    DefaultScalarType normdist = fabs(v3_diff.dot(plane_normal));

    if (normdist > 0.05)
      normdist = 0.05;
    sum_error += normdist*est.weight(normdist);
  }

  return std::make_pair(inlier_indices, sum_error/inlier_indices.size());
}
}  // namespace slick
