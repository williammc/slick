#pragma once
#include "slick/geometry/plane3d.h"
#include "slick/util/mestimator.h"

namespace slick {

template <typename Scalar>
inline std::pair<std::vector<int>, Scalar>
LeastSquareFitPlane(const std::vector<Eigen::Matrix<Scalar, 3, 1>> &points,
                    Plane3DBase<Scalar> &plane) {
  Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> v4_points(points.size(),
                                                                  4);
  for (auto pt : points) {
    v4_points.row(i) = Eigen::Matrix<Scalar, 4, 1>(pt[0], pt[1], pt[2], 1);
  }
  Eigen::Matrix<Scalar, 4, 4> A = v4_points.transpose() * v4_points;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Scalar, 4, 4>> sym(A);
  Scalar smallest = 1.e+9;
  int sm_id = 0;
  for (int i = 0; i < 3; i++) {
    if (smallest > sym.eigenvalues()[i]) {
      smallest = sym.eigenvalues()[i];
      sm_id = i;
    }
  }
  Eigen::Matrix<Scalar, 4, 1> v4_plane = sym.eigenvectors().col(sm_id);

  // calculates error
  std::vector<Scalar> v_errsq;
  Tukey<Scalar> est;
  for (int i = 0; i < points.size(); i++) {
    const Scalar error = v4_points.row(i).dot(v4_plane);
    v_errsq.push_back(error * error);
  }
  est.compute_sigma_squared(v_errsq);

  std::vector<int> inlier_indices;
  Scalar sum_error = 0.0;
  for (int i = 0; i < points.size(); i++) {
    const Scalar error = v4_points.row(i).dot(v4_plane);
    const Scalar w = est.weight(error * error);
    sum_error += error * error * w;
    if (w > 0.7)
      inlier_indices.push_back(i);
  }
  plane = Plane3DBase<Scalar>(v4_plane);
  return std::make_pair(inlier_indices, sum_error / inlier_indices.size());
}

template <typename Scalar>
inline std::pair<std::vector<int>, Scalar>
RobustFitPlane(const std::vector<Eigen::Matrix<Scalar, 3, 1>> &points,
                         Plane3DBase<Scalar> &plane,
                         Eigen::Matrix<Scalar, 3, 1> *anchor_point1,
                         Eigen::Matrix<Scalar, 3, 1> *anchor_point2) {
  int npoints = points.size();
  Eigen::Matrix<Scalar, 3, 1> v3_bestmean, v3_bestnormal;
  Scalar best_distsq = 1.e+9;
  int nransac = std::min(100, int(points.size()));
  Eigen::Matrix<Scalar, 3, 1> point2, point3;
  for (int i = 0; i < nransac; ++i) {
    int id1 = rand() % npoints;
    if (anchor_point1 == nullptr) {
      int id2 = id1;
      while (id2 == id1)
        id2 = rand() % npoints;
      point2 = points[id2];
      int id3 = id2;
      while (id3 == id2 || id3 == id1)
        id3 = rand() % npoints;
      point3 = points[id3];
    } else {
      point2 = *anchor_point1;
      if (anchor_point2 == nullptr) {
        int id3 = id1;
        while (id3 == id1)
          id3 = rand() % npoints;
        point3 = points[id3];
      } else {
        point3 = *anchor_point2;
      }
    }

    Eigen::Matrix<Scalar, 3, 1> v3_mean =
        0.33333333 * (points[id1] + point2 + point3);
    Eigen::Matrix<Scalar, 3, 1> v3_ca = point3 - points[id1];
    Eigen::Matrix<Scalar, 3, 1> v3_ba = point2 - points[id1];
    Eigen::Matrix<Scalar, 3, 1> v3_normal = v3_ca.cross(v3_ba).normalized();

    if (v3_normal.norm() == 0)
      continue;
    v3_normal.normalized();

    std::vector<Scalar> v_errsq;
    Tukey<Scalar> est;
    for (int i = 0; i < npoints; i++) {
      Vector3d v3_diff = (points[i] - v3_mean).normalized();
      Scalar distsq = v3_diff.norm();
      if (distsq == 0.0)
        continue;
      v_errsq.push_back(fabs(v3_diff.dot(v3_normal)));
    }
    est.compute_sigma_squared(v_errsq);

    Scalar sum_error = 0.0;
    for (int i = 0; i < npoints; i++) {
      Vector3d v3_diff = (points[i] - v3_mean).normalized();
      Scalar distsq = v3_diff.norm();
      if (distsq == 0.0)
        continue;
      Scalar normdist = fabs(v3_diff.dot(v3_normal));

      if (normdist > 0.05)
        normdist = 0.05;
      sum_error += normdist * est.weight(normdist);
    }
    if (sum_error < best_distsq) {
      best_distsq = sum_error;
      v3_bestmean = v3_mean;
      v3_bestnormal = v3_normal;
    }
  }
  std::vector<int> inlier_indices;
  for (int i = 0; i < npoints; i++) {
    Vector3d v3_diff = points[i] - v3_bestmean;
    Scalar distsq = v3_diff.norm();
    if (distsq == 0.0)
      continue;
    Scalar normdist = fabs(v3_diff.dot(v3_bestnormal));
    if (normdist < 0.05)
      inlier_indices.push_back(i);
  }

  // With these inliers, calculate mean and cov
  Eigen::Matrix<Scalar, 3, 1> plane_point = Vector3d::Zero();
  for (int i = 0; i < inlier_indices.size(); i++)
    plane_point += points[inlier_indices[i]];
  plane_point *= (1.0 / inlier_indices.size());

  Eigen::Matrix<Scalar, 3, 3> m3_cov = Eigen::Matrix<Scalar, 3, 3>::Zero();
  for (int i = 0; i < inlier_indices.size(); ++i) {
    Vector3d v3_diff = points[inlier_indices[i]] - plane_point;
    m3_cov += v3_diff * v3_diff.transpose();
  };

  // find the principal component with the minimal variance: this is the plane
  // normal
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Scalar, 3, 3>> sym(m3_cov);
  Scalar smallest = 1.e+9;
  int sm_id = 0;
  for (int i = 0; i < 3; i++) {
    if (smallest > sym.eigenvalues()[i]) {
      smallest = sym.eigenvalues()[i];
      sm_id = i;
    }
  }
  Eigen::Matrix<Scalar, 3, 1> plane_normal = sym.eigenvectors().col(sm_id);
  plane = Plane3d(plane_point, plane_normal);

  // calculates error
  std::vector<Scalar> v_errsq;
  Tukey<Scalar> est;
  for (int i = 0; i < inlier_indices.size(); i++) {
    int id2 = i;
    while (id2 == i)
      id2 = rand() % npoints;
    point2 = points[id2];
    int id3 = id2;
    while (id3 == id2 || id3 == i)
      id3 = rand() % npoints;
    point3 = points[id3];
    Eigen::Matrix<Scalar, 3, 1> v3_mean =
        0.33333333 * (points[inlier_indices[i]] + point2 + point3);
    Eigen::Matrix<Scalar, 3, 1> v3_diff =
        (points[inlier_indices[i]] - v3_mean).normalized();
    Scalar distsq = v3_diff.norm();
    if (distsq == 0.0)
      continue;
    v_errsq.push_back(fabs(v3_diff.dot(plane_normal)));
  }
  est.compute_sigma_squared(v_errsq);

  Scalar sum_error = 0.0;
  for (int i = 0; i < inlier_indices.size(); i++) {
    int id2 = i;
    while (id2 == i)
      id2 = rand() % npoints;
    point2 = points[id2];
    int id3 = id2;
    while (id3 == id2 || id3 == i)
      id3 = rand() % npoints;
    point3 = points[id3];
    Eigen::Matrix<Scalar, 3, 1> v3_mean =
        0.33333333 * (points[inlier_indices[i]] + point2 + point3);
    Eigen::Matrix<Scalar, 3, 1> v3_diff =
        (points[inlier_indices[i]] - v3_mean).normalized();
    Scalar distsq = v3_diff.norm();
    if (distsq == 0.0)
      continue;
    Scalar normdist = std::fabs(v3_diff.dot(plane_normal));

    if (normdist > 0.05)
      normdist = 0.05;
    sum_error += normdist * est.weight(normdist);
  }

  return std::make_pair(inlier_indices, sum_error / inlier_indices.size());
}

} // namespace slick