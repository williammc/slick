// Copyright 2012, 2013, 2014 The Look3D Authors. All rights reserved.
#include "geometry/triangulation.h"

namespace look3d {

template<typename Scalar>
std::pair<Eigen::Matrix<Scalar, 3, 1>,Scalar> triangulate(
    const Eigen::Matrix<Scalar, 2, 1>& point1, const Eigen::Matrix<Scalar, 2, 1>& point2,
    SE3Group<Scalar> pose1,  SE3Group<Scalar> pose2) {
  Eigen::Matrix<Scalar, 3, 4> m34_pose1;
  m34_pose1.template block<3, 3>(0, 0) = pose1.get_rotation().get_matrix();
  m34_pose1.template block<3, 1>(0, 3) = pose1.get_translation();
  Eigen::Matrix<Scalar, 3, 4> m34_pose2;
  m34_pose2.template block<3, 3>(0, 0) = pose2.get_rotation().get_matrix();
  m34_pose2.template block<3, 1>(0, 3) = pose2.get_translation();
  Eigen::Matrix<Scalar, 4, 4> m44A;
  m44A.row(0) = point1[0]*m34_pose1.row(2) - m34_pose1.row(0);
  m44A.row(1) = point1[1]*m34_pose1.row(2) - m34_pose1.row(1);
  m44A.row(2) = point2[0]*m34_pose2.row(2) - m34_pose2.row(0);
  m44A.row(3) = point2[1]*m34_pose2.row(2) - m34_pose2.row(1);
  Eigen::JacobiSVD<Eigen::Matrix<Scalar, 4, 4> > svdA(m44A, Eigen::ComputeFullU
                                                      | Eigen::ComputeFullV);
  const Eigen::Matrix<Scalar, 4, 1> t = svdA.matrixV().col(3);
  Eigen::Matrix<Scalar, 3, 1> v3 = look3d::project(t);
  return std::make_pair(look3d::project(t),
                        ((point1-look3d::project(pose1*v3)).norm() +
                         (point2-look3d::project(pose2*v3)).norm())/2);
}

// instantiate =================================================================
template std::pair<Eigen::Matrix<DefaultScalarType, 3, 1>,DefaultScalarType> triangulate(
    const Eigen::Matrix<DefaultScalarType, 2, 1>& point1, const Eigen::Matrix<DefaultScalarType, 2, 1>& point2,
    SE3Group<DefaultScalarType> pose1,  SE3Group<DefaultScalarType> pose2);
template std::pair<Eigen::Matrix<float, 3, 1>,float> triangulate(
    const Eigen::Matrix<float, 2, 1>& point1, const Eigen::Matrix<float, 2, 1>& point2,
    SE3Group<float> pose1,  SE3Group<float> pose2);
}  // namespace look3d
