// Copyright 2014 The Slick Authors. All rights reserved.
#pragma once
#include <Eigen/Core>
#include "slick/math/se3.h"
#include "slick/util/common.h"

namespace slick {

// DLT method in MVG book(p.312)
// @param	point1[in]	2D point on camera plane(z=1) from 1st KeyFrame
// @param	point2[in]	2D point on camera plane(z=1) from 2nd KeyFrame
// @param	pose1[in] camera pose of the 1st frame
// @param	pose2[in] pose of the 2nd frame
// @return	(world point, reprojection error)
template<typename Scalar>
inline std::pair<Eigen::Matrix<Scalar, 3, 1>,Scalar> triangulate(
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
  Eigen::Matrix<Scalar, 3, 1> v3 = slick::project(t);
  return std::make_pair(slick::project(t),
                        ((point1-slick::project(pose1*v3)).norm() +
                         (point2-slick::project(pose2*v3)).norm())/2);
}
}  // namespace slick
