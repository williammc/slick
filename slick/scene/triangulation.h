// Copyright 2014 The Slick Authors. All rights reserved.
#pragma once
#include <Eigen/Dense>
#include "slick/math/se3.h"
#include "slick/util/common.h"
#include "slick/slick_api.h"

namespace slick {

// DLT method in MVG book(p.312)
// @param	point1[in]	2D point on camera plane(z=1) from 1st KeyFrame
// @param	point2[in]	2D point on camera plane(z=1) from 2nd KeyFrame
// @param	pose1[in] camera pose of the 1st frame
// @param	pose2[in] pose of the 2nd frame
// @return	(world point, reprojection error)
template<typename Scalar>
std::pair<Eigen::Matrix<Scalar, 3, 1>,Scalar> triangulate(
    const Eigen::Matrix<Scalar, 2, 1>& point1, const Eigen::Matrix<Scalar, 2, 1>& point2,
    SE3Group<Scalar> pose1,  SE3Group<Scalar> pose2);

template SLICK_API
std::pair<Eigen::Matrix<float, 3, 1>,float> triangulate(
    const Eigen::Matrix<float, 2, 1>& point1, const Eigen::Matrix<float, 2, 1>& point2,
    SE3Group<float> pose1,  SE3Group<float> pose2);
template SLICK_API

std::pair<Eigen::Matrix<double, 3, 1>,double> triangulate(
    const Eigen::Matrix<double, 2, 1>& point1, const Eigen::Matrix<double, 2, 1>& point2,
    SE3Group<double> pose1,  SE3Group<double> pose2);
}  // namespace slick
