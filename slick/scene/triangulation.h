// Copyright 2012, 2013, 2014 The Look3D Authors. All rights reserved.
#ifndef LOOK3D_GEOMETRY_TRIANGULATE_H_
#define LOOK3D_GEOMETRY_TRIANGULATE_H_
#include <Eigen/Dense>
#include "math/se3.h"
#include "math/utilities.h"
#include "geometry/look3d_geometry_api.h"

namespace look3d {

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

template LOOK3D_GEOMETRY_API
std::pair<Eigen::Matrix<float, 3, 1>,float> triangulate(
    const Eigen::Matrix<float, 2, 1>& point1, const Eigen::Matrix<float, 2, 1>& point2,
    SE3Group<float> pose1,  SE3Group<float> pose2);
template LOOK3D_GEOMETRY_API

std::pair<Eigen::Matrix<double, 3, 1>,double> triangulate(
    const Eigen::Matrix<double, 2, 1>& point1, const Eigen::Matrix<double, 2, 1>& point2,
    SE3Group<double> pose1,  SE3Group<double> pose2);
}  // namespace look3d
#endif  // LOOK3D_GEOMETRY_TRIANGULATE_H_
