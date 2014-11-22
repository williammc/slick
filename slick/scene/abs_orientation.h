// Copyright 2014 The Slick Authors. All rights reserved.
#pragma once
#include <vector>
#include "slick/math/se3.h"
#include "slick/datatypes.h"
#include "slick/slick_api.h"

namespace slick {

// @defgroup absorient Absolute Orientation
// contains various functions to calculate rotations and rigid transformations
// between sets of 3D correspondences.

// computes the rotation between two sets of points maximizing b * Ta
// this function is part of the absolute orientation algorithm after Horn
// and is used in @ref computeAbsoluteOrientation function.
// @param[in] a vector of 3D points
// @param[in] b vector of 3D points
// @return SO3 containing the rotation such that b = T a
// @ingroup absorient
SLICK_API SO3 ComputeOrientation(const std::vector<Eigen::Matrix<DefaultScalarType, 3, 1> > & a,
                         const std::vector<Eigen::Matrix<DefaultScalarType, 3, 1> > & b);

// computes rotation between two pairs of rays in space maximizing b * T a
// its 8x faster then using the general ComputeOrientation for 2 correspondences
// @param[in] a1 first input vector
// @param[in] b1 first output vector
// @param[in] a2 second input vector
// @param[in] b2 second output vector
// @return SO3 containing the rotation such that b = T a
// @ingroup absorient
SLICK_API SO3 ComputeOrientation(const Eigen::Matrix<DefaultScalarType, 3, 1> & a1,
                         const Eigen::Matrix<DefaultScalarType, 3, 1> & b1,
                         const Eigen::Matrix<DefaultScalarType, 3, 1> & a2,
                         const Eigen::Matrix<DefaultScalarType, 3, 1> & b2);

// computes rigid transformation between two corresponding point sets after Horn
// result is an SE3 that maps points from vector a to points from
// vector b : b[i] = SE3 * a[i]
// @param[in] a vector of 3D points
// @param[in] b vector of 3D points
// @return SE3 containing the transformation such that b = T a
// @ingroup absorient
SLICK_API SE3 ComputeAbsoluteOrientation(const std::vector<Eigen::Matrix<DefaultScalarType, 3, 1> > & a,
                                 const std::vector<Eigen::Matrix<DefaultScalarType, 3, 1> > & b);

// computes the mean rotation of a set of rotations.
// This is the rotation R such that R^{-1} * R_i is minimal for all R_i.
// @param[in] r a vector of rotations
// @return SO3 mean rotation of all input rotations
// @ingroup absorient
SLICK_API SO3 ComputeMeanOrientation(const std::vector<SO3 > & r);

// computes a rotation matrix corresponding to a unit quaternion. The quaternion
// is in the format (q0,qx,qy,qz) to fit the absolute orientation algorithm.This
// is a helper function for the @ref computeOrientation function.
// @param[in] q unit quaternion as (q0,qx,qy,qz)
// @return a 3x3 rotation matrix corresponding to the quaternion
// @ingroup absorient
SLICK_API Eigen::Matrix<DefaultScalarType, 3, 3> QuaternionToMatrix(
  const Eigen::Matrix<DefaultScalarType, 4, 1> & q);

}  // namespace slick
