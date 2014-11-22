// Copyright 2012, 2013, 2014 The Look3D Authors. All rights reserved.
#ifndef LOOK3D_GEOMETRY_THREE_POINT_POSE_H_
#define LOOK3D_GEOMETRY_THREE_POINT_POSE_H_
#include <utility>
#include <vector>
#include <Eigen/Eigen>

#include "math/se3.h"
#include "math/mestimator.h"
#include "math/utilities.h"
#include "geometry/look3d_geometry_api.h"

namespace look3d {

// A function to evaluate x^4 + Bx^3 + Cx^2 + Dx + E
inline DefaultScalarType eval_quartic(DefaultScalarType B, DefaultScalarType C, DefaultScalarType D, DefaultScalarType E, DefaultScalarType x) {
  return E + x*(D + x*(C + x*(B + x)));
}

// A function that performs one iteration of Newton's method
// on the quartic x^4 + Bx^3 + Cx^2 + Dx + E
template<typename Scalar>
inline Scalar eval_newton_quartic(Scalar B, Scalar C, Scalar D, Scalar E, Scalar x) {
  Scalar fx = E + x*(D + x*(C + x*(B + x)));
  Scalar dx = D + x*(2*C + x*(3*B + x*4));
  return x - fx/dx;
}

template<typename Scalar>
static int DepressedCubicRealRoots(Scalar P, Scalar Q, Scalar r[]);

template<typename Scalar>
int FindQuarticRealRoots( Scalar B, Scalar C, Scalar D, Scalar E, Scalar r[]);

template<typename Scalar>
inline static Scalar square(Scalar x) { return x*x; }

template<typename Scalar>
static look3d::SE3Group<Scalar> GetAbsoluteOrientationFrom3Points(
    const Eigen::Matrix<Scalar, 3, 1> x[],
    const Eigen::Matrix<Scalar, 3, 1> y[]);

template LOOK3D_GEOMETRY_API SE3Group<float> GetAbsoluteOrientationFrom3Points(
const Eigen::Matrix<float, 3, 1> x[],
const Eigen::Matrix<float, 3, 1> y[]);

template LOOK3D_GEOMETRY_API SE3Group<double> GetAbsoluteOrientationFrom3Points(
const Eigen::Matrix<double, 3, 1> x[],
const Eigen::Matrix<double, 3, 1> y[]);

// The function for pose estimation from three 2D - 3D point correspondences.
// It implements the algorithm given by the Fischer and Bolles RANSAC paper,1980
// It assumes that the three points are in general position (not collinear).
// Input is an array of 3D cartesian positions and an array of
// 2D vectors that are the perspective projections of the points.
// Ouput is up to four poses satisfying the input with positive depths
// (points in front of the camera).
// @param[in] x an array containing at least 3 points
// @param[in] z an array containing the perspective projections
// of the points given by x in the current pose
// @param[out] poses the vector onto which any valid poses are appended
// @return the number of  poses appended to the vector
// @ingroup algorithms
template<typename Scalar>
int CalcThreePointPoses(const Eigen::Matrix<Scalar, 3, 1> xi[],
                        const Eigen::Matrix<Scalar, 2, 1> zi[],
                        std::vector<look3d::SE3Group<Scalar> >& poses);
template LOOK3D_GEOMETRY_API
int CalcThreePointPoses(const Eigen::Matrix<float, 3, 1> xi[],
                        const Eigen::Matrix<float, 2, 1> zi[],
                        std::vector<look3d::SE3Group<float> >& poses);
template LOOK3D_GEOMETRY_API
int CalcThreePointPoses(const Eigen::Matrix<double, 3, 1> xi[],
                        const Eigen::Matrix<double, 2, 1> zi[],
                        std::vector<look3d::SE3Group<double> >& poses);


template <typename Scalar>
std::pair<bool, look3d::SE3Group<Scalar> > ComputeRobustAbsolutPoseRANSACMLE(
    std::vector<std::pair<Eigen::Matrix<Scalar, 4, 1>,
                        Eigen::Matrix<Scalar, 2, 1> > > &observations);

template LOOK3D_GEOMETRY_API
std::pair<bool, look3d::SE3Group<float> > ComputeRobustAbsolutPoseRANSACMLE(
    std::vector<std::pair<Eigen::Matrix<float, 4, 1>,
                        Eigen::Matrix<float, 2, 1> > > &observations);

template LOOK3D_GEOMETRY_API
std::pair<bool, look3d::SE3Group<double> > ComputeRobustAbsolutPoseRANSACMLE(
    std::vector<std::pair<Eigen::Matrix<double, 4, 1>,
                        Eigen::Matrix<double, 2, 1> > > &observations);
}  // namespace look3d
#endif  // LOOK3D_GEOMETRY_THREE_POINT_POSE_H_
