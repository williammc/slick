// Copyright 2012, 2013, 2014 The Look3D Authors. All rights reserved.
#ifndef LOOK3D_MATH_DERIVATIVE_H_
#define LOOK3D_MATH_DERIVATIVE_H_
#include <math/so3.h>

namespace look3d {

template<typename Deri31>
Eigen::Matrix<typename Deri31::Scalar, 3, 3> dSO3(
    Eigen::MatrixBase<Deri31> const& v3) {
  Eigen::Matrix<typename Deri31::Scalar, 3, 3> derivative;
  derivative.col(0) = SO3Group<typename Deri31::Scalar>::generator_field(0, v3);
  derivative.col(1) = SO3Group<typename Deri31::Scalar>::generator_field(1, v3);
  derivative.col(2) = SO3Group<typename Deri31::Scalar>::generator_field(2, v3);
  return derivative;
}
}  // namespace look3d
#endif  // LOOK3D_MATH_DERIVATIVE_H_
