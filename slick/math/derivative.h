// Copyright 2014 The Slick Authors. All rights reserved.
#pragma once
#include "slick/math/so3.h"

namespace slick {

template<typename Deri31> Eigen::Matrix<typename Deri31::Scalar, 3, 3> dSO3(
    Eigen::MatrixBase<Deri31> const& v3) {
  Eigen::Matrix<typename Deri31::Scalar, 3, 3> derivative;
  derivative.col(0) = SO3Group<typename Deri31::Scalar>::generator_field(0, v3);
  derivative.col(1) = SO3Group<typename Deri31::Scalar>::generator_field(1, v3);
  derivative.col(2) = SO3Group<typename Deri31::Scalar>::generator_field(2, v3);
  return derivative;
}
}  // namespace slick
