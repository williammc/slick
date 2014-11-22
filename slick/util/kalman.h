// Copyright 2014 The Slick Authors. All rights reserved.
#pragma once

#include <Eigen/Cholesky>

namespace slick {
template<typename DeriSMxSM, typename DeriMMx1,
         typename DeriMMxSM, typename DeriMMxMM>
Eigen::Matrix<typename DeriSMxSM::Scalar, DeriSMxSM::RowsAtCompileTime, 1>
kalman_update(Eigen::MatrixBase<DeriSMxSM> & P,
              const Eigen::MatrixBase<DeriMMx1> &innovation,
              const Eigen::MatrixBase<DeriMMxSM> &H,
              const Eigen::MatrixBase<DeriMMxMM> &measurement_noise) {
  const Eigen::Matrix<typename DeriSMxSM::Scalar, DeriSMxSM::RowsAtCompileTime,
      DeriMMxMM::RowsAtCompileTime> P12 = P * H.transpose();
#if 0  // this only works for R*R' = chol(X), but TooN Cholesky is L * D * L':/
  TooN::Cholesky<MEASDIM, T> R(H * P12 + measurement_noise);
  const TooN::Matrix<STATEDIM, MEASDIM, T> U = R.backsub(P12.T()).T();
  P = P - U * U.T();
  return U * R.backsub(innovation);
#endif
  Eigen::LLT<Eigen::Matrix<typename DeriMMxMM::Scalar,
      DeriMMxMM::RowsAtCompileTime,
      DeriMMxMM::RowsAtCompileTime> > R(H * P12 + measurement_noise);
  P = P - P12 * R.solve(P12.transpose());
  return P12 * R.solve(innovation);
}
}  // namespace slick
