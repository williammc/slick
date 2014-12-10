// Copyright 2014 The Slick Authors. All rights reserved.
#include "slick/scene/atan_camera.h"

namespace slick {

template<typename Scalar>
AtanCamera<Scalar>::AtanCamera(AtanCamera const &cam_in) {
  image_size_ = cam_in.image_size_;
  params_ = cam_in.params_;
}

template<typename Scalar>
AtanCamera<Scalar>::AtanCamera(const Eigen::VectorXf& params) {
  Init(params);
}

template<typename Scalar>
void AtanCamera<Scalar>::Init(const Eigen::VectorXf& params) {
  image_size_ = params.head<2>().cast<Scalar>();
  params_ = params.segment(2, 5).cast<Scalar>();
}

template<typename Scalar>
Eigen::Matrix<Scalar, 2, 2> AtanCamera<Scalar>::GetProjectionDerivatives(
    const Eigen::Matrix<Scalar, 2, 1>& v2_cam) const {
  const Scalar r = v2_cam.norm();
  const Scalar d = 2.0 * tan(params_[4] / 2.0);
  const Scalar rf = radial_factor(r);

  const Scalar& x = v2_cam[0];
  const Scalar& y = v2_cam[1];

  Scalar g_frac_dx = 0.0;
  Scalar g_frac_dy = 0.0;
  if (r >= 0.01) {
    const Scalar r2 = r * r;
    const Scalar t = d / (r2*(1 + d*d*r2));
    g_frac_dx = x / params_[4] * t - x * rf / r2;
    g_frac_dy = y / params_[4] * t - y * rf / r2;
  }
  Eigen::Matrix<Scalar, 2, 2> derivs;
  derivs(0, 0) = params_[0] * (g_frac_dx * x + rf);
  derivs(0, 1) = params_[0] * g_frac_dy * x;
  derivs(1, 0) = params_[1] * g_frac_dx * y;
  derivs(1, 1) = params_[1] * (g_frac_dy * y + rf);
  return derivs;
}

template<typename Scalar>
void AtanCamera<Scalar>::SetParameters(const Eigen::Matrix<Scalar, 5, 1>& params) {
  params_ = params;
}


// instantiate =================================================================
template AtanCamera<double>::AtanCamera(AtanCamera<double> const &cam);
template AtanCamera<float>::AtanCamera(AtanCamera<float> const &cam);

template AtanCamera<double>::AtanCamera(const Eigen::VectorXf&);
template AtanCamera<float>::AtanCamera(const Eigen::VectorXf&);

template Eigen::Matrix<double, 2, 2> AtanCamera<double>::GetProjectionDerivatives(
    const Eigen::Matrix<double, 2, 1>& v2_cam) const;
template Eigen::Matrix<float, 2, 2> AtanCamera<float>::GetProjectionDerivatives(
    const Eigen::Matrix<float, 2, 1>& v2_cam) const;

template Eigen::Matrix<double, 2, 5> AtanCamera<double>::GetParameterDerivs(
    const Eigen::Matrix<double, 2, 1>& v2_cam);
template Eigen::Matrix<float, 2, 5> AtanCamera<float>::GetParameterDerivs(
    const Eigen::Matrix<float, 2, 1>& v2_cam);

template void AtanCamera<double>::SetParameters(
    const Eigen::Matrix<double, 5, 1>& vCPs);
template void AtanCamera<float>::SetParameters(
    const Eigen::Matrix<float, 5, 1>& vCPs);

}  // end namespace slick