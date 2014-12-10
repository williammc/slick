// Copyright 2014 The Slick Authors. All rights reserved.
#include "slick/scene/poli_camera.h"

namespace slick {

template<typename Precision>
PoliCamera<Precision>::PoliCamera(PoliCamera const &cam_in) {
  Precision params[param_n_] = {cam_in.cam_params_[0],
                                 cam_in.cam_params_[1],
                                 cam_in.cam_params_[2],
                                 cam_in.cam_params_[3],
                                 cam_in.cam_params_[4],
                                 cam_in.cam_params_[5]};
  init(cam_in.width_, cam_in.height_, &params[0]);
}

template<typename Precision>
PoliCamera<Precision>::PoliCamera(
    int width, int height, const Precision* params) {
  init(width, height, params);
}

template<typename Precision>
PoliCamera<Precision>::PoliCamera(const Eigen::VectorXf& params) {
  width_ = params[0];
  height_ = params[1];
  cam_params_ = params.segment(2, 6).cast<Precision>();
  // update camera params to internal params
  UpdateInternalParams();
}

template<typename Precision>
void PoliCamera<Precision>::init(int width, int height,
                                        const Precision* params) {

  width_ = width;
  height_ = height;
  if (params != NULL) {
    for (int i = 0; i < param_n_; i++)
      cam_params_[i] = params[i];
  } else {  // default params
    cam_params_ << Precision(width_), Precision(width_), Precision(width_/2), Precision(height_/2), 0.0, 0.0;
  }

  // update camera params to internal params
  UpdateInternalParams();
}

template<typename Precision>
Eigen::Matrix<Precision, 2, 2> PoliCamera<Precision>::GetProjectionDerivatives(
    const Eigen::Matrix<Precision, 2, 1>& v2_camplane) const {
  Eigen::Matrix<Precision, 2, 2> m22Derivs =
      Eigen::Matrix<Precision, 2, 2>::Identity();
  Precision dTemp1 = v2_camplane.squaredNorm();
  Precision dTemp2 = k2_*dTemp1;
  m22Derivs *= 1+dTemp1*(k1_+dTemp2);
  m22Derivs += (2*(k1_+2*dTemp2)*v2_camplane)*v2_camplane.transpose();
  m22Derivs.row(0) *= fx_;
  m22Derivs.row(1) *= fy_;
  return m22Derivs;
}

template<typename Precision>
void PoliCamera<Precision>::UpdateInternalParams() {
  fx_ = cam_params_[0];
  fy_ = cam_params_[1];
  cx_ = cam_params_[2];
  cy_ = cam_params_[3];
  k1_ = cam_params_[4];
  k2_ = cam_params_[5];
  // find view angles
  Eigen::Matrix<Precision, 2, 1> v2, v2Left;
  //// find 1 pixel length
  //v2[0] = Precision(width_/2);
  //v2[1] = Precision(0);
  //v2 = UnProject(v2);
  //v2Left[0] = Precision(width_/2 - 1);  // 1 unit pixel different
  //v2Left[1] = Precision(0);
  //v2Left = UnProject(v2Left);
  //unit_pixel_distance_ = (v2Left - v2).norm();
}

template<typename Precision>
void PoliCamera<Precision>::set_parameters(
    const Eigen::Matrix<Precision, PoliCamera<Precision>::param_n_, 1>& vCPs) {
  cam_params_ = vCPs;
  UpdateInternalParams();
}


// instantiate =================================================================
template PoliCamera<double>::PoliCamera(PoliCamera<double> const &cam);
template PoliCamera<float>::PoliCamera(PoliCamera<float> const &cam);

template PoliCamera<double>::PoliCamera(const Eigen::VectorXf&);
template PoliCamera<float>::PoliCamera(const Eigen::VectorXf&);

template PoliCamera<double>::PoliCamera(int width, int height, const double* params);
template PoliCamera<float>::PoliCamera(int width, int height, const float* params);

template void PoliCamera<double>::init(int width, int height, const double* params);
template void PoliCamera<float>::init(int width, int height, const float* params);

template Eigen::Matrix<double, 2, 2> PoliCamera<double>::GetProjectionDerivatives(
  const Eigen::Matrix<double, 2, 1>& v2_camplane) const;
template Eigen::Matrix<float, 2, 2> PoliCamera<float>::GetProjectionDerivatives(
    const Eigen::Matrix<float, 2, 1>& v2_camplane) const;

template Eigen::Matrix<double, 2, PoliCamera<double>::param_n_> PoliCamera<double>::GetParameterDerivs(
  const Eigen::Matrix<double, 2, 1>& v2_camplane) const;
template Eigen::Matrix<float, 2, PoliCamera<float>::param_n_> PoliCamera<float>::GetParameterDerivs(
    const Eigen::Matrix<float, 2, 1>& v2_camplane) const;

template void PoliCamera<double>::UpdateInternalParams();
template void PoliCamera<float>::UpdateInternalParams();

template void PoliCamera<double>::set_parameters(
  const Eigen::Matrix<double, PoliCamera<double>::param_n_, 1>& vCPs);
template void PoliCamera<float>::set_parameters(
    const Eigen::Matrix<float, PoliCamera<float>::param_n_, 1>& vCPs);

}  // end namespace slick

