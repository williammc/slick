// Copyright 2014 The Slick Authors. All rights reserved.
// I'm camera man, in mathematical aspect.
// You can know me more through references & description below.
// Welcome suggestions & comments about camera's models & calibration techniques.
#pragma once
#include <memory>
#include <Eigen/Eigen>
#include "slick/math/se3.h"
#include "slick/util/common.h"
#include "slick/slick_api.h"

namespace slick {

/// Polynomial Camera model is Pin Hole projection model
/// applying Taylor series for radial distortion
/// (u,v) = ProjectionMatrix*(R_/R)*(x,y)
/// Pixel coordinate (u,v) start from Left-Top corner,
/// with vector u is along hirozontal left-right, vector v is vertical top-down
/// Responsibilities: Projecttion(3D to 2D) &
/// Unprojection(2D to 3D, at camera plane where Zc=1.)
template<typename Precision = DefaultScalarType>
class SLICK_API PoliCamera {
 public:
  typedef std::shared_ptr<PoliCamera<Precision> > Ptr;
  typedef std::shared_ptr<const PoliCamera<Precision> > ConstPtr;
  typedef Precision ScalarType;
  /// Camera params number (fx,fy,cx,cy,k1,k2)
  static const int param_n_ = 6;
  // Constructors  =============================================================
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  /// explicit default constructor
  PoliCamera() {
    init(640, 480);
  }

  /// convenient constructor
  PoliCamera(PoliCamera const &cam_in);

  /// constructor with width, height, & camera params
  PoliCamera(int width, int height, const Precision* params = NULL);

  /// Key Methods  ==============================================================
  /// Do projection.
  /// @param v3Cam[in] 3D position in camera coordinate frame
  /// @return 2D point on image plane (on screen)
  template<typename OtherDerived>
  Eigen::Matrix<Precision, 2, 1> Project(
      const Eigen::MatrixBase<OtherDerived>& v2_camplane) const;

  /// Do unprojection
  /// @param v2_implane[in] 2D point on image plane (on screen)
  /// @return 2D point on camera plane (Zc=1)
  template<typename OtherDerived>
  Eigen::Matrix<Precision, 2, 1> UnProject(
      const Eigen::MatrixBase<OtherDerived>& v2_implane) const;

  /// Project linearly (only apply PinHole part, no distortion)
  /// 3D point into image plane
  /// @param  v2_camplane[in] 3D point
  /// @return projected 2D point without distortion
  template<typename OtherDerived>
  Eigen::Matrix<Precision, 2, 1> ProjectLinear(
      const Eigen::MatrixBase<OtherDerived>& v2_camplane) const;

  /// Unproject an undistorting 2D into 3D point in image plane with z=1
  /// @param vIn2DPoint undistorting 2D point
  /// @return 3D point at image plane that z=1
  template<typename OtherDerived>
  Eigen::Matrix<Precision, 2, 1> UnProjectLinear(
      const Eigen::MatrixBase<OtherDerived>& vIn2DPoint) const;

  Eigen::Matrix<Precision, 2, 2> GetProjectionDerivatives(
      const Eigen::Matrix<Precision, 2, 1>& v2_camplane) const;

  Precision unit_pixel_distance() {
    return unit_pixel_distance_;
  }

  /// return (width,height) resolution of the camera
  void get_resolution(int& width, int& height) const {
    width = width_;
    height = height_;
  }

  /// reset image resolution for this camera model
  void set_resolution(int w, int h) {
    int old_w = width_, old_h = height_;
    width_ = w;
    height_ = h;
    fx_ *= static_cast<DefaultScalarType>(w)/old_w;
    fy_ *= static_cast<DefaultScalarType>(h)/old_h;
    cx_ *= static_cast<DefaultScalarType>(w)/old_w;
    cy_ *= static_cast<DefaultScalarType>(h)/old_h;
    UpdateInternalParams();
  }

  int width() const {
    return width_;
  }

  int height() const {
    return height_;
  }

  // For Calibration purposes  =================================================
  Eigen::Matrix<Precision, 2, PoliCamera<Precision>::param_n_> GetParameterDerivs(
      const Eigen::Matrix<Precision, 2, 1>& v2_camplane) const {
    Eigen::Matrix<Precision, 2, PoliCamera<Precision>::param_n_> result;
    Precision r2 = v2_camplane.squaredNorm();
    Precision r4 = r2 * r2;
    Eigen::Matrix<Precision, 2, 1> mod_camframe =
        v2_camplane * (1+ r2 * (cam_params_[4] + r2 * cam_params_[5]));

    result(0, 0) = mod_camframe[0];
    result(0, 1) = 0;
    result(0, 2) = 1;
    result(0, 3) = 0;
    result(0, 4) = cam_params_[0]*v2_camplane[0]*r2;
    result(0, 5) = cam_params_[0]*v2_camplane[0]*r4;

    result(1, 0) = 0;
    result(1, 1) = mod_camframe[1];
    result(1, 2) = 0;
    result(1, 3) = 1;
    result(1, 4) = cam_params_[1]*v2_camplane[1]*r2;
    result(1, 5) = cam_params_[1]*v2_camplane[1]*r4;
    return result;
  }

  const Eigen::Matrix<Precision, PoliCamera<Precision>::param_n_, 1>& parameters() const {
    return cam_params_;
  }
  void set_parameters(
      const Eigen::Matrix<Precision, PoliCamera<Precision>::param_n_, 1>& vCPs);

 protected:
  /// common method for initializing internal params
  void init(int width, int height, const Precision* params = NULL);
  /// update internal params
  void UpdateInternalParams();

  /// The intrinsic parameters in the camera
  Eigen::Matrix<Precision, PoliCamera::param_n_, 1> cam_params_;
  int width_;  /// the width of the camera image
  int height_;  /// the height of the camera image
  Precision unit_pixel_distance_;  /// dist of 1-unit pixel on camera plane (z=1)
  /// Camera's instrinsic parameters
  Precision fx_, fy_;  /// focal length
  Precision cx_, cy_;  /// principal point
  Precision k1_, k2_;  /// Distortion Coefficients
};

// Implementation ============================ =================================
/// const int iParamN = PoliCamera<>::iParamN;
/// explicit declared to advoide compiler's specific configuration

template<typename Precision>
template<typename OtherDerived>
inline Eigen::Matrix<Precision, 2, 1> PoliCamera<Precision>::Project(
    const Eigen::MatrixBase<OtherDerived>& v2_camplane) const {
  Eigen::Matrix<Precision, 2, 1> v2Point;
  Precision dR2, dR4, dFactor;
  v2Point[0] = v2_camplane[0];
  v2Point[1] = v2_camplane[1];
  dR2 = v2Point[0]*v2Point[0] + v2Point[1]*v2Point[1];
  dR4 = dR2*dR2;
  dFactor = 1 + k1_*dR2 + k2_*dR4;
  v2Point[0] = fx_*v2Point[0]*dFactor + cx_;
  v2Point[1] = fy_*v2Point[1]*dFactor + cy_;
  return v2Point;
}

template<typename Precision>
template<typename OtherDerived>
inline Eigen::Matrix<Precision, 2, 1> PoliCamera<Precision>::UnProject(
    const Eigen::MatrixBase<OtherDerived>& v2_implane) const {
  Eigen::Matrix<Precision, 2, 1> v2Point;
  v2Point = v2_implane;
  v2Point[0] = v2_implane(0, 0);
  v2Point[1] = v2_implane(1, 0);
  v2Point[0] = (v2Point[0]-cx_)/fx_;
  v2Point[1] = (v2Point[1]-cy_)/fy_;
  /// employ libcvd's implementation
  /// first guess
  DefaultScalarType scale = v2Point.squaredNorm();

  /// iterations of Newton-Rapheson
  for (int i = 0; i < 3 ; ++i) {  /// 3 iteration is good enough in certain range
    DefaultScalarType temp = 1 + scale*(k1_+k2_*scale);
    DefaultScalarType error = v2Point.squaredNorm() - scale*temp*temp;
    DefaultScalarType deriv = temp*(temp+2*scale*(k1_ + 2*k2_*scale));
    scale += error/deriv;
  }
  v2Point = v2Point/static_cast<Precision>((1 + scale*(k1_ + k2_*scale)));
  return v2Point;
}

template<typename Precision>
template<typename OtherDerived>
inline Eigen::Matrix<Precision, 2, 1> PoliCamera<Precision>::ProjectLinear(
    const Eigen::MatrixBase<OtherDerived>& v2_camplane) const {
  Eigen::Matrix<Precision, 2, 1> v2Point;
  v2Point[0] = fx_*(v2_camplane[0]) + cx_;
  v2Point[1] = fy_*(v2_camplane[1]) + cy_;
  return v2Point;
}

template<typename Precision>
template<typename OtherDerived>
inline Eigen::Matrix<Precision, 2, 1> PoliCamera<Precision>::UnProjectLinear(
    const Eigen::MatrixBase<OtherDerived>& v2_implane) const {
  Eigen::Matrix<Precision, 2, 1> v2Point;
  /// Transfer from Pixel coordinate to (u,v) coordinate
  v2Point[0] = v2_implane[0];
  v2Point[1] = v2_implane[1];
  v2Point[0] = (v2Point[0]-cx_)/fx_;
  v2Point[1] = (v2Point[1]-cy_)/fy_;
  return v2Point;
}
}  // end namespace slick
