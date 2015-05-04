// Copyright 2014 The Slick Authors. All rights reserved.
// I'm camera man, in mathematical aspect.
// You can know me more through references & description below.
// Welcome suggestions & comments about camera's models & calibration techniques.
#pragma once
#include <memory>
#include <Eigen/Core>
#include "slick/math/se3.h"
#include "slick/util/common.h"

namespace slick {

/// Polynomial Camera model is Pin Hole projection model
/// applying Taylor series for radial distortion
/// (u,v) = ProjectionMatrix*(R_/R)*(x,y)
/// Pixel coordinate (u,v) start from Left-Top corner,
/// with vector u is along hirozontal left-right, vector v is vertical top-down
/// Responsibilities: Projecttion(3D to 2D) &
/// Unprojection(2D to 3D, at camera plane where Zc=1.)
template<typename Scalar = SlickScalar>
class PoliCamera {
  typedef Eigen::Matrix<Scalar, 2, 1> Vec2_t;
  typedef Eigen::Matrix<Scalar, 3, 1> Vec3_t;
  typedef Eigen::Matrix<Scalar, 6, 1> Vec6_t;
  typedef Eigen::Matrix<Scalar, 8, 1> Vec8_t;
  typedef Eigen::Matrix<Scalar, 2, 2> Mat2_t;
 public:
   typedef std::shared_ptr<PoliCamera> Ptr;
  typedef Scalar ScalarType;
  /// Camera params number (fx,fy,cx,cy,k1,k2)
  static const int param_n_ = 6;
  // Constructors --------------------------------------------------------------
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  /// explicit default constructor
  PoliCamera() {
    Eigen::VectorXf v(8);
    v << 640, 480, 500, 500, 320, 240, 0, 0;
    Init(v);
  }

  /// convenient constructor
  PoliCamera(PoliCamera const &cam) {
    image_size_ = cam.image_size_;
    params_ = cam.params_;
  }

  PoliCamera(const Eigen::VectorXf& params) {
    Init(params);
  }

  /// Key Methods --------------------------------------------------------------
  /// Projection of 2D/3D point on camera plane (z=1) to image plane
  template<typename OtherDerived>
  Vec2_t Project(const Eigen::MatrixBase<OtherDerived>& pt) const{
    Vec2_t loc;
    Scalar r2, r4, factor;
    loc[0] = pt[0];
    loc[1] = pt[1];
    r2 = loc[0]*loc[0] + loc[1]*loc[1];
    r4 = r2*r2;
    factor = 1 + params_[4]*r2 + params_[5]*r4;
    loc[0] = params_[0]*loc[0]*factor + params_[2];
    loc[1] = params_[1]*loc[1]*factor + params_[3];
    return loc;
  }
  /// Unprojection 2D loc from image plane to point on camera plane (Z=1)
  template<typename OtherDerived>
  Vec2_t UnProject(const Eigen::MatrixBase<OtherDerived>& loc) const {
    Vec2_t pt;
    pt = loc;
    pt[0] = loc(0, 0);
    pt[1] = loc(1, 0);
    pt[0] = (pt[0] - params_[2]) / params_[0];
    pt[1] = (pt[1] - params_[3]) / params_[1];
    /// employ libcvd's implementation
    /// first guess
    Scalar scale = pt.squaredNorm();

    /// iterations of Newton-Rapheson
    for (int i = 0; i < 3 ; ++i) {  /// 3 iteration is good enough in certain range
      Scalar temp = 1 + scale*(params_[4] + params_[5]*scale);
      Scalar error = pt.squaredNorm() - scale*temp*temp;
      Scalar deriv = temp*(temp+2*scale*(params_[4] + 2*params_[5]*scale));
      scale += error/deriv;
    }
    pt = pt/static_cast<Scalar>((1 + scale*(params_[4] + params_[5]*scale)));
    return pt;
  }

  /// Project linearly (only apply PinHole part, no distortion)
  /// 3D point into image plane
  /// @param  pt[in] 3D point
  /// @return projected 2D point without distortion
  template<typename OtherDerived>
  Vec2_t LinearProject(const Eigen::MatrixBase<OtherDerived>& pt) const {
    Vec2_t loc;
    loc[0] = params_[0]*pt[0] + params_[2];
    loc[1] = params_[1]*pt[1] + params_[3];
    return loc;
  }

  /// Unproject an undistorting 2D into 3D point in image plane with z=1
  /// @param vIn2DPoint undistorting 2D point
  /// @return 3D point at image plane that z=1
  template<typename OtherDerived>
  Vec2_t LinearUnProject(const Eigen::MatrixBase<OtherDerived>& loc) const {
    Vec2_t pt;
    /// Transfer from Pixel coordinate to (u,v) coordinate
    pt[0] = loc[0];
    pt[1] = loc[1];
    pt[0] = (pt[0] - params_[2]) / params_[0];
    pt[1] = (pt[1] - params_[3]) / params_[1];
    return pt;
  }

  Eigen::Matrix<Scalar, 2, 2> GetProjectionDerivatives(const Vec2_t& pt) const {
    Eigen::Matrix<Scalar, 2, 2> derivs = Eigen::Matrix<Scalar, 2, 2>::Identity();
    Scalar r2 = pt.squaredNorm();
    Scalar t = params_[5]*r2;
    derivs *= 1 + r2*(params_[4] + t);
    derivs += (2*(params_[4] + 2*t)*pt) * pt.transpose();
    derivs.row(0) *= params_[0];
    derivs.row(1) *= params_[1];
    return derivs;
  }

  Scalar unit_pixel_distance() {
    return unit_pixel_distance_;
  }

  /// return (width,height) resolution of the camera
  void get_resolution(int& width, int& height) const {
    width = image_size_[0];
    height = image_size_[1];
  }

  /// reset image resolution for this camera model
  void set_resolution(int w, int h) {
    int old_w = image_size_[0], old_h = image_size_[1];
    image_size_[0] = w;
    image_size_[1] = h;
    params_[0] *= static_cast<SlickScalar>(w)/old_w;
    params_[1] *= static_cast<SlickScalar>(h)/old_h;
    params_[2] *= static_cast<SlickScalar>(w)/old_w;
    params_[3] *= static_cast<SlickScalar>(h)/old_h;
  }

  Vec2_t ImageSize() const {
    return image_size_;
    return Vec2_t(image_size_[0], image_size_[1]);
  }

  /// reset image resolution for this camera model
  void SetImageSize(const Vec2_t& size) {
    image_size_[0] = size[0];
    image_size_[1] = size[1];
    int old_w = image_size_[0], old_h = image_size_[1];
    image_size_[0] = size[0];
    image_size_[1] = size[1];
    params_[0] *= static_cast<SlickScalar>(size[0]) / old_w;
    params_[1] *= static_cast<SlickScalar>(size[1]) / old_h;
    params_[2] *= static_cast<SlickScalar>(size[0]) / old_w;
    params_[3] *= static_cast<SlickScalar>(size[1]) / old_h;
  }

  // For Calibration purposes --------------------------------------------------
  Eigen::Matrix<Scalar, 2, PoliCamera<Scalar>::param_n_> GetParameterDerivs(
      const Vec2_t& pt) const {
    Eigen::Matrix<Scalar, 2, PoliCamera<Scalar>::param_n_> result;
    Scalar r2 = pt.squaredNorm();
    Scalar r4 = r2 * r2;
    Vec2_t mod_camframe =
        pt * (1+ r2 * (params_[4] + r2 * params_[5]));

    result(0, 0) = mod_camframe[0];
    result(0, 1) = 0;
    result(0, 2) = 1;
    result(0, 3) = 0;
    result(0, 4) = params_[0]*pt[0]*r2;
    result(0, 5) = params_[0]*pt[0]*r4;

    result(1, 0) = 0;
    result(1, 1) = mod_camframe[1];
    result(1, 2) = 0;
    result(1, 3) = 1;
    result(1, 4) = params_[1]*pt[1]*r2;
    result(1, 5) = params_[1]*pt[1]*r4;
    return result;
  }

  const Eigen::Matrix<Scalar, PoliCamera<Scalar>::param_n_, 1>& parameters() const {
    return params_;
  }

  void SetParameters(
    const Eigen::Matrix<Scalar, PoliCamera<Scalar>::param_n_, 1>& params) {
    params_ = params;
  }
  
  void set_resolution_and_parameters(
    const int w, const int h, const Eigen::Matrix<Scalar, param_n_, 1>& params) {
    image_size_[0] = w;
    image_size_[1] = h;
    params_ = params;
  }

 protected:/// common method for initializing internal params
  void Init(const Eigen::VectorXf& params) {
    image_size_ = params.head<2>().cast<Scalar>();
    params_ = params.segment<6>(2).cast<Scalar>();
  }

  Vec2_t image_size_;  ///< [width, height] of video image 
  /// fx, fy, cx, cy, k1, k2
  Eigen::Matrix<Scalar, PoliCamera::param_n_, 1> params_;
};
}  // end namespace slick