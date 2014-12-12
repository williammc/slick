// Copyright 2014 The Slick Authors. All rights reserved.
// I'm camera man, in mathematical aspect.
// You can know me more through references & description below.
// Welcome suggestions & comments about camera's models & calibration techniques.
#pragma once
#include <memory>
#include <Eigen/Eigen>
#include "slick/util/common.h"

namespace slick {

/// Atan Camera model is Pin Hole projection model
/// applying atan for radial distortion
/// @ref: Devernay and Faugeras, "Straight lines have to be straight" for fish-eye lenses.
/// Pixel coordinate (u,v) start from Left-Top corner,
/// with vector u is along hirozontal left-right, vector v is vertical top-down
/// Responsibilities: Projecttion(3D to 2D) &
/// Unprojection(2D to 3D, at camera plane where Zc = 1.)
template<typename Scalar = SlickScalar>
class AtanCamera {
  typedef Eigen::Matrix<Scalar, 2, 1> Vec2_t;
  typedef Eigen::Matrix<Scalar, 3, 1> Vec3_t;
  typedef Eigen::Matrix<Scalar, 5, 1> Vec5_t;
  typedef Eigen::Matrix<Scalar, 8, 1> Vec8_t;
  typedef Eigen::Matrix<Scalar, 2, 2> Mat2_t;
public:
  typedef std::shared_ptr<AtanCamera> Ptr;
  typedef Scalar ScalarType;
  enum ParamNumber { param_n_ = 5 };
  /// Camera params number (width, height, fx, fy, cx, cy, w)
  // Constructors --------------------------------------------------------------
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  /// explicit default constructor
  AtanCamera() {
    Eigen::VectorXf v(7);
    v << 640, 480, 500, 500, 320, 240, 0;
    Init(v);
  }

  /// convenient constructor
  AtanCamera(AtanCamera const &cam_in) {
    image_size_ = cam_in.image_size_;
    params_ = cam_in.params_;
  }

  /// constructor with width, height, & camera params
  AtanCamera(const Eigen::VectorXf& params) {
    Init(params);
  }

  /// Key Methods --------------------------------------------------------------
  /// Projection of @v_cam (2D or 3D point on camera plane) to image plane.
  template<typename OtherDerived>
  Vec2_t Project(const Eigen::MatrixBase<OtherDerived>& v_cam) const {
    return LinearProject(radial_distort(v_cam) * v_cam);
  }

  /// Unprojection of 2D point on image plane to 3D point on camera plane (Z=1)
  template<typename OtherDerived>
  Vec2_t UnProject(const Eigen::MatrixBase<OtherDerived>& loc) const {
    Vec2_t v2_cam = LinearUnProject(loc);
    auto const r2 = v2_cam.dot(v2_cam);
    Scalar w2 = params_[4] * params_[4];
    Scalar k3 = -w2 / 3.0;
    Scalar k5 = w2*w2 / 5.0;
    Scalar k7 = -w2*w2*w2 / 7.0;

    // 3 iterations of Newton-Raphson
    auto scale = r2;
    for (int i = 0; i < 3; ++i) {
      Scalar t = 1 + scale * (k3 + scale*(k5 + scale*k7));
      Scalar error = r2 - scale*t*t;
      Scalar deriv = t*(t + 2*scale*(k3 + 2*scale*(k5 + 1.5*scale*k7)));
      scale += error / deriv;
    }
    
    v2_cam = v2_cam / (1 + scale*(k3 + scale*(k5 + scale*k7)));
    return v2_cam;
  }

  /// Linear projection of 2D or 3D point on camera plane to image plane
  template<typename OtherDerived>
  Vec2_t LinearProject(const Eigen::MatrixBase<OtherDerived>& v_cam) const {
    Vec2_t loc;
    loc[0] = params_[0]*v_cam[0] + params_[2];
    loc[1] = params_[1]*v_cam[1] + params_[3];
    return loc;
  }

  /// Unproject an undistorting 2D @loc into 3D point in image plane with z=1
  template<typename OtherDerived>
  Vec2_t LinearUnProject(const Eigen::MatrixBase<OtherDerived>& loc) const {
    Vec2_t v2_cam;
    /// Transfer from Pixel coordinate to (u,v) coordinate
    v2_cam[0] = (loc[0] - params_[2]) / params_[0];
    v2_cam[1] = (loc[1] - params_[3]) / params_[1];
    return v2_cam;
  }

  /// Derivative w.r.t camera frame at point @pt
  Mat2_t GetProjectionDerivatives(const Vec2_t& pt) const {
    Eigen::Matrix<Scalar, 2, 2> derivs = Eigen::Matrix<Scalar, 2, 2>::Identity();
    auto const r2 = pt.dot(pt);
    Scalar w2 = params_[4] * params_[4];
    Scalar k3 = -w2 / 3.0;
    Scalar k5 = w2*w2 / 5.0;
    Scalar k7 = -w2*w2*w2 / 7.0;
    derivs *= (1 + k3*r2 + k5*r2*r2 + k7*r2*r2*r2);
    derivs += ((2*k3 + 4*k5*r2 + 6*k7*r2*r2) * pt) * pt.transpose();
    derivs.row(0) *= params_[0];
    derivs.row(1) *= params_[1];
    return derivs;
  }

  /// return (width,height) resolution of the camera
  Vec2_t ImageSize() const { return image_size_; }

  /// reset image resolution for this camera model
  void SetImageSize(const Vec2_t& size) { 
    auto old_size = image_size_;
    image_size_ = size; 
    params_[0] *= image_size_[0] / old_size[0];
    params_[1] *= image_size_[1] / old_size[1];
    params_[2] *= image_size_[0] / old_size[0];
    params_[3] *= image_size_[1] / old_size[1];
  }

  void SetImageSizeAndParameters(const Vec8_t& v8) {
    Init(v8);
  }

  // For Calibration purposes --------------------------------------------------
  /// Employing numerical automatic differenciation (George Klein's approach)
  Eigen::Matrix<Scalar, 2, 5> GetParameterDerivs(const Vec2_t& pt) {
    auto const v2 = radial_distort(pt)*pt;
    auto const r2 = pt.dot(pt);
    auto const r4 = r2*r2;
    auto const r6 = r2*r4;

    auto const w = params_[4];
    auto const w3 = w*w*w;
    auto const w5 = w*w*w3;

    auto const k1 = -(2.0/3.0)*w*r2;
    auto const k2 = (4.0/5.0)*w3*r4;
    auto const k3 = -(6.0/7.0)*w5*r6;

    Eigen::Matrix<Scalar, 2, 5> derivs;
    derivs(0, 0) = v2[0];
    derivs(0, 1) = 0;
    derivs(0, 2) = 1;
    derivs(0, 3) = 0;
    derivs(0, 4) = params_[0] * (k1 + k2 + k3) * pt[0];

    derivs(1, 0) = 0;
    derivs(1, 1) = v2[1];
    derivs(1, 2) = 0;
    derivs(1, 3) = 1;
    derivs(1, 4) = params_[1] * (k1 + k2 + k3) * pt[1];
    return derivs;
  }

  const Vec5_t& parameters() const { return params_; }
  void SetParameters(const Vec5_t& params) {
    params_ = params;
  }

 protected:
  /// common method for initializing internal params
  void Init(const Eigen::VectorXf& params) {
    image_size_ = params.head<2>().cast<Scalar>();
    params_ = params.segment<5>(2).cast<Scalar>();
  }

  /// radial factor : distorted / undistorted radius
  Scalar radial_distort(const Vec2_t& v_cam) const {
    auto const r2 = v_cam.dot(v_cam);
    auto const w2 = params_[4]*params_[4];
    auto const fac = w2 * r2;
    Scalar cons[3] = {-1.0/3.0, 1.0/5.0, -1.0/7.0};
    Scalar term = 1.0;
    Scalar scale = term;
    for (int i = 0; i < 3; ++i) {
      term *= fac;
      scale += term * cons[i];
    }
    return scale;
  }

  Vec2_t image_size_;  ///< the width, height of the camera image
  /// Camera's instrinsic parameters
  Vec5_t params_;  ///< fx, fy, cx, cy, k1, k2
};

typedef std::shared_ptr<AtanCamera<double>> AtanCameradPtr;
typedef std::shared_ptr<AtanCamera<float>> AtanCamerafPtr;
typedef std::shared_ptr<const AtanCamera<double>> AtanCameradConstPtr;
typedef std::shared_ptr<const AtanCamera<float>> AtanCamerafConstPtr;
}  // end namespace slick