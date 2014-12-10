// Copyright 2014 The Slick Authors. All rights reserved.
// I'm camera man, in mathematical aspect.
// You can know me more through references & description below.
// Welcome suggestions & comments about camera's models & calibration techniques.
#pragma once
#include <memory>
#include <Eigen/Eigen>
#include "slick/util/common.h"
#include "slick/slick_api.h"

namespace slick {

/// Atan Camera model is Pin Hole projection model
/// applying atan for radial distortion
/// (u,v) = ProjectionMatrix*(R_/R)*(x,y)
/// Pixel coordinate (u,v) start from Left-Top corner,
/// with vector u is along hirozontal left-right, vector v is vertical top-down
/// Responsibilities: Projecttion(3D to 2D) &
/// Unprojection(2D to 3D, at camera plane where Zc = 1.)
template<typename Scalar = SlickScalar>
class SLICK_API AtanCamera {
  typedef Eigen::Matrix<Scalar, 2, 1> Vec2_t;
  typedef Eigen::Matrix<Scalar, 3, 1> Vec3_t;
  typedef Eigen::Matrix<Scalar, 5, 1> Vec5_t;
  typedef Eigen::Matrix<Scalar, 8, 1> Vec8_t;
  typedef Eigen::Matrix<Scalar, 2, 2> Mat2_t;
 public:
  typedef Scalar ScalarType;
  enum ParamNumber { param_n_ = 5 };
  /// Camera params number (width, height, fx, fy, cx, cy, w)
  // Constructors --------------------------------------------------------------
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  /// explicit default constructor
  AtanCamera() {
    Vec8_t v8;
    v8 << 640, 480, 500, 500, 320, 240, 0;
    Init(v8);
  }

  /// convenient constructor
  AtanCamera(AtanCamera const &cam_in);

  /// constructor with width, height, & camera params
  AtanCamera(const Eigen::VectorXf& params);

  /// Key Methods --------------------------------------------------------------
  /// Do projection.
  /// @param v3Cam[in] 3D position in camera coordinate frame
  /// @return 2D point on image plane (on screen)
  template<typename OtherDerived>
  Vec2_t Project(const Eigen::MatrixBase<OtherDerived>& v2_cam) const;

  /// Do unprojection
  /// @param loc[in] 2D point on image plane (on screen)
  /// @return 2D point on camera plane (Zc=1)
  template<typename OtherDerived>
  Vec2_t UnProject(const Eigen::MatrixBase<OtherDerived>& loc) const;

  /// Project linearly (only apply PinHole part, no distortion)
  /// 3D point into image plane
  /// @param  v2_cam[in] 3D point
  /// @return projected 2D point without distortion
  template<typename OtherDerived>
  Vec2_t ProjectLinear(const Eigen::MatrixBase<OtherDerived>& v2_cam) const;

  /// Unproject an undistorting 2D into 3D point in image plane with z=1
  /// @param vIn2DPoint undistorting 2D point
  /// @return 3D point at image plane that z=1
  template<typename OtherDerived>
  Vec2_t UnProjectLinear(const Eigen::MatrixBase<OtherDerived>& vIn2DPoint) const;

  Mat2_t GetProjectionDerivatives(const Vec2_t& v2_cam) const;

  Scalar unit_pixel() { return unit_pixel_distance_; }

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
  Eigen::Matrix<Scalar, 2, 5> GetParameterDerivs(const Vec2_t& v2_cam) {
    Eigen::Matrix<Scalar, 2, 5> jacobian;
    /// TODO
    Vec5_t params = params_;
    auto loc = Project(v2_cam);
    for (int i = 0; i < 5; ++i) {
      if (i == 4 && params_[4] == 0.0)
        continue;
      Vec5_t up = Vec5_t::Zero();
      up[i] += 0.001;
      params_ += up;
      jacobian.col(i) = Project(v2_cam - loc) / 0.001;
      params_ = params;
    }
    if (params_[4] == 0.0)
      jacobian.col(4).setZero();
    return jacobian;
  }

  const Vec5_t& parameters() const { return params_; }
  void SetParameters(const Vec5_t& params);

 protected:
  /// common method for initializing internal params
   void Init(const Eigen::VectorXf& params);
  // Radial factor : distorted / undistorted radius.
  Scalar radial_factor(Scalar r) const {
    if (r < 0.001 || params_[4] == 0.0)
      return 1.0;
    else {
      const Scalar d = 2.0 * tan(params_[4] / 2.0);
      return (1.0f / params_[4] * atan(r * d) / r);
    }
  }

  Scalar inverse_radial_factor(Scalar r) const {
    if (params_[4] == 0.0)
      return r;
    const Scalar d = 2.0 * tan(params_[4] / 2.0);

    return  tan(r * params_[4]) / d;
  }


  Vec2_t image_size_;  ///< the width, height of the camera image
  /// Camera's instrinsic parameters
  Vec5_t params_;  ///< fx, fy, cx, cy, k1, k2
};

typedef std::shared_ptr<AtanCamera<double>> AtanCameradPtr;
typedef std::shared_ptr<AtanCamera<float>> AtanCamerafPtr;
typedef std::shared_ptr<const AtanCamera<double>> AtanCameradConstPtr;
typedef std::shared_ptr<const AtanCamera<float>> AtanCamerafConstPtr;
// Implementation ==============================================================
/// const int iParamN = AtanCamera<>::iParamN;
/// explicit declared to advoide compiler's specific configuration

template<typename Scalar>
template<typename OtherDerived>
inline Eigen::Matrix<Scalar, 2, 1> AtanCamera<Scalar>::Project(
    const Eigen::MatrixBase<OtherDerived>& v2_cam) const {
  const Scalar r = v2_cam.norm();
  const Scalar rf = radial_factor(r);
  const Scalar dist_r = rf * r;
  auto dist_cam = rf * v2_cam;
  Vec2_t loc;
  loc[0] = params_[2] + params_[0] * dist_cam[0];
  loc[1] = params_[3] + params_[1] * dist_cam[1];
  return loc;
}

template<typename Scalar>
template<typename OtherDerived>
inline Eigen::Matrix<Scalar, 2, 1> AtanCamera<Scalar>::UnProject(
    const Eigen::MatrixBase<OtherDerived>& loc) const {
  Vec2_t v2_cam;
  v2_cam[0] = (loc[0] - params_[2]) / params_[0];
  v2_cam[1] = (loc[1] - params_[3]) / params_[1];
  const Scalar r = v2_cam.norm();
  const Scalar r1 = inverse_radial_factor(r);
  Scalar factor = (r > 0.01) ? r1 / r : 1.0;
  return factor * v2_cam;
}

template<typename Scalar>
template<typename OtherDerived>
inline Eigen::Matrix<Scalar, 2, 1> AtanCamera<Scalar>::ProjectLinear(
    const Eigen::MatrixBase<OtherDerived>& v2_cam) const {
  Vec2_t loc;
  loc[0] = params_[0]*v2_cam[0] + params_[2];
  loc[1] = params_[1]*v2_cam[1] + params_[3];
  return loc;
}

template<typename Scalar>
template<typename OtherDerived>
inline Eigen::Matrix<Scalar, 2, 1> AtanCamera<Scalar>::UnProjectLinear(
    const Eigen::MatrixBase<OtherDerived>& loc) const {
  Vec2_t v2_cam;
  /// Transfer from Pixel coordinate to (u,v) coordinate
  v2_cam[0] = loc[0];
  v2_cam[1] = loc[1];
  v2_cam[0] = (v2_cam[0] - params_[2]) / params_[0];
  v2_cam[1] = (v2_cam[1] - params_[3]) / params_[1];
  return v2_cam;
}
}  // end namespace slick