// Copyright 2014 The Slick Authors. All rights reserved.
// I'm camera man, in mathematical aspect.
// You can know me more through references & description below.
// Welcome suggestions & comments about camera's models & calibration techniques.
#pragma once
#include <memory>
#include <Eigen/Eigen>
#include "slick/math/se3.h"
#include "slick/util/common.h"

namespace slick {

// Cylindrical panorama camera used to project images onto a cylindrical
// panorama geometry.
// seam is at -Z
// image X is in positive rotation direction
template<typename Scalar = float>
class Cylinder {
public:
  Cylinder (const Eigen::Matrix<Scalar, 2, 1> & size,
            const Scalar & height_angle) {
    set_parameters( size, height_angle );
  }

  Cylinder () {
    set_parameters( Eigen::Matrix<Scalar, 2, 1>(1024.0, 256.0), 1.0 );
  }

  void set_parameters(const Eigen::Matrix<Scalar, 2, 1> & s, 
                      const Scalar & height_angle) {
    size_ = s;
    tan_height_angle_ = tan(height_angle*M_PI/180/2);
    twopi_inverse_ = 1.0/(2.0*M_PI);
    inv_tan_heightangle_= 1.0/(2.0*tan_height_angle_);
  }

  Eigen::Matrix<Scalar, 2, 1> project(const Eigen::Matrix<Scalar, 3, 1> & v) const {
    const Eigen::Matrix<Scalar, 3, 1> dir = v.normalized();

    const Scalar rot = atan2(dir[0],dir[2]);
    const Scalar up = asin(dir[1]);
    Eigen::Matrix<Scalar, 2, 1> nMapCoord;

    nMapCoord[0] = (rot + M_PI)*twopi_inverse_*size_[0];
    nMapCoord[1] = tan(up)*inv_tan_heightangle_*size_[1]+size_[1]*0.5;
    return nMapCoord;
  }

  Eigen::Matrix<Scalar, 3, 1> unproject(const Eigen::Matrix<Scalar, 2, 1> & p) const {
    const Scalar rot = p[0] * 2 * M_PI / size_[0] - M_PI;
    const Scalar up = atan((p[1] - size_[1] * 0.5) * 2 * tan_height_angle_ / size_[1]);
    return Eigen::Matrix<Scalar, 3, 1>(cos(up) * sin(rot), sin(up), cos(up) * cos(rot));
  }

  const Eigen::Matrix<Scalar, 2, 1> & size()const { return size_; }

protected:
  // image size
  Eigen::Matrix<Scalar, 2, 1> size_;
  // tan of the angle out of the horizon plane
  Scalar tan_height_angle_;
  Scalar twopi_inverse_;
  Scalar inv_tan_heightangle_;
};

}  // namespace slick