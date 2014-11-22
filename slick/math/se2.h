// Copyright 2012, 2013, 2014 The Look3D Authors. All rights reserved.
#ifndef LOOK3D_MATH_SE2_H_
#define LOOK3D_MATH_SE2_H_
#include "math/common.h"
#include "math/so2.h"
#include "math/look3d_math_api.h"

namespace look3d {

template <typename Precision>
class SE2Group;

typedef SE2Group<DefaultScalarType> SE2;
typedef SE2Group<double> SE2d;
typedef SE2Group<float> SE2f;

/// Represent a two-dimensional Euclidean transformation (rotation &
/// translation)
/// This transformation is a member of the Special Euclidean Lie group SE2Group.
/// These can be parameterised with
/// three numbers (in the space of the Lie Algebra).
/// In this class, the first two parameters are a translation vector
/// while the third is the amount of rotation in the plane as for SO2Group.
template <typename Precision = DefaultScalarType>
class SE2Group {
  typedef Eigen::Matrix<Precision, 2, 3> MatrixType;

 public:
  // Constructors  =============================================================
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  /// initialises the rotation to zero(the identity) and the translation to zero
  SE2Group() { translation_ = Eigen::Matrix<Precision, 2, 1>::Zero(); }
  template <int O>
  SE2Group(const SO2Group<Precision>& R,
           const Eigen::Matrix<Precision, 2, 1, O>& T) {
    rotation_ = R;
    translation_ = T;
  }

  template <class P, int S, int O>
  SE2Group(const Eigen::Matrix<P, S, 1, O>& v) {
    *this = this->exp(v);
  }

  // Operators  ================================================================
  /// Right-multiply by another SE2Group (concatenate the two transformations)
  /// @param rhs The multipier
  template <typename P>
  SE2Group<typename Eigen::internal::scalar_product_traits<
      Precision, P>::ReturnType> operator*(const SE2Group<P>& rhs) const {
    Eigen::Matrix<Precision, 2, 1> v2 = translation_ + rotation_*rhs.get_translation();
    return SE2Group<typename Eigen::internal::scalar_product_traits<
        Precision, P>::ReturnType>(rotation_ * rhs.get_rotation(), v2);
  }

  /// Self right-multiply by another SE2Group (concatenate the two
  /// transformations)
  /// @param rhs The multipier
  inline SE2Group& operator*=(const SE2Group& rhs) {
    *this = *this * rhs;
    return *this;
  }

  /// Right-multiply by another matrix
  template <typename OtherDerived>
  typename Eigen::Matrix<Precision, OtherDerived::RowsAtCompileTime,
                         OtherDerived::ColsAtCompileTime> operator*(
      const Eigen::MatrixBase<OtherDerived>& other) const {
    Eigen::Matrix<Precision, OtherDerived::RowsAtCompileTime,
                  OtherDerived::ColsAtCompileTime> res = other;
    if (other.rows() == 3) {
      res.block(0, 0, 2, other.cols()) =
          rotation_ * other.block(0, 0, 2, other.cols()) +
          translation_ * other.block(2, 0, 1, other.cols());
      res.block(2, 0, 1, other.cols()) = other.block(2, 0, 1, other.cols());
    } else {
      res.block(0, 0, 2, other.cols()) =
          rotation_ * other.block(0, 0, 2, other.cols()) +
          translation_;
    }
    return res;
  }

  // SE2Group Specific  ========================================================
  /// Exponentiate a Vector in the Lie Algebra to generate a new SE2Group.
  /// See the Detailed Description for details of this vector.
  /// @param vect The Vector to exponentiate
  template <typename P, int O>
  static inline SE2Group exp(const Eigen::Matrix<P, 3, 1, O>& vect);

  /// Take the logarithm of the matrix,
  /// generating the corresponding vector in the Lie Algebra.
  /// See the Detailed Description for details of this vector.
  static inline Eigen::Matrix<Precision, 3, 1> ln(const SE2Group& SE2Group);
  // @overload
  Eigen::Matrix<Precision, 3, 1> ln() const { return SE2Group::ln(*this); }

  /// compute the inverse of the transformation
  SE2Group inverse() const {
    const SO2Group<Precision>& rinv = this->rotation_.inverse();
    Eigen::Matrix<Precision, 2, 1> v2 = -(rinv * translation_);
    return SE2Group(rinv, v2);
  }

  /// returns the generators for the Lie group. These are a set of matrices that
  /// form a basis for the vector space of the Lie algebra.
  /// - 0 is translation in x
  /// - 1 is translation in y
  /// - 2 is rotation in the plane
  static inline Eigen::Matrix<Precision, 3, 3> generator(int i) {
    Eigen::Matrix<Precision, 3, 3> result(
        Eigen::Matrix<Precision, 3, 3>::Zero());
    if (i < 2) {
      result(i, 2) = 1;
      return result;
    }
    result(0, 1) = -1;
    result(1, 0) = 1;
    return result;
  }

  /// transfers a vector in the Lie algebra, from one coord frame to another
  /// so that exp(adjoint(vect)) = (*this) * exp(vect) * (this->inverse())
  template <int Options>
  Eigen::Matrix<Precision, 3, 1> adjoint(
      const Eigen::Matrix<Precision, 3, 1, Options>& vect) const {
    Eigen::Matrix<Precision, 3, 1> result;
    result[2] = vect[2];
    result.template block<0, 2>() =
        this->block<2, 2>(0, 0) * vect.template segment<2>(0);
    result[0] += vect[2] * this(1, 2);
    result[1] -= vect[2] * this(0, 2);
    return result;
  }

  template <int Options>
  Eigen::Matrix<Precision, 3, 3> adjoint(
      const Eigen::Matrix<Precision, 3, 3, Options>& M) const {
    Eigen::Matrix<Precision, 3, 3> result;
    for (int i = 0; i < 3; ++i) result.col(i) = M.col(i).adjoint();
    for (int i = 0; i < 3; ++i) result.row(i) = result.row(i).adjoint();
    return result;
  }

  const SO2Group<Precision>& get_rotation() const { return rotation_; }

  const Eigen::Matrix<Precision, 2, 1>& get_translation() const {
    return translation_;
  }

  Eigen::Matrix<Precision, 2, 3> get_matrix() const {
    Eigen::Matrix<Precision, 2, 3> res;
    res.template block<2, 2>(0 ,0) = rotation_.get_matrix();
    res.template block<2, 1>(0, 2) = translation_;
    return res;
  }

 protected:
  SO2Group<Precision> rotation_;
  Eigen::Matrix<Precision, 2, 1> translation_;
};  // end class SE2Group

template <typename Precision>
template <typename PV, int O>
inline SE2Group<Precision> SE2Group<Precision>::exp(
    const Eigen::Matrix<PV, 3, 1, O>& mu) {
  static const Precision one_6th = 1.0 / 6.0;
  static const Precision one_20th = 1.0 / 20.0;

  const Precision theta = mu[2];
  const Precision theta_sq = theta * theta;

  const Eigen::Matrix<Precision, 2, 1> cross =
      Eigen::Matrix<Precision, 2, 1>(-theta * mu[1], theta * mu[0]);
  SO2Group<Precision> res_rot = SO2Group<Precision>::exp(theta);
  Eigen::Matrix<Precision, 2, 1> res_t;

  if (theta_sq < 1e-8) {
    res_t = mu.template segment<2>(0) + 0.5 * cross;
  } else {
    Precision A, B;
    if (theta_sq < 1e-6) {
      A = 1.0 - theta_sq * one_6th * (1.0 - one_20th * theta_sq);
      B = 0.5 - 0.25 * one_6th * theta_sq;
    } else {
      const Precision inv_theta = (1.0 / theta);
      const Precision sine = res_rot.get_matrix()(1, 0);
      const Precision cosine = res_rot.get_matrix()(0, 0);
      A = sine * inv_theta;
      B = (1 - cosine) * (inv_theta * inv_theta);
    }
    res_t = A * mu.template segment<2>(0) + B * cross;
  }
  return SE2Group<Precision>(res_rot, res_t);
}

template <typename Precision>
inline Eigen::Matrix<Precision, 3, 1> SE2Group<Precision>::ln(
    const SE2Group<Precision>& SE2Group) {
  const Precision theta = SE2Group.get_rotation().ln();

  Precision shtot = 0.5;
  if (fabs(theta) > 0.00001) shtot = sin(theta / 2) / theta;

  const SO2Group<Precision> halfrotator(theta * -0.5);
  Eigen::Matrix<Precision, 3, 1> result;
  result.template segment<2>(0) =
      (halfrotator * SE2Group.get_translation()) / (2 * shtot);
  result[2] = theta;
  return result;
}

// External Operators  =========================================================
/// Write an SE2Group to a stream
/// @relates SE2Group
template <class Precision>
inline std::ostream& operator<<(std::ostream& os,
                                const SE2Group<Precision>& rhs) {
  std::streamsize fw = os.width();
  for (int i = 0; i < 2; i++) {
    os.width(fw);
    os << rhs.get_rotation().get_matrix().row(i);
    os.width(fw);
    os << " " << rhs.get_translation()[i] << std::endl;
  }
  return os;
}

/// Read an SE2Group from a stream
/// @relates SE2Group
template <class Precision>
    inline std::istream& operator>>(std::istream& is,
                                    SE2Group<Precision>& rhs) {
  for (int i = 0; i < 2; i++)
    is >> const_cast<Eigen::Matrix<Precision, 2, 2> >(
              rhs.get_rotation().get_matrix()).row(i) >>
        const_cast<Eigen::Matrix<Precision, 2, 1> >(rhs.get_translation())[i];
  const_cast<SO2Group<Precision> >(rhs.get_rotation()).coerce();
  return is;
}

// Left-multiply by a Matrix
// @relates SO2Group
template <typename OtherDerived, typename Precision>
inline Eigen::Matrix<typename Eigen::internal::scalar_product_traits<
                         typename OtherDerived::Scalar, Precision>::ReturnType,
                     OtherDerived::RowsAtCompileTime, 3> operator*(
    const Eigen::MatrixBase<OtherDerived>& lhs,
    const SE2Group<Precision>& rhs) {
  Eigen::Matrix<typename Eigen::internal::scalar_product_traits<
                    typename OtherDerived::Scalar, Precision>::ReturnType,
                OtherDerived::RowsAtCompileTime, 3> res;
  res.block(0, 0, lhs.rows(), 2) =
      lhs.block(0, 0, lhs.rows(), 2) * rhs.get_rotation();
  res.block(0, 2, lhs.rows(), 1) =
      lhs.block(0, 2, lhs.rows(), 1) +
      lhs.block(0, 0, lhs.rows(), 2) * rhs.get_translation();
  return res;
}

/// Multiply a SO2Group with and SE2Group
/// @relates SE2Group
/// @relates SO2Group
template <typename Precision1, typename Precision2>
inline SE2Group<typename Eigen::internal::scalar_product_traits<
    Precision1, Precision2>::ReturnType> operator*(
    const SO2Group<Precision1>& lhs, const SE2Group<Precision2>& rhs) {
  Eigen::Matrix<typename Eigen::internal::scalar_product_traits<
                    Precision1, Precision2>::ReturnType,
                2, 1> v2 = lhs * rhs.get_translation();
  return SE2Group<typename Eigen::internal::scalar_product_traits<
      Precision1, Precision2>::ReturnType>(lhs * rhs.get_rotation(), v2);
}
}       // namespace look3d
#endif  // LOOK3D_MATH_SE2_H_
