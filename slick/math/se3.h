// Copyright 2014 The Slick Authors. All rights reserved.
#pragma once
#include "slick/math/so3.h"
#include "slick/slick_api.h"

namespace slick {

template <typename Precision> class SE3Group;

typedef SE3Group<SlickScalar> SE3;
typedef SE3Group<double> SE3d;
typedef SE3Group<float> SE3f;

/// Represent a three-dimensional Euclidean transformation (a rotation and a translation).
/// This can be represented by a \f$3\times\f$4 matrix operating on a homogeneous co-ordinate,
/// so that a vector \f$\underline{x}\f$ is transformed to a new location \f$\underline{x}'\f$
/// by
/// \f[\begin{aligned}\underline{x}' &= E\times\underline{x}\\ \begin{bmatrix}x'\\y'\\z'\end{bmatrix} &= \begin{pmatrix}r_{11} & r_{12} & r_{13} & t_1\\r_{21} & r_{22} & r_{23} & t_2\\r_{31} & r_{32} & r_{33} & t_3\end{pmatrix}\begin{bmatrix}x\\y\\z\\1\end{bmatrix}\end{aligned}\f]
///
/// This transformation is a member of the Special Euclidean Lie group SE3Group. These can be parameterised
/// six numbers (in the space of the Lie Algebra). In this class, the first three parameters are a
/// translation vector while the second three are a rotation vector, whose direction is the axis of rotation
/// and length the amount of rotation (in radians), as for SO3Group
/// @ingroup math
template <typename Precision = SlickScalar>
class SLICK_API SE3Group {
  typedef Eigen::Matrix<Precision, 4, 4> MatrixType;

 public:
  // Constructors  =============================================================
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  /// Default constructor. Initialises the the rotation to zero (the identity) and the translation to zero
  SE3Group() { translation_ = Eigen::Matrix<Precision, 3, 1>::Zero(); }

  template <typename Deri31>
  SE3Group(const SO3Group<typename Deri31::Scalar> & R,
      const Eigen::MatrixBase<Deri31>& T) {
    rotation_ = R;
    translation_ = T;
  }

  /// v = tx, ty, tz, rx, ry, rz
  template <typename DeriVector6>
  explicit SE3Group(const Eigen::MatrixBase<DeriVector6> & v) { *this = exp(v); }

  // Operators  ================================================================
  /// Right-multiply by another SE3Group (concatenate the two transformations)
  /// @param rhs The multipier
  SE3Group<Precision>& operator *= (const SE3Group& rhs) {
    translation_ += get_rotation() * rhs.get_translation();
    rotation_ *= rhs.get_rotation();
    return *this;
  }

  /// Right-multiply by another SE3Group (concatenate the two transformations)
  /// @param rhs The multipier
  template<typename P>
  SE3Group<typename Eigen::internal::scalar_product_traits<Precision, P>::ReturnType>
    operator * (const SE3Group<P>& rhs) const {
    SO3Group<typename Eigen::internal::scalar_product_traits<Precision, P>::ReturnType>
        SO3GroupTemp = rotation_*rhs.get_rotation();
    Eigen::Matrix<Precision, 3, 1> v3Temp =
        translation_ + rotation_*rhs.get_translation();
    return SE3Group<typename Eigen::internal::scalar_product_traits<Precision, P>::ReturnType>(SO3GroupTemp, v3Temp);
  }

  /// Required for Eigen inheritance
  template<typename OtherDerived>
  SE3Group& operator = (const Eigen::MatrixBase <OtherDerived>& other) {
    rotation_ = other.template block<3, 3>(0, 0);  /// SO3Group object is constructed
    translation_ = other.template block<3, 1>(0, 3);
    return *this;
  }

  /// Right-multiply by another matrix
  template<typename OtherDerived>
  typename Eigen::Matrix<Precision, OtherDerived::RowsAtCompileTime, OtherDerived::ColsAtCompileTime> operator *(const Eigen::MatrixBase <OtherDerived>& other) const {
    Eigen::Matrix<Precision, OtherDerived::RowsAtCompileTime, OtherDerived::ColsAtCompileTime> res = other;
    if (other.rows() == 4) {
      res.block(0, 0, 3, other.cols()) =
          rotation_*other.block(0, 0, 3, other.cols()) +
          translation_*other.block(3, 0, 1, other.cols());
      res.block(3, 0, 1, other.cols()) = other.block(3, 0, 1, other.cols());
    } else {
      res.block(0, 0, 3, other.cols()) =
          rotation_*other.block(0, 0, 3, other.cols());
      for (int i = 0; i < other.cols(); ++i) {
        res.block(0, i, 3, 1) += translation_;
      }
    }
    return res;
  }

  // SE3Group Specific  ========================================================
  /// Exponentiate a Vector in the Lie Algebra to generate a new SE3Group.
  /// See the Detailed Description for details of this vector.
  /// @param vect The Vector to exponentiate
  template <typename OtherDerived>
  static inline SE3Group<Precision> exp(const Eigen::MatrixBase<OtherDerived>& vect);

  /// Take the logarithm of the matrix, generating the corresponding vector in the Lie Algebra.
  /// See the Detailed Description for details of this vector.
  static inline Eigen::Matrix<Precision, 6, 1> ln(const SE3Group& SE3Group);
  /// @overload
  Eigen::Matrix<Precision, 6, 1> ln() const { return SE3Group::ln(*this); }

  SE3Group<Precision> inverse() const {
    const SO3Group<Precision> rinv = this->rotation_.inverse();
    Eigen::Matrix<Precision, 3, 1> v = -(rinv*translation_);
    return SE3Group<Precision>(rinv, v);
  }

  static inline Eigen::Matrix<Precision, 4, 4> generator(int i) {
    Eigen::Matrix<Precision, 4, 4> result(
          Eigen::Matrix<Precision, 4, 4>::Zero());
    if (i < 3) {
      result(i, 3)=1;
      return result;
    }
    result((i+1)%3, (i+2)%3) = -1;
    result((i+2)%3, (i+1)%3) = 1;
    return result;
  }

  /// Returns the i-th generator times pos
  template<int Base>
  static Eigen::Matrix<Precision, 4, 1>
    generatorField(int i, const Eigen::Matrix<Precision, 4, 1, Base>& pos) {
    Eigen::Matrix<Precision, 4, 1> result(
          Eigen::Matrix<Precision, 4, 1>::Zeros());
    if (i < 3) {
      result[i]=pos[3];
      return result;
    }
    result[(i+1)%3] = - pos[(i+2)%3];
    result[(i+2)%3] = pos[(i+1)%3];
    return result;
  }

  /// Transfer a matrix in the Lie Algebra from one
  /// co-ordinate frame to another. This is the operation such that for a matrix
  /// \f$ B \f$,
  /// \f$ e^{\text{Adj}(v)} = Be^{v}B^{-1} \f$
  /// @param M The Matrix to transfer
  template<int S, typename P2, int Options>
  Eigen::Matrix<Precision, 6, 1> adjoint(
      const Eigen::Matrix<P2, S, 1, Options>& vect)const;

  /// Transfer covectors between frames (using the transpose of the inverse of the adjoint)
  /// so that trinvadjoint(vect1) * adjoint(vect2) = vect1 * vect2
  template<int S, typename P2, int Options>
  Eigen::Matrix<Precision, 6, 1> trinvadjoint(
      const Eigen::Matrix<P2, S, 1, Options>& vect)const;

  /// @overload
  template <int R, int C, typename P2, int Options>
  Eigen::Matrix<Precision, 6, 6> adjoint(
      const Eigen::Matrix<P2, R, C, Options>& M)const;

  /// @overload
  template <int R, int C, typename P2, int Options>
  Eigen::Matrix<Precision, 6, 6> trinvadjoint(
      const Eigen::Matrix<P2, R, C, Options>& M)const;

  const SO3Group<Precision>& get_rotation() const { return rotation_; }

  const Eigen::Matrix<Precision, 3, 1>& get_translation() const {
    return translation_;
  }

  Eigen::Matrix<Precision, 3, 4> get_matrix() const {
    Eigen::Matrix<Precision, 3, 4> res;
    res.template block<3, 3>(0 ,0) = rotation_.get_matrix();
    res.template block<3, 1>(0, 3) = translation_;
    return res;
  }

 protected:
  SO3Group<Precision> rotation_;
  Eigen::Matrix<Precision, 3, 1> translation_;
};  /// class SE3Group

/// transfers a vector in the Lie algebra
/// from one coord frame to another
/// so that exp(adjoint(vect)) = (*this) * exp(vect) * (this->inverse())
template<typename Precision>
template<int S, typename P2, int Options>
inline Eigen::Matrix<Precision, 6, 1> SE3Group<Precision>::
  adjoint(const Eigen::Matrix<P2, S, 1, Options>& vect) const {
  assert(S == vect.rows());
  Eigen::Matrix<Precision, 6, 1> result;
  result.template segment<3>(3) = rotation_ * vect.template segment<3>(3);
  result.template segment<3>(0) = rotation_ * vect.template segment<3>(0);
  result.template segment<3>(0) +=
      translation_.cross(result.template segment<3>(0));
  return result;
}

/// transfers covectors between frames
/// (using the transpose of the inverse of the adjoint)
/// so that trinvadjoint(vect1) * adjoint(vect2) = vect1 * vect2
template<typename Precision>
template<int S, typename P2, int Options>
inline Eigen::Matrix<Precision, 6, 1> SE3Group<Precision>::
  trinvadjoint(const Eigen::Matrix<P2, S, 1, Options>& vect) const {
  assert(S == vect.rows());
  Eigen::Matrix<Precision, 6, 1> result;
  result.template segment<3>(3) = rotation_ * vect.template segment<3>(3);
  result.template segment<3>(0) = rotation_ * vect.template segment<3>(0);
  result.template segment<3>(3) +=
      translation_.cross(result.template segment<3>(0));
  return result;
}

template<typename Precision>
template<int R, int C, typename P2, int Options>
inline Eigen::Matrix<Precision, 6, 6> SE3Group<Precision>::
  adjoint(const Eigen::Matrix<P2, R, C, Options>& M) const {
  assert(R == M.cols() && C == M.rows());

  Eigen::Matrix<Precision, 6, 6> result;
  for (int i = 0; i < 6; i++) {
    result.col(i) = adjoint(M.col(i));
  }
  for (int i = 0; i < 6; i++) {
    result.row(i) = adjoint(result.row(i));
  }
  return result;
}

template<typename Precision>
template<int R, int C, typename P2, int Options>
inline Eigen::Matrix<Precision, 6, 6> SE3Group<Precision>::
  trinvadjoint(const Eigen::Matrix<P2, R, C, Options>& M) const {
  assert(R == M.cols() && C == M.rows());
  Eigen::Matrix<Precision, 6, 6> result;
  for (int i = 0; i < 6; i++) {
    result.col(i) = trinvadjoint(M.col(i));
  }
  for (int i = 0; i < 6; i++) {
    result.row(i) = trinvadjoint(result.row(i));
  }
  return result;
}

template <typename Precision>
template <typename OtherDerived>
inline SE3Group<Precision> SE3Group<Precision>::exp(
    const Eigen::MatrixBase<OtherDerived>& mu) {
  static const Precision one_6th = 1.0/6.0;
  static const Precision one_20th = 1.0/20.0;

  Eigen::Matrix<Precision, 3, 1> vec3;

  const Eigen::Matrix<Precision, 3, 1> w = mu.template block<3, 1>(3, 0);
  const Precision theta_sq = w.transpose()*w;
  const Precision theta = sqrt(theta_sq);
  Precision A, B;

  const Eigen::Matrix<Precision, 3, 1> cross =
      w.cross(mu.template block<3, 1>(0, 0));
  if (theta_sq < 1e-8) {
    A = 1.0 - one_6th * theta_sq;
    B = 0.5;
    vec3 = mu.template segment<3>(0) + 0.5 * cross;
  } else {
    Precision C;
    if (theta_sq < 1e-6) {
      C = one_6th*(1.0 - one_20th * theta_sq);
      A = 1.0 - theta_sq * C;
      B = 0.5 - 0.25 * one_6th * theta_sq;
    } else {
      const Precision inv_theta = 1.0/theta;
      A = sin(theta) * inv_theta;
      B = (1 - cos(theta)) * (inv_theta * inv_theta);
      C = (1 - A) * (inv_theta * inv_theta);
    }
    vec3 = mu.template segment<3>(0) + B*cross + C*(w.cross(cross));
  }
  SO3Group<Precision> rot;
  rodrigues_so3_group_exp(w, A, B, rot.matrix_);
  return SE3Group<Precision>(rot, vec3);
}

template <typename Precision>
inline Eigen::Matrix<Precision, 6, 1> SE3Group<Precision>::ln(
    const SE3Group<Precision>& SE3Group) {
  SO3Group<Precision> SO3GroupTemp = SE3Group.get_rotation();
  Eigen::Matrix<Precision, 3, 1> rot = SO3GroupTemp.ln();
  const Precision theta = std::sqrt(rot.transpose()*rot);

  Precision shtot = 0.5;
  if (theta > 0.00001)
    shtot = sin(theta/2)/theta;

  /// now do the rotation
  Eigen::Matrix<Precision, 3, 1> v3Rot = rot*-0.5;
  const SO3Group<Precision> halfrotator = SO3Group<Precision>::exp(v3Rot);
  Eigen::Matrix<Precision, 3, 1> rottrans = halfrotator * SE3Group.get_translation();

  if (theta > 0.001) {
    rottrans -= rot *
        ((SE3Group.get_translation().transpose() * rot) *
         (1-2*shtot) / (rot.transpose()*rot));
  } else {
    rottrans -= rot * ((SE3Group.get_translation().transpose() * rot)/Precision(24.0));
  }
  rottrans /= (2 * shtot);

  Eigen::Matrix<Precision, 6, 1> result;
  result.template segment<3>(0)=rottrans;
  result.template segment<3>(3)=rot;
  return result;
}

// External Operators  =========================================================
/// Write an SE3Group to a stream
/// @relates SE3Group
template <typename Precision>
inline std::ostream& operator << (std::ostream& os, const SE3Group<Precision>& rhs) {
  std::streamsize fw = os.width();
  for (int i = 0; i < 3; i++) {
    os.width(fw);
    os << rhs.get_rotation().get_matrix().row(i);
    os.width(fw);
    os << " " << rhs.get_translation()[i] << '\n';
  }
  return os;
}

/// Reads an SE3Group from a stream
/// @relates SE3Group
template <typename Precision>
inline std::istream& operator >> (std::istream& is, SE3Group<Precision>& rhs) {
  for (int r = 0; r < 3; r++) {
    for (int c = 0; c < 3; c++) {
      is >> const_cast<Eigen::Matrix<Precision, 3, 3>& >(
            rhs.get_rotation().get_matrix())(r, c);
    }
    is >> const_cast<Eigen::Matrix<Precision, 3, 1>& >(rhs.get_translation())[r];
  }
  const_cast<SO3Group<Precision>& >(rhs.get_rotation()).coerce();
  return is;
}

/// Left-multiply by a Matrix
/// @relates SE3Group
template <typename OtherDerived, typename Precision>
inline Eigen::Matrix<typename Eigen::internal::scalar_product_traits<typename OtherDerived::Scalar, Precision>::ReturnType, OtherDerived::RowsAtCompileTime, 4>
  operator * (
    const Eigen::MatrixBase<OtherDerived>& lhs,
    const SE3Group<Precision>& rhs) {
  Eigen::Matrix<typename Eigen::internal::scalar_product_traits<typename OtherDerived::Scalar, Precision>::ReturnType,OtherDerived::RowsAtCompileTime,4> res;
  res.block(0, 0, lhs.rows(), 3) = lhs.block(0, 0, lhs.rows(), 3)*rhs.get_rotation();
  res.block(0, 3, lhs.rows(), 1) = lhs.block(0, 3, lhs.rows(), 1) +
      lhs.block(0, 0, lhs.rows(), 3)*rhs.get_translation();
  return res;
}

/// Multiply a SO2Group with and SE2Group
/// @relates SE3Group
/// @relates SO3Group
template <typename Precision1, typename Precision2>
inline SE3Group<typename Eigen::internal::scalar_product_traits<Precision1, Precision2>::ReturnType >
  operator * (const SO3Group<Precision1> & lhs, const SE3Group<Precision2>& rhs) {
  Eigen::Matrix<typename Eigen::internal::scalar_product_traits<Precision1, Precision2>::ReturnType , 3, 1> v3 = lhs*rhs.get_translation();
  return SE3Group<typename Eigen::internal::scalar_product_traits<Precision1, Precision2>::ReturnType >( lhs*rhs.get_rotation(), v3);
}
}  // namespace slick
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(slick::SE3f)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(slick::SE3d)
