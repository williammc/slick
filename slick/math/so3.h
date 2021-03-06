// Copyright 2014 The Slick Authors. All rights reserved.
#pragma once
#include <cassert>
#include <math.h>

#include <Eigen/Core>
#include "slick/datatypes.h"

namespace slick {

template <typename Scalar>
class SO3Group;
template <typename Scalar>
class SE3Group;

using SO3 = SO3Group<SlickScalar>;
using SO3d = SO3Group<double>;
using SO3f = SO3Group<float>;

/// Class to represent a three-dimensional rotation matrix.
/// matrices are members of the Special Orthogonal Lie group SO3Group.
/// This group can be parameterised three numbers
/// (a vector in the space of the Lie Algebra).
/// In this class, the three parameters are the
/// finite rotation vector, i.e. a three-dimensional vector
/// whose direction is the axis of rotation
/// and whose length is the angle of rotation in radians.
/// Exponentiating this vector gives the matrix,
/// and the logarithm of the matrix gives this vector.
template <typename Scalar = SlickScalar>
class SO3Group {
  using MatrixType = Eigen::Matrix<Scalar, 3, 3>;

 public:
  friend class SE3Group<Scalar>;
  // Constructors  =============================================================
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  /// Default constructor. Initialises the matrix to the identity (no rotation)
  SO3Group() { matrix_.setIdentity(); }

  /// Construct from the axis of rotation (and angle given by the magnitude).
  template <typename P, int O>
  explicit SO3Group(const Eigen::Matrix<P, 3, 1, O>& v) {
    *this = this->exp(v);
  }

  /// Construct from a rotation matrix.
  template <typename P, int O>
  SO3Group(const Eigen::Matrix<P, 3, 3, O>& rhs) {
    *this = rhs;
  }

  /// creates an SO3Group as a rotation that takes Vector a into the direction
  /// of Vector b
  /// with the rotation axis along a ^ b. If |a ^ b| == 0, it creates the
  /// identity rotation.
  /// An assertion will fail if Vector a and Vector b are in exactly opposite
  /// directions.
  /// @param a source Vector
  /// @param b target Vector
  template <typename P1, typename P2, int O1, int O2>
  SO3Group(const Eigen::Matrix<P1, 3, 1, O1>& a,
           const Eigen::Matrix<P2, 3, 1, O2>& b) {
    Eigen::Matrix<Scalar, 3, 1> n = a.cross(b);
    if (n.norm() == 0) {
      /// check that the vectors are in the same direction if cross product is
      /// 0.
      /// If not,
      /// this means that the rotation is 180 degrees,
      /// which leads to an ambiguity in the rotation axis.
      assert((a.transpose() * b) >= 0);
      matrix_.setIdentity();
      return;
    }
    n.normalize();
    Eigen::Matrix<Scalar, 3, 3> R1;
    R1.col(0) = a.normalized();
    R1.col(1) = n;
    R1.col(2) = n.cross(R1.col(0));
    matrix_.col(0) = b.normalized();
    matrix_.col(1) = n;
    matrix_.col(2) = n.cross(matrix_.col(0));
    matrix_ *= R1.transpose();
  }

  // Operators  ================================================================
  /// Assignment operator from a general matrix. This also calls coerce()
  /// to make sure that the matrix is a valid rotation matrix.
  template <typename Deri33>
  SO3Group& operator=(const Eigen::MatrixBase<Deri33>& rhs) {
    matrix_ = rhs;
    coerce();
    return *this;
  }

  /// Right-multiply by another rotation matrix
  SO3Group& operator*=(const SO3Group& rhs) {
    matrix_ *= rhs.matrix_;
    return *this;
  }

  /// Right-multiply by another rotation matrix
  template <typename P>
  SO3Group<
      typename Eigen::internal::scalar_product_traits<Scalar, P>::ReturnType>
  operator*(const SO3Group<P>& rhs) const {
    SO3Group<typename Eigen::internal::scalar_product_traits<
        Scalar, P>::ReturnType> SO3GroupTemp;
    SO3GroupTemp = matrix_ * rhs.get_matrix();
    return SO3GroupTemp;
  }

  /// Returns the i-th generator times pos
  template <typename Deri31>
  static Eigen::Matrix<typename Deri31::Scalar, 3, 1> generator_field(
      int i, const Eigen::MatrixBase<Deri31>& pos) {
    Eigen::Matrix<typename Deri31::Scalar, 3, 1> result;
    result[i] = 0;
    result[(i + 1) % 3] = -pos[(i + 2) % 3];
    result[(i + 2) % 3] = pos[(i + 1) % 3];
    return result;
  }

  /// Transfer a vector in the Lie Algebra from one
  /// co-ordinate frame to another such that for a matrix
  /// \f$ M \f$, the adjoint \f$Adj()\f$ obeys
  /// \f$ e^{\text{Adj}(v)} = Me^{v}M^{-1} \f$
  template <int S, int O>
  Eigen::Matrix<Scalar, 3, 1> adjoint(
      const Eigen::Matrix<Scalar, S, 1, O>& vect) const {
    assert(S == vect.rows());
    return *this * vect;
  }

  template <typename PA, typename PB>
  inline SO3Group(const SO3Group<PA>& a, const SO3Group<PB>& b) {
    *this = a * b;
  }

  /// Right-multiply by another matrix
  template <typename OtherDerived>
  typename Eigen::ProductReturnType<MatrixType, OtherDerived>::Type operator*(
      const Eigen::MatrixBase<OtherDerived>& other) const {
    return matrix_ * other;
  }

  // SO3Group Specific  ========================================================
  /// Modifies the matrix to make sure it is a valid rotation matrix.
  void coerce() {
    matrix_.row(0).normalize();
    matrix_.row(1) -=
        (matrix_.row(0) * (matrix_.row(1).transpose())) * matrix_.row(0);
    matrix_.row(1).normalize();
    matrix_.row(2) -=
        (matrix_.row(0) * (matrix_.row(2).transpose())) * matrix_.row(0);
    matrix_.row(2) -=
        (matrix_.row(1) * (matrix_.row(2).transpose())) * matrix_.row(1);
    matrix_.row(2).normalize();
    /// check for positive determinant <=> correct rotation matrix should have
    /// determinant == 1
    if ((matrix_.row(0).cross(matrix_.row(1))).dot(matrix_.row(2)) <=
        Scalar(0)) {
      throw("Invalid Rotation Matrix");
    }
  }

  /// Exponentiate a vector in the Lie algebra to generate a new SO3Group.
  /// See the Detailed Description for details of this vector.
  template <typename Deri31>
  static SO3Group exp(const Eigen::MatrixBase<Deri31>& vect);

  /// Take the logarithm of the matrix, generating the corresponding vector in
  /// the Lie Algebra.
  /// See the Detailed Description for details of this vector.
  inline Eigen::Matrix<Scalar, 3, 1> ln() const;

  /// Returns the inverse of this matrix (=the transpose, so this is a fast
  /// operation)
  SO3Group inverse() const {
    SO3Group<Scalar> SO3GroupInv;
    SO3GroupInv = matrix_.transpose();
    return SO3GroupInv;
  }
  const MatrixType& get_matrix() const { return matrix_; }

  template <typename P>
  SO3Group<P> cast() {
    SO3Group<P> SO3Group;
    SO3Group.matrix_ = matrix_.template cast<P>();
    return SO3Group;
  }

 protected:
  MatrixType matrix_;
};  // class SO3Group

/// Compute a rotation exponential using the Rodrigues Formula.
/// The rotation axis is given by \f$\vec{w}\f$, and the rotation angle must
/// be computed using \f$ \theta = |\vec{w}|\f$. This is provided as a separate
/// function primarily to allow fast and rough matrix exponentials using fast
/// and rough approximations to \e A and \e B.
///
/// @param w Vector about which to rotate.
/// @param A \f$\frac{\sin \theta}{\theta}\f$
/// @param B \f$\frac{1 - \cos \theta}{\theta^2}\f$
/// @param R Matrix to hold the return value.
/// @relates SO3Group
template <typename Scalar, typename Deri31, typename Deri33>
inline void rodrigues_so3_group_exp(const Eigen::MatrixBase<Deri31>& w,
                                    const Scalar A, const Scalar B,
                                    Eigen::MatrixBase<Deri33>& R) {
  const Scalar wx2 = (Scalar)w(0, 0) * w(0, 0);
  const Scalar wy2 = (Scalar)w(1, 0) * w(1, 0);
  const Scalar wz2 = (Scalar)w(2, 0) * w(2, 0);

  R(0, 0) = 1.0 - B * (wy2 + wz2);
  R(1, 1) = 1.0 - B * (wx2 + wz2);
  R(2, 2) = 1.0 - B * (wx2 + wy2);
  {
    const Scalar a = A * w(2, 0);
    const Scalar b = B * (w(0, 0) * w(1, 0));
    R(0, 1) = b - a;
    R(1, 0) = b + a;
  }
  {
    const Scalar a = A * w(1, 0);
    const Scalar b = B * (w(0, 0) * w(2, 0));
    R(0, 2) = b + a;
    R(2, 0) = b - a;
  }
  {
    const Scalar a = A * w(0, 0);
    const Scalar b = B * (w(1, 0) * w(2, 0));
    R(1, 2) = b - a;
    R(2, 1) = b + a;
  }
}

/// Perform the exponential of the matrix \f$ \sum_i w_iG_i\f$
/// @param w Weightings of the generator matrices.
template <typename Scalar>
template <typename Deri31>
inline SO3Group<Scalar> SO3Group<Scalar>::exp(
    const Eigen::MatrixBase<Deri31>& w) {
  using std::sqrt;
  using std::sin;
  using std::cos;

  static const Scalar one_6th = Scalar(1.0) / Scalar(6.0);
  static const Scalar one_20th = Scalar(1.0) / Scalar(20.0);

  const Scalar theta_sq = w.col(0).squaredNorm();
  const Scalar theta = sqrt(theta_sq);
  Scalar A, B;
  /// Use a Taylor series expansion near zero. This is required for
  /// accuracy, since sin t / t and (1-cos t)/t^2 are both 0/0.
  if (theta_sq < Scalar(1e-8)) {
    A = Scalar(1.0) - one_6th * theta_sq;
    B = Scalar(0.5);
  } else {
    if (theta_sq < Scalar(1e-6)) {
      B = Scalar(0.5) - Scalar(0.25) * one_6th * theta_sq;
      A = Scalar(1.0) -
          theta_sq * one_6th * (Scalar(1.0) - one_20th * theta_sq);
    } else {
      const Scalar inv_theta = Scalar(1.0) / theta;
      A = sin(theta) * inv_theta;
      B = (Scalar(1) - cos(theta)) * (inv_theta * inv_theta);
    }
  }
  MatrixType m;
  m.setIdentity();
  rodrigues_so3_group_exp(w, A, B, m);

  SO3Group<Scalar> res;
  res = m;
  return res;
}

template <typename Scalar>
inline Eigen::Matrix<Scalar, 3, 1> SO3Group<Scalar>::ln() const {
  Eigen::Matrix<Scalar, 3, 1> result;

  const Scalar cos_angle =
      (matrix_(0, 0) + matrix_(1, 1) + matrix_(2, 2) - Scalar(1.0)) *
      Scalar(0.5);
  result[0] = (matrix_(2, 1) - matrix_(1, 2)) / Scalar(2);
  result[1] = (matrix_(0, 2) - matrix_(2, 0)) / Scalar(2);
  result[2] = (matrix_(1, 0) - matrix_(0, 1)) / Scalar(2);
  const Scalar t = result.squaredNorm();
  Scalar sin_angle_abs = (t > 0.0) ? sqrt(t) : Scalar(0.0);
  if (cos_angle > Scalar(M_SQRT1_2)) {  /// [0 - Pi/4] use asin
    if (sin_angle_abs > Scalar(0)) {
      result *= asin(sin_angle_abs) / sin_angle_abs;
    }
  } else if (cos_angle > -Scalar(M_SQRT1_2)) {  /// [Pi/4 - 3Pi/4] use acos,
                                                   /// but antisymmetric part
    const Scalar angle = acos(cos_angle);
    result *= angle / sin_angle_abs;
  } else {  /// rest use symmetric part
    /// antisymmetric part vanishes, but still large rotation, need information
    /// from symmetric part
    const Scalar angle = Scalar(M_PI) - asin(sin_angle_abs);
    const Scalar d0 = matrix_(0, 0) - cos_angle,
                    d1 = matrix_(1, 1) - cos_angle,
                    d2 = matrix_(2, 2) - cos_angle;
    Eigen::Matrix<Scalar, 3, 1> r2;
    if (d0 * d0 > d1 * d1 &&
        d0 * d0 > d2 * d2) {  /// first is largest, fill with first column
      r2[0] = d0;
      r2[1] = (matrix_(1, 0) + matrix_(0, 1)) / Scalar(2);
      r2[2] = (matrix_(0, 2) + matrix_(2, 0)) / Scalar(2);
    } else if (d1 * d1 >
               d2 * d2) {  /// second is largest, fill with second column
      r2[0] = (matrix_(1, 0) + matrix_(0, 1)) / Scalar(2);
      r2[1] = d1;
      r2[2] = (matrix_(2, 1) + matrix_(1, 2)) / Scalar(2);
    } else {  /// third is largest, fill with third column
      r2[0] = (matrix_(0, 2) + matrix_(2, 0)) / Scalar(2);
      r2[1] = (matrix_(2, 1) + matrix_(1, 2)) / Scalar(2);
      r2[2] = d2;
    }
    /// flip, if we point in the wrong direction!
    if (r2.dot(result) < Scalar(0)) r2 *= Scalar(-1);
    result = angle * r2.normalized();
  }
  return result;
}

/// External Operators  ========================================================
/// Write an SO3Group to a stream
/// @relates SO3Group
template <typename Scalar>
inline std::ostream& operator<<(std::ostream& os,
                                const SO3Group<Scalar>& rhs) {
  return os << rhs.get_matrix();
}

/// Read from a stream to SO3Group
/// @relates SO3Group
template <typename Scalar>
inline std::istream& operator>>(std::istream& is, SO3Group<Scalar>& rhs) {
  Eigen::Matrix<Scalar, 3, 3> m;
  return is >> m;
  rhs = m;
}

/// Left-multiply by a Matrix
/// @relates SO3Group
template <typename OtherDerived, typename Scalar>
inline typename Eigen::ProductReturnType<OtherDerived,
                                         Eigen::Matrix<Scalar, 3, 3> >::Type
operator*(const Eigen::MatrixBase<OtherDerived>& lhs,
          const SO3Group<Scalar>& rhs) {
  return lhs * rhs.get_matrix();
}
}  // end namespace slick
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(slick::SO3f)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(slick::SO3d)
