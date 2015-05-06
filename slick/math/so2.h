// Copyright 2014 The Slick Authors. All rights reserved.
#pragma once
#include <Eigen/Core>
#include <type_traits>
#include <Eigen/Eigen>
#include "slick/datatypes.h"

namespace slick {

template <typename Scalar> class SO2Group;
template <typename Scalar> class SE2Group;

using SO2 = SO2Group<SlickScalar>;
using SO2d = SO2Group<double>;
using SO2f = SO2Group<float>;

/// Class to represent two-dimensional rotation matrix.
/// matrices are members of the Special Orthogonal Lie group SO2Group.
/// This group can be parameterised with one number (the rotation angle).
template <typename Scalar = SlickScalar> class SO2Group {
  using MatrixType = Eigen::Matrix<Scalar, 2, 2>;

public:
  friend class SE2Group<Scalar>;
  // Constructors ==============================================================
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  /// Default constructor. Initialises the matrix to the identity (no rotation)
  SO2Group() { matrix_.setIdentity(); }

  /// Construct from a rotation matrix.
  template <typename Derived> SO2Group(const Eigen::MatrixBase<Derived> &rhs) {
    *this = rhs;
  }

  /// Construct from an angle.
  explicit SO2Group(const Scalar l) { *this = this->exp(l); }

  // this strictly require C++11 enabled
  /// Construct from list of 4 elements in row-major order of 2x2 matrix.
  /// Currently, giving a list of less than 4 elements causes std::cerr.
  SO2Group(std::initializer_list<Scalar> l) {
    // this can be solved with static_assert in C++14
    assert(l.size() == 4);
    Scalar *t = const_cast<Scalar *>(l.begin());
    matrix_(0, 0) = *t++;
    matrix_(0, 1) = *t++;
    matrix_(1, 0) = *t++;
    matrix_(1, 1) = *t;
    coerce();
  }

  // Operators  ================================================================
  /// Right-multiply by another rotation matrix
  template <typename P>
  SO2Group<Scalar> operator*(const SO2Group<P> &rhs) const {
    //    static_assert(std::is_same<P, Scalar>::value, "SO2::operator*
    //    expecting same scalar type");
    SO2Group<Scalar> res;
    res.matrix_ = matrix_ * rhs.matrix_;
    return res;
  }

  /// Self right-multiply by another rotation matrix
  template <typename P> SO2Group &operator*=(const SO2Group<P> &rhs) {
    matrix_ *= rhs.matrix_;
    return *this;
  }

  /// Required for Eigen inheritance
  template <typename OtherDerived>
  SO2Group &operator=(const Eigen::MatrixBase<OtherDerived> &other) {
    matrix_ = other;
    coerce();
    return *this;
  }

  /// Right-multiply by another matrix
  template <typename OtherDerived>
  typename Eigen::ProductReturnType<MatrixType, OtherDerived>::Type
  operator*(const Eigen::MatrixBase<OtherDerived> &other) const {
    return matrix_ * other;
  }

  // SO2Group Specific  ========================================================
  /// returns generator matrix
  static Eigen::Matrix<Scalar, 2, 2> generator() {
    Eigen::Matrix<Scalar, 2, 2> result;
    result.row(0) = Eigen::Matrix<Scalar, 2, 1>(0, -1);
    result.row(1) = Eigen::Matrix<Scalar, 2, 1>(1, 0);
    return result;
  }

  /// Modifies the matrix to make sure it is a valid rotation matrix.
  void coerce() {
    matrix_.row(0).normalize();
    matrix_.row(1) -=
        (matrix_.row(0) * matrix_.row(1).transpose()) * matrix_.row(0);
    matrix_.row(1).normalize();
  }

  /// Exponentiate an angle in the Lie algebra to generate a new SO2Group.
  inline static SO2Group exp(const Scalar &d) {
    SO2Group<Scalar> res;
    res.matrix_(0, 0) = res.matrix_(1, 1) = std::cos(d);
    res.matrix_(1, 0) = std::sin(d);
    res.matrix_(0, 1) = -res.matrix_(1, 0);
    return res;
  }

  /// extracts the rotation angle from the SO2Group
  Scalar ln() const { return atan2(matrix_(1, 0), matrix_(0, 0)); }

  /// Returns the inverse of this matrix (=the transpose, so this is a fast
  /// operation)
  SO2Group<Scalar> inverse() const {
    return SO2Group<Scalar>(matrix_.transpose());
  }

  const MatrixType &get_matrix() const { return matrix_; }

protected:
  MatrixType matrix_;
};

// External Operators  =========================================================
/// Write an SO2Group to a stream
/// @relates SO2Group
template <typename Scalar>
inline std::ostream &operator<<(std::ostream &os, const SO2Group<Scalar> &rhs) {
  return os << rhs.get_matrix();
}

/// Read from a stream to SO2Group
/// @relates SO2Group
template <typename Scalar>
inline std::istream &operator>>(std::istream &is, SO2Group<Scalar> &rhs) {
  Eigen::Matrix<Scalar, 2, 2> m;
  return is >> m;
  rhs = m;
}

/// Left-multiply by a Matrix
/// @relates SO2Group
template <typename OtherDerived, typename Scalar>
inline typename Eigen::ProductReturnType<OtherDerived,
                                         Eigen::Matrix<Scalar, 2, 2>>::Type
operator*(const Eigen::MatrixBase<OtherDerived> &lhs,
          const SO2Group<Scalar> &rhs) {
  return lhs * rhs.get_matrix();
}
} // namespace slick
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(slick::SO2f)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(slick::SO2d)
