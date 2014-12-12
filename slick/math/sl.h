// Copyright 2014 The Slick Authors. All rights reserved.
#pragma once

#include <cassert>

#include <Eigen/Eigen>
#include <Eigen/MatrixFunctions>

#include "slick/datatypes.h"

namespace slick {

template <typename P, int N> class SL;

// represents an element from the group SL(n), NxN matrices M with det(M) = 1.
// This can be used to estimate homographies on n-1 dimentional spaces.
// The implementation uses the matrix exponential function @ref exp for
// exponentiation from an element in the Lie algebra and LU to compute inverse.
//
// The Lie algebra are the NxN matrices M with trace(M) = 0.
// The N*N-1 generators used to represent this vector space are the following:
// - diag(...,1,-1,...), n-1 along the diagonal
// - symmetric generators for every pair of off-diagonal elements
// - anti-symmetric generators for every pair of off-diagonal elements
// This choice represents the fact that SL(n) can be interpreted as the product
// of all symmetric matrices with det() = 1 times SO(n).
template <typename Precision = SlickScalar, int N = Eigen::Dynamic>
class SL {
  typedef Eigen::Matrix<Precision, N, N> MatrixType;

 public:
  static const int size = N;  // size of the matrices represented by SL<N>
  static const int dim = N*N - 1;  // dimension of the vector space SL<N>

  // ======================== Constructors =========================///
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  // default constructor, creates identity element
  SL() { matrix_.setIdentity(); }

  // exp constructor, creates element through exponentiation of Lie algebra
  template <int S, typename P, int B>
  SL(const Eigen::Matrix<P, S, 1, B> & v ) { *this = exp(v); }

  // copy constructor from a matrix, coerces matrix to be of determinant = 1
  template <int R, int C, typename P, int A>
  SL(const Eigen::Matrix<P, R, C, A>& M) { *this = M; coerce(); }

  // ========================= Operators =================================///
  // multiplies to SLs together by multiplying the underlying matrices
  template<typename P>
  SL<typename Eigen::internal::scalar_product_traits<Precision, P>::ReturnType, N>
    operator * (const SL<P, N> & rhs) const {
    SL<typename Eigen::internal::scalar_product_traits<Precision, P>::ReturnType, N> sl_temp;
    sl_temp = matrix_*rhs.getMatrix();
    return sl_temp;
  }
  // right multiplies this SL with another one
  SL operator *= (const SL & rhs) {
    *this = *this*rhs;
    return *this;
  }

  // Required for Eigen inheritance
  template<typename OtherDerived>
  SL& operator = (const Eigen::MatrixBase <OtherDerived>& other) {
    matrix_ = other;
    coerce();
    return *this;
  }

  // Right-multiply by another matrix
  template<typename OtherDerived>
  typename Eigen::ProductReturnType< MatrixType, OtherDerived >::Type
    operator * (const Eigen::MatrixBase <OtherDerived>& other) const {
    return matrix_*other;
  }

  // =========================== SLN Specific ==============================///
  // returns the inverse using LU
  SL inverse() const {
    MatrixType m = matrix_.inverse();
    SL<Precision, N> slTemp = m;
    return slTemp;
  }

  // exponentiates a vector in the Lie algebra to compute the corresponding element
  // @arg v a vector of dimension SL::dim
  template <int S, typename P, int B>
  static inline SL exp(const Eigen::Matrix<P, S, 1, B> &);

  // returns one generator of the group. see SL for a detailed description of
  // the generators used.
  // @arg i number of the generator between 0 and SL::dim -1 inclusive
  static inline Eigen::Matrix<Precision, N, N> generator(int i);

  inline const MatrixType& getMatrix() const { return matrix_;}

 private:
  SL(const SL & a, const SL & b) {*this = a*b;}

  void coerce() {
    using std::abs;
    Precision det = matrix_.determinant();
    assert(abs(det) > 0);
        using std::pow;
    matrix_ /= pow(det, 1.0/N);
  }

  // these constants indicate which parts of the parameter vector
  // map to which generators
  static const int COUNT_DIAG = N - 1;
  static const int COUNT_SYMM = (dim - COUNT_DIAG)/2;
  static const int COUNT_ASYMM = COUNT_SYMM;
  static const int DIAG_LIMIT = COUNT_DIAG;
  static const int SYMM_LIMIT = COUNT_SYMM + DIAG_LIMIT;
  MatrixType matrix_;
};  // class SL

template <typename Precision, int N>
template <int S, typename P, int B>
inline SL<Precision, N> SL<Precision, N>::exp(
    const Eigen::Matrix<P, S, 1, B> & v) {
  assert(S == v.size());
  Eigen::Matrix<Precision, N, N> t= Eigen::Matrix<Precision, N, N>::Zero();
  for (int i = 0; i < dim; ++i)
    t += generator(i) * v[i];
  SL<Precision, N> result;
  Eigen::Matrix<Precision, N, N> mN;
  // from MatrixBase::exp() => (struc) MatrixExponentialReturnValue
  // (t.exp()).evalTo(mN);
  // from unsupported Eigen class
  Eigen::MatrixExponential<Eigen::Matrix<Precision, N, N> > mExp(t);
  mExp.compute(mN);
  result = mN;
  return result;
}

template <typename Precision, int N>
inline Eigen::Matrix<Precision, N, N> SL<Precision, N>::generator(int i) {
  assert(i > -1 && i < dim);
  Eigen::Matrix<Precision, N, N> result = Eigen::Matrix<Precision, N, N>::Zero();
  if (i < DIAG_LIMIT) {  // first ones are the diagonal ones
    result(i, i) = 1;
    result(i+1, i+1) = -1;
  } else if (i < SYMM_LIMIT) {  // then the symmetric ones
    int row = 0, col = i - DIAG_LIMIT + 1;
    while (col > (N - row - 1)) {
      col -= (N - row - 1);
      ++row;
    }
    col += row;
    result(row, col) = result(col, row) = 1;
  } else {  // finally the antisymmetric ones
    int row = 0, col = i - SYMM_LIMIT + 1;
    while (col > N - row - 1) {
      col -= N - row - 1;
      ++row;
    }
    col += row;
    result(row, col) = -1;
    result(col, row) = 1;
  }
  return result;
}
// ========================= External Operators ========================///
// Write an SL to a stream
// @relates SL
template <typename Precision, int N>
inline std::ostream& operator << (
    std::ostream& os, const SL<Precision, N> & rhs) {
  return os << rhs.getMatrix();
}

// Read from a stream to SL
// @relates SL
template <typename Precision, int N>
inline std::istream& operator >> (std::istream& is, SL<Precision, N>& rhs) {
  Eigen::Matrix<Precision, N, N> m;
  return is >> m;
  rhs = m;
}

// Left-multiply by a Matrix
// @relates SL
template <typename OtherDerived, typename Precision, int N>
inline typename Eigen::ProductReturnType< OtherDerived, Eigen::Matrix<Precision, N, N> >::Type
  operator * (const Eigen::MatrixBase<OtherDerived>& lhs, const SL<Precision, N>& rhs) {
  return lhs * rhs.getMatrix();
}
}  // namespace slick
