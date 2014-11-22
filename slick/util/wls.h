// Copyright 2014 The Slick Authors. All rights reserved.
#pragma once
#include <cmath>
#include <Eigen/Dense>
#include "slick/datatypes.h"
#include "slick/slick_api.h"

namespace slick {

// Performs Gauss-Newton weighted least squares computation.
// @param Size The number of dimensions in the system
// @param Scalar The numerical Scalar used (DefaultScalarType, float etc)
// @param Decomposition The class used to invert the inverse Covariance matrix
// (must have one integer size and one typename Scalar template arguments)
// this is Cholesky by default, but could also be SQSVD
// @ingroup gEquations
template <class Scalar = DefaultScalarType, int Size = Eigen::Dynamic,
    class Decomposition = Eigen::LDLT<Eigen::Matrix<Scalar, Size, Size> > >
class SLICK_API WLS {
public:
  // Default constructor or construct with the number of dimensions for the Dynamic case
  explicit WLS(int size = Size) {
    A = Eigen::Matrix<Scalar, Size, Size>(size, size);
    b = Eigen::Matrix<Scalar, Size, 1>(size);
    decomposition = Decomposition(size);
    solution = Eigen::Matrix<Scalar, Size, 1>(size);
    clear();
  }

  // Clear all the measurements and apply a constant regularisation term.
  void clear() {
    A.setZero(Size, Size);
    b.setZero(Size, 1);
  }

  // Applies a constant regularisation term.
  // Equates to a prior that says all the parameters are zero with \f$\sigma^2 = \frac{1}{\text{val}}\f$.
  // @param val The strength of the prior
  void add_prior(Scalar val) {
    for (int i = 0; i < A.rows(); i++) {
      A(i, i) += val;
    }
  }

  // Applies a regularisation term with a different strength for each parameter value.
  // Equates to a prior that says all the parameters are zero with \f$\sigma_i^2 = \frac{1}{\text{v}_i}\f$.
  // @param v The vector of priors
  template<int B2>
  void add_prior(const Eigen::Matrix<Scalar, Size, 1, B2>& v) {
    assert(A.rows() == v.size());
    for (int i = 0; i < A.num_rows(); i++)
      A(i, i) += v[i];
  }

  // Applies a whole-matrix regularisation term.
  // This is the same as adding the \f$m\f$ to the inverse covariance matrix.
  // @param m The inverse covariance matrix to add
  template<typename DeriSizeSize>
  void add_prior(const Eigen::MatrixBase<DeriSizeSize>& m) {
    A.noalias() += m;
  }

  // Add a single measurement
  // @param m The value of the measurement
  // @param J The Jacobian for the measurement \f$\frac{\partial\text{m}}{\partial\text{param}_i}\f$
  // @param weight The inverse variance of the measurement (default = 1)
  template<typename DeriSize1>
  void add_mJ(Scalar m, const Eigen::MatrixBase<DeriSize1>& J, Scalar weight = 1) {
    // Upper right triangle only, for speed
    for (int r = 0; r < A.rows(); r++) {
      DefaultScalarType Jw = weight * J(r, 0);
      b[r] += m * Jw;
      for (int c = r; c < A.rows(); c++)
        A(r, c) += Jw * J(c, 0);
    }
  }

  // Add multiple measurements at once (much more efficiently)
  // @param m The measurements to add
  // @param J The Jacobian matrix \f$\frac{\partial\text{m}_i}{\partial\text{param}_j}\f$
  // @param invcov The inverse covariance of the measurement values
  template<typename DeriN1, typename DeriSizeN, typename DeriNN>
  void add_mJ(const Eigen::MatrixBase<DeriN1>& m,
             const Eigen::MatrixBase<DeriSizeN>& J,
             const Eigen::MatrixBase<DeriNN>& invcov) {
    const Eigen::Matrix<Scalar, Size, DeriNN::RowsAtCompileTime> temp =  J * invcov;
    A.noalias() += temp * J.transpose();
    b.noalias() += temp * m;
  }

  // Add multiple measurements at once (much more efficiently)
  // @param m The measurements to add
  // @param J The Jacobian matrix \f$\frac{\partial\text{m}_i}{\partial\text{param}_j}\f$
  // @param invcov The inverse covariance of the measurement values
  template<typename DeriN1, typename DeriNSize, typename DeriNN>
  void add_mJ_rows(const Eigen::MatrixBase<DeriN1>& m,
             const Eigen::MatrixBase<DeriNSize>& J,
             const Eigen::MatrixBase<DeriNN>& invcov) {
    const Eigen::Matrix<Scalar, Size, DeriNN::RowsAtCompileTime> temp =
        J.transpose() * invcov;
    A.noalias() += temp * J;
    b.noalias() += temp * m;
  }

  // Add a single measurement at once with a sparse Jacobian (much, much more efficiently)
  // @param m The measurements to add
  // @param J1 The first block of the Jacobian matrix \f$\frac{\partial\text{m}_i}{\partial\text{param}_j}\f$
  // @param index1 starting index for the first block
  // @param invcov The inverse covariance of the measurement values
  template<typename DeriN1>
  void add_sparse_mJ(const Scalar m,
                     const Eigen::MatrixBase<DeriN1>& J1,
                     const int index1,
                     const Scalar weight) {
    // Upper right triangle only, for speed
    for (int r = 0; r < J1.size(); r++) {
      DefaultScalarType Jw = weight * J1[r];
      b[r+index1] += m * Jw;
      for (int c = r; c < J1.size(); c++)
        A(r+index1, c+index1).noalias() += Jw * J1[c];
    }
  }

  // Add multiple measurements at once with a sparse Jacobian (much, much more efficiently)
  // @param m The measurements to add
  // @param J1 The first block of the Jacobian matrix \f$\frac{\partial\text{m}_i}{\partial\text{param}_j}\f$
  // @param index1 starting index for the first block
  // @param invcov The inverse covariance of the measurement values
  template<typename DeriN1, typename DeriNS1, typename DeriNN>
  void add_sparse_mJ_rows(const Eigen::MatrixBase<DeriN1>& m,
             const Eigen::MatrixBase<DeriNS1>& J1, const int index1,
             const Eigen::MatrixBase<DeriNN>& invcov) {
    const Eigen::Matrix<Scalar, DeriNS1::ColsAtCompileTime, DeriNS1::RowsAtCompileTime> temp1 = J1.transpose() * invcov;
    const int size1 = J1.cols();
    A.block(index1, index1, size1, size1).noalias() += temp1 * J1;
    b.block(index1, size1).noalias() += temp1 * m;
  }

  // Add multiple measurements at once with a sparse Jacobian (much, much more efficiently)
  // @param m The measurements to add
  // @param J1 The first block of the Jacobian matrix \f$\frac{\partial\text{m}_i}{\partial\text{param}_j}\f$
  // @param index1 starting index for the first block
  // @param J2 The second block of the Jacobian matrix \f$\frac{\partial\text{m}_i}{\partial\text{param}_j}\f$
  // @param index2 starting index for the second block
  // @param invcov The inverse covariance of the measurement values
  template<typename DeriN1, typename DeriNS1, typename DeriNS2, typename DeriNN>
  void add_sparse_mJ_rows(const Eigen::MatrixBase<DeriN1>& m,
             const Eigen::MatrixBase<DeriNS1>& J1, const int index1,
             const Eigen::MatrixBase<DeriNS2>& J2, const int index2,
             const Eigen::MatrixBase<DeriNN>& invcov) {
    const Eigen::Matrix<Scalar, DeriNS1::ColsAtCompileTime, DeriNS1::RowsAtCompiletime> temp1 = J1.transpose() * invcov;
    const Eigen::Matrix<Scalar, DeriNS2::ColsAtCompileTime, DeriNS2::RowsAtCompiletime> temp2 = J2.transpose() * invcov;
    const Eigen::Matrix<Scalar, DeriNS1::ColsAtCompileTime, DeriNS2::ColsAtCompiletime> mixed = temp1 * J2;
    const int size1 = J1.cols();
    const int size2 = J2.cols();
    A.block(index1, index1, size1, size1).noalias() += temp1 * J1;
    A.block(index2, index2, size2, size2).noalias() += temp2 * J2;
    A.block(index1, index2, size1, size2).noalias() += mixed;
    A.block(index2, index1, size2, size1).noalias() += mixed.transpose();
    b.block(index1, size1).noalias() += temp1 * m;
    b.block(index2, size2).noalias() += temp2 * m;
  }

  // Process all the measurements and compute the weighted least squares set of parameter values
  // stores the result internally which can then be accessed by calling get_mu()
  void compute() {
    // Copy the upper right triangle to the empty lower-left.
    A.template triangularView<Eigen::UnitLower>() =
        A.transpose();
    decomposition.compute(A);
    solution = decomposition.solve(b);
  }

  // Combine measurements from two WLS systems
  // @param meas The measurements to combine with
  void operator += (const WLS& meas) {
    b.noalias() += meas.b;
    A.noalias() += meas.A;
  }

  // Returns the inverse covariance matrix
  Eigen::Matrix<Scalar, Size, Size>& get_C_inv() {
    return A;
  }
  // Returns the inverse covariance matrix
  const Eigen::Matrix<Scalar, Size, Size>& get_C_inv() const {
    return A;
  }
  // Returns the update. With no prior, this is the result of \f$J^\dagger e\f$.
  Eigen::Matrix<Scalar, Size, 1>& get_mu() {
    return solution;
  }
  // Returns the update. With no prior, this is the result of \f$J^\dagger e\f$.
  const Eigen::Matrix<Scalar, Size, 1>& get_mu() const {
    return solution;
  }
  // Returns the  vector \f$J^{\mathsf T} e\f$
  Eigen::Matrix<Scalar, Size, 1>& get_vector() { return b; }
  // Returns the  vector \f$J^{\mathsf T} e\f$
  const Eigen::Matrix<Scalar, Size, 1>& get_vector() const {
    return b;
  }
  // Return the decomposition object used to compute \f$(J^{\mathsf T}  J + P)^{-1}\f$
  Decomposition& get_decomposition() { return decomposition; }
  // Return the decomposition object used to compute \f$(J^{\mathsf T}  J + P)^{-1}\f$
  const Decomposition& get_decomposition() const { return decomposition; }

  Eigen::Matrix<Scalar, Size, Size> A;
  Eigen::Matrix<Scalar, Size, 1> b;
  Decomposition decomposition;
  Eigen::Matrix<Scalar, Size, 1> solution;
};  // end class WLS
}  // namespace slick
