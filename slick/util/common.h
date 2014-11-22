// Copyright 2012, 2013, 2014 The Look3D Authors. All rights reserved.
#ifndef LOOK3D_MATH_UTILITIES_H_
#define LOOK3D_MATH_UTILITIES_H_

#include <algorithm>
#include <Eigen/Eigen>

#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif

// Various utility routines not fitting elsewhere.
namespace look3d {
// For a vector \e v of length \e i, return \f$[v_1, v_2, \cdots, v_{i-1}] / v_i \f$
// @param v \e v
// @ingroup gLinAlg
template<typename OtherDerived>
inline Eigen::Matrix<typename OtherDerived::Scalar, (OtherDerived::RowsAtCompileTime == Eigen::Dynamic ? Eigen::Dynamic : OtherDerived::RowsAtCompileTime - 1) + 0, 1>
  project(const Eigen::MatrixBase<OtherDerived> & v) {
  //  static const int Len = (OtherDerived::RowsAtCompileTime==Eigen::Dynamic?Eigen::Dynamic:OtherDerived::RowsAtCompileTime-1);
  return (v.segment(0, v.size() - 1) / v[v.size() - 1]);
}

// This should probably be done with an operator to prevent an extra new[] for dynamic vectors.
// For a vector \e v of length \e i, return \f$[v_1, v_2, \cdots, v_{i}, 1]\f$
// @param v \e v
// @ingroup gLinAlg
template<typename OtherDerived>
inline Eigen::Matrix<typename OtherDerived::Scalar, (OtherDerived::RowsAtCompileTime == Eigen::Dynamic ? Eigen::Dynamic : OtherDerived::RowsAtCompileTime + 1) + 0, 1>
  unproject(const Eigen::MatrixBase<OtherDerived> & v) {
  Eigen::Matrix<typename OtherDerived::Scalar, (OtherDerived::RowsAtCompileTime == Eigen::Dynamic ? Eigen::Dynamic : OtherDerived::RowsAtCompileTime + 1), 1> result(v.size() + 1);
  //  static const int Len = (OtherDerived::RowsAtCompileTime==Eigen::Dynamic?Eigen::Dynamic:OtherDerived::RowsAtCompileTime);
  result.segment(0, v.size()) = v;
  result[v.size()] = 1;
  return result;
}

// Compute the \f$L_\infty\f$ norm of \e v
// @param v \e v
// @ingroup gLinAlg
template<typename OtherDerived>
inline typename OtherDerived::Scalar inf_norm(
  const Eigen::MatrixBase<OtherDerived>& v) {
  using std::abs;
  using std::max;
  typename OtherDerived::Scalar n = 0;
  n = abs(v[0]);

  for (int i = 1; i < v.size(); i++)
    n = max(n, abs(v[i]));
  return n;
}

template<typename Derived>
inline std::istream & operator >> (std::istream& in, Eigen::MatrixBase<Derived>& m) {
  for (int r = 0; r < m.rows(); ++r) {
    for (int c = 0; c < m.cols(); ++c) {
      in >> m(r, c);
    }
  }
  return in;
}
}  // namespace look3d
#endif  // LOOK3D_MATH_UTILITIES_H_