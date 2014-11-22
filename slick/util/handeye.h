// Copyright 2012, 2013, 2014 The Look3D Authors. All rights reserved.
#ifndef LOOK3D_GEOMETRY_HANDEYE_H_
#define LOOK3D_GEOMETRY_HANDEYE_H_
#include <utility>
#include <vector>
#include <Eigen/Cholesky>

#include "geometry/abs_orientation.h"
#include "math/so3.h"

namespace look3d {

static inline Eigen::Matrix<DefaultScalarType, 3, 1> GetRotationVector(const SO3  & r) {
  const Eigen::Matrix<DefaultScalarType, 3, 3> & my_matrix = r.get_matrix();
  Eigen::Matrix<DefaultScalarType, 3, 1> result;
  result[0] = (my_matrix(2, 1)-my_matrix(1, 2));
  result[1] = (my_matrix(0, 2)-my_matrix(2, 0));
  result[2] = (my_matrix(1, 0)-my_matrix(0, 1));
  result = result.normalized();
  return result;
}

static inline Eigen::Matrix<DefaultScalarType, 3, 1> GetRotationVector(const SE3  & t) {
  return GetRotationVector(t.get_rotation());
}

template<class T>
static inline SO3  SolveXABX(const std::vector<T> & A,
                               const std::vector<T> & B) {
  std::vector<Eigen::Matrix<DefaultScalarType, 3, 1> > va(A.size()), vb(A.size());
  for (unsigned int i = 0; i < A.size(); ++i) {
    va[i] = GetRotationVector(A[i]);
    vb[i] = GetRotationVector(B[i]);
  }
  return ComputeOrientation(va, vb);
}

static inline Eigen::Matrix<DefaultScalarType, 3, 3> EyeMinus(const Eigen::Matrix<DefaultScalarType, 3, 3> & m) {
  Eigen::Matrix<DefaultScalarType, 3, 3> result = Eigen::Matrix<DefaultScalarType, 3, 3>::Identity();
  result -= m;
  return result;
}

SE3  ComputeHandEyeSingle(const std::vector<SE3 > & AB,
                            const std::vector<SE3 > & CD) {
  std::vector<SE3 > A(AB.size()-1), B(AB.size()-1);
  for (unsigned int i = 0; i < AB.size() - 1; ++i) {
    A[i] = CD[i] * CD[i+1].inverse();
    B[i] = AB[i].inverse() * AB[i+1];
  }
  SO3 R = SolveXABX(A, B);
  Eigen::Matrix<DefaultScalarType, 3, 3> JTJ = Eigen::Matrix<DefaultScalarType, 3, 3>::Zero();
  Eigen::Matrix<DefaultScalarType, 3, 1> JTE = Eigen::Matrix<DefaultScalarType, 3, 1>::Zero();
  for (unsigned int i = 0; i < A.size(); ++i) {
    Eigen::Matrix<DefaultScalarType, 3, 3> m =  EyeMinus(B[i].get_rotation().get_matrix());
    JTJ += m.transpose() * m;
    JTE += m.transpose() * (B[i].get_translation() - R * A[i].get_translation());
  }
  Eigen::LLT<Eigen::Matrix<DefaultScalarType, 3, 3> > chol(JTJ);
  Eigen::Matrix<DefaultScalarType, 3, 1> t = chol.solve(JTE);
  return SE3 (R, t);
}

std::pair<SE3 , SE3 > ComputeHandEye(const std::vector<SE3 > & AB,
                                         const std::vector<SE3 > & CD) {
  assert(AB.size() == CD.size() && AB.size() > 2);
  return std::make_pair(ComputeHandEyeSingle(AB, CD),
                        ComputeHandEyeSingle(CD, AB));
}

SO3  ComputeHandEyeSingle(const std::vector<SO3 > &AB,
                            const std::vector<SO3 > & CD) {
  std::vector<SO3 > A(AB.size()-1), B(AB.size()-1);
  for (unsigned int i = 0; i < AB.size() - 1; ++i) {
    A[i] = CD[i] * CD[i+1].inverse();
    B[i] = AB[i].inverse() * AB[i+1];
  }
  return SolveXABX(A, B);
}

std::pair<SO3 , SO3 > ComputeHandEye(const std::vector<SO3 > &AB,
                                         const std::vector<SO3 > & CD) {
  assert(AB.size() == CD.size() && AB.size() > 2);
  return std::make_pair(ComputeHandEyeSingle(AB, CD),
                        ComputeHandEyeSingle(CD, AB));
}
}  // namespace look3d
#endif  // LOOK3D_GEOMETRY_HANDEYE_H_
