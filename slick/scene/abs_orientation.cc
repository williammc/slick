// Copyright 2014 The Slick Authors. All rights reserved.
#include "slick/scene/abs_orientation.h"
#include <scene/abs_orientation.h>
#include <Eigen/Eigenvalues>

using std::vector;
namespace slick {
Eigen::Matrix<DefaultScalarType, 3, 3> QuaternionToMatrix(const Eigen::Matrix<DefaultScalarType, 4, 1> & q) {
  Eigen::Matrix<DefaultScalarType, 3, 3> result;
  const int w = 0, x = 1, y = 2, z = 3;
  result(0, 0) = q[w]*q[w] + q[x]*q[x] - q[y]*q[y] - q[z]*q[z];
  result(0, 1) = 2*(q[x]*q[y] - q[w]*q[z]);
  result(1, 0) = 2*(q[x]*q[y] + q[w]*q[z]);
  result(0, 2) = 2*(q[x]*q[z] + q[w]*q[y]);
  result(2, 0) = 2*(q[x]*q[z] - q[w]*q[y]);
  result(1, 1) = q[w]*q[w] - q[x]*q[x] + q[y]*q[y] - q[z]*q[z];
  result(1, 2) = 2*(q[y]*q[z] - q[w]*q[x]);
  result(2, 1) = 2*(q[y]*q[z] + q[w]*q[x]);
  result(2, 2) = q[w]*q[w] - q[x]*q[x] - q[y]*q[y] + q[z]*q[z];
  return result;
}

SO3 ComputeOrientation(const std::vector<Eigen::Matrix<DefaultScalarType, 3, 1> > & a,
                              const std::vector<Eigen::Matrix<DefaultScalarType, 3, 1> > & b) {
  const size_t N = a.size();
  // compute cross correlations
  const int x = 0, y = 1, z = 2;
  Eigen::Matrix<DefaultScalarType, 3, 3> s = Eigen::Matrix<DefaultScalarType, 3, 3>::Zero();
  for (unsigned int i = 0; i < N; i++) {
    s += a[i] * b[i].transpose();
  }

  // create symmetric M for eigenvalue analysis
  Eigen::Matrix4d M;
  M(0, 0) = s(x, x) + s(y, y) + s(z, z);
  M(1, 0) = M(0, 1) = s(y, z) - s(z, y);
  M(2, 0) = M(0, 2) = s(z, x) - s(x, z);
  M(3, 0) = M(0, 3) = s(x, y) - s(y, x);
  M(1, 1) = s(x, x) - s(y, y) - s(z, z);
  M(2, 1) = M(1, 2) = s(x, y) + s(y, x);
  M(3, 1) = M(1, 3) = s(z, x) + s(x, z);
  M(2, 2) = - s(x, x) + s(y, y) - s(z, z);
  M(3, 2) = M(2, 3) = s(y, z) + s(z, y);
  M(3, 3) = - s(x, x) - s(y, y) + s(z, z);

  // eigenvalue decomposition to find eigenvector to largest eigenvalue
  Eigen::EigenSolver<Eigen::Matrix4d> ev(M);
  Eigen::Matrix<DefaultScalarType, 4, 1> evals = ev.eigenvalues().real();
  int index = 0;
  for (unsigned int i = index+1; i < 4; i++)
    if (evals[i] > evals[index])
      index = i;
  Eigen::Matrix<DefaultScalarType, 4, 1> evec = ev.eigenvectors().col(index).real();
  SO3  result;
  result = QuaternionToMatrix(evec);
  return result;
}

// computes the orientation from (e1,e2,e3) -> (a,(a^b)^a,a^b),
// which means that b the second vector is in the a, b plane
static inline SO3 canonicalOrientation(const Eigen::Matrix<DefaultScalarType, 3, 1> & a,
                                         const Eigen::Matrix<DefaultScalarType, 3, 1> & b) {
  Eigen::Matrix<DefaultScalarType, 3, 1> n = a.cross(b);
  if (n.squaredNorm() < 1e-30)
    return SO3();
  Eigen::Matrix<DefaultScalarType, 3, 3> result;
  result.col(0) = a.normalized();
  result.col(2) = n.normalized();
  result.col(1) = result.col(2).cross(result.col(0));
  return SO3 (result);
}

SO3 ComputeOrientation(const Eigen::Matrix<DefaultScalarType, 3, 1> & a1,
                              const Eigen::Matrix<DefaultScalarType, 3, 1> & b1,
                              const Eigen::Matrix<DefaultScalarType, 3, 1> & a2,
                              const Eigen::Matrix<DefaultScalarType, 3, 1> & b2) {
  SO3  r1 = canonicalOrientation(a1, a2);
  SO3  r2 = canonicalOrientation(b1, b2);
  const SO3  rAB = r2 * r1.inverse();
  r1 = canonicalOrientation(a2, a1);
  r2 = canonicalOrientation(b2, b1);
  const SO3  rBA = r2 * r1.inverse();
  const SO3  diff = rBA * rAB.inverse();
  return SO3 ::exp(diff.ln() * 0.5) * rAB;
}

SE3 ComputeAbsoluteOrientation(const std::vector<Eigen::Matrix<DefaultScalarType, 3, 1> > & a,
                                      const std::vector<Eigen::Matrix<DefaultScalarType, 3, 1> > & b) {
  // std::assert(a.size() <= b.size());
  const size_t N = a.size();
  Eigen::Matrix<DefaultScalarType, 3, 1> ma = Eigen::Matrix<DefaultScalarType, 3, 1>::Zero(), mb = Eigen::Matrix<DefaultScalarType, 3, 1>::Zero();
  // compute centroids
  for (unsigned int i = 0; i < N; i++) {
    ma += a[i];
    mb += b[i];
  }
  ma /= N;
  mb /= N;

  // compute shifted locations
  std::vector<Eigen::Matrix<DefaultScalarType, 3, 1> > ap(N), bp(N);
  for (unsigned int i = 0; i < N; i++) {
    ap[i] = a[i] - ma;
    bp[i] = b[i] - ma;
  }

  // put resulting transformation together
  SO3 rot = ComputeOrientation(ap, bp);
  Eigen::Matrix<DefaultScalarType, 3, 1> trans = mb - rot * ma;
  return SE3(rot, trans);
}

SO3 ComputeMeanOrientation(const std::vector<SO3 > & r) {
  const size_t N = r.size();
  std::vector<SO3 > rt(N);
  SO3  base = r.front();
  SO3  baseInv = base.inverse();
  Eigen::Matrix<DefaultScalarType, 3, 1> center = Eigen::Matrix<DefaultScalarType, 3, 1>::Zero();
  for (unsigned int i = 0; i < N; i++) {
    rt[i] = r[i] * baseInv;
    center += rt[i].ln();
  }
  center /= N;
  SO3 mean(center);
  do {
    center.setZero();
    for (unsigned int i = 0; i < N; i++) {
      SO3  diff = rt[i] * mean.inverse();
      center += diff.ln();
    }
    center /= N;
    mean = SO3::exp(center) * mean;
  } while (center.norm() > 1e-12);

  return mean * base;
}

}  // namespace slick
