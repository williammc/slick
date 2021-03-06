// Copyright 2014 The Slick Authors. All rights reserved.
#include <iostream>
#include <utility>
#include <vector>
#include <Eigen/Core>

#include "slick/math/se3.h"
#include "slick/util/mestimator.h"
#include "slick/util/common.h"

using ObservationPairFloat =
    std::pair<Eigen::Matrix<float, 4, 1>, Eigen::Matrix<float, 2, 1>>;
using ObservationPairDouble =
    std::pair<Eigen::Matrix<double, 4, 1>, Eigen::Matrix<double, 2, 1>>;

EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(ObservationPairFloat)
// EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(ObservationPairDouble)

namespace slick {

// A function to evaluate x^4 + Bx^3 + Cx^2 + Dx + E
inline SlickScalar eval_quartic(SlickScalar B, SlickScalar C, SlickScalar D,
                                SlickScalar E, SlickScalar x) {
  return E + x * (D + x * (C + x * (B + x)));
}

// A function that performs one iteration of Newton's method
// on the quartic x^4 + Bx^3 + Cx^2 + Dx + E
template <typename Scalar>
inline Scalar eval_newton_quartic(Scalar B, Scalar C, Scalar D, Scalar E,
                                  Scalar x) {
  Scalar fx = E + x * (D + x * (C + x * (B + x)));
  Scalar dx = D + x * (2 * C + x * (3 * B + x * 4));
  return x - fx / dx;
}

template <typename Scalar> inline static Scalar square(Scalar x) {
  return x * x;
}

#if defined(_WIN32)
inline SlickScalar cbrt(SlickScalar x) {
  if (fabs(x) < DBL_EPSILON)
    return 0.0;
  if (x > 0.0)
    return pow(x, 1.0 / 3.0);
  return -pow(-x, 1.0 / 3.0);
}
#endif

template <typename Scalar>
inline int DepressedCubicRealRoots(Scalar P, Scalar Q, Scalar r[]) {
  static const Scalar third = Scalar(1.0 / 3.0);
  Scalar third_P = third * P;
  Scalar disc = 4 * (third_P * third_P * third_P) + Q * Q;
  if (disc >= 0) {
    Scalar root_disc = sqrt(disc);
    Scalar cube =
        Q < 0 ? Scalar(-0.5) * (Q - root_disc) : Scalar(-0.5) * (Q + root_disc);
    Scalar u = static_cast<Scalar>(cbrt(cube));
    Scalar x = u - third_P / u;
    r[0] = x;
    return 1;
  } else {
    Scalar y3_re = Scalar(-0.5) * Q;
    Scalar y3_im = Scalar(0.5) * sqrt(-disc);
    // y = cube root (y3)
    Scalar theta = atan2(y3_im, y3_re) * third;
    Scalar radius = sqrt(-third_P);
    Scalar y_re = radius * cos(theta);
    Scalar x = 2 * y_re;
    Scalar root_disc = sqrt(-3 * x * x - 4 * P);
    if (x > 0) {
      r[0] = Scalar(-0.5) * (x + root_disc);
      r[1] = Scalar(-0.5) * (x - root_disc);
      r[2] = x;
    } else {
      r[0] = x;
      r[1] = Scalar(-0.5) * (x + root_disc);
      r[2] = Scalar(-0.5) * (x - root_disc);
    }
    return 3;
  }
}

template <typename Scalar>
inline int FindQuarticRealRoots(Scalar B, Scalar C, Scalar D, Scalar E,
                                Scalar r[]) {
  static const Scalar third = Scalar(1.0 / 3.0);
  Scalar alpha = C - Scalar(0.375) * B * B;
  Scalar half_B = Scalar(0.5) * B;
  Scalar beta = half_B * (half_B * half_B - C) + D;
  Scalar q_B = Scalar(0.25) * B;
  Scalar four_gamma = B * (q_B * (C - 3 * q_B * q_B) - D) + 4 * E;

  if (beta * beta < Scalar(1e-18)) {
    Scalar disc = alpha * alpha - four_gamma;
    if (disc < 0)
      return 0;
    Scalar root_disc = sqrt(disc);
    {
      Scalar disc2 = -alpha + root_disc;
      if (disc2 < 0)
        return 0;
      Scalar root_disc2 = sqrt(Scalar(0.5) * disc2);
      r[0] = -q_B - root_disc2;
      r[1] = -q_B + root_disc2;
    }
    {
      Scalar disc2 = -alpha - root_disc;
      if (disc2 < 0)
        return 2;
      Scalar root_disc2 = sqrt(Scalar(0.5) * disc2);
      r[2] = -q_B - root_disc2;
      r[3] = -q_B + root_disc2;
    }
    return 4;
  }

  Scalar third_alpha = alpha * third;
  Scalar P = Scalar(-0.25) * (alpha * third_alpha + four_gamma);
  Scalar Q =
      third_alpha * (Scalar(0.25) * (four_gamma - third_alpha * third_alpha)) -
      beta * beta * Scalar(0.125);

  Scalar dcr[3];
  int ndcr = DepressedCubicRealRoots(P, Q, dcr);
  Scalar disc2 = 2 * (dcr[ndcr - 1] - third_alpha);
  Scalar W = sqrt(disc2);

  Scalar disc_base = 2 * alpha + disc2; // 3*alpha + 2*y;
  Scalar disc_add = 2 * beta / W;

  int count = 0;
  if (disc_base + disc_add <= 0) {
    Scalar sr = sqrt(-disc_base - disc_add);
    r[count++] = -q_B + Scalar(0.5) * (W - sr);
    r[count++] = -q_B + Scalar(0.5) * (W + sr);
  }
  if (disc_base - disc_add <= 0) {
    Scalar sr = sqrt(-disc_base + disc_add);
    r[count++] = -q_B - Scalar(0.5) * (W + sr);
    r[count++] = -q_B - Scalar(0.5) * (W - sr);
  }
  return count;
}

template <typename Scalar>
inline slick::SE3Group<Scalar>
GetAbsoluteOrientationFrom3Points(const Eigen::Matrix<Scalar, 3, 1> x[],
                                  const Eigen::Matrix<Scalar, 3, 1> y[]) {
  Eigen::Matrix<Scalar, 3, 3> D, D1;

  D.transpose().col(0) = x[1] - x[0];
  D.transpose().col(1) = x[2] - x[0];
  D.transpose().col(2) = D.transpose().col(1).cross(D.transpose().col(0));

  D1.transpose().col(0) = y[1] - y[0];
  D1.transpose().col(1) = y[2] - y[0];
  D1.transpose().col(2) = D1.transpose().col(1).cross(D1.transpose().col(0));

  Eigen::Matrix<Scalar, 3, 3> m3 = (D.inverse() * D1).transpose();
  slick::SO3Group<Scalar> so3(m3);

  Eigen::Matrix<Scalar, 3, 1> T = y[0] - so3 * x[0];
  return slick::SE3Group<Scalar>(so3, T);
}

// The function for pose estimation from three 2D - 3D point correspondences.
// It implements the algorithm given by the Fischer and Bolles RANSAC paper,1980
// It assumes that the three points are in general position (not collinear).
// Input is an array of 3D cartesian positions and an array of
// 2D vectors that are the perspective projections of the points.
// Ouput is up to four poses satisfying the input with positive depths
// (points in front of the camera).
// @param[in] x an array containing at least 3 points
// @param[in] z an array containing the perspective projections
// of the points given by x in the current pose
// @param[out] poses the vector onto which any valid poses are appended
// @return the number of  poses appended to the vector
// @ingroup algorithms
template <typename Scalar>
inline int CalcThreePointPoses(const Eigen::Matrix<Scalar, 3, 1> xi[],
                               const Eigen::Matrix<Scalar, 2, 1> zi[],
                               std::vector<slick::SE3Group<Scalar>> &poses) {
  Scalar ab_sq = (xi[1] - xi[0]).squaredNorm();
  Scalar ac_sq = (xi[2] - xi[0]).squaredNorm();
  Scalar bc_sq = (xi[2] - xi[1]).squaredNorm();
  typedef Eigen::Matrix<Scalar, 2, 1> Vec2;
  const Vec2 t1 = zi[0];
  const Vec2 t2 = zi[1];
  const Vec2 t3 = zi[2];
  Eigen::Matrix<Scalar, 3, 1> za = unproject(t1);
  Eigen::Matrix<Scalar, 3, 1> zb = unproject(t2);
  Eigen::Matrix<Scalar, 3, 1> zc = unproject(t3);
  za.normalize();
  zb.normalize();
  zc.normalize();

  Scalar cos_ab = za.transpose() * zb;
  Scalar cos_ac = za.transpose() * zc;
  Scalar cos_bc = zb.transpose() * zc;

  Scalar K1 = bc_sq / ac_sq;
  Scalar K2 = bc_sq / ab_sq;

  Scalar p[5];
  {
    Scalar K12 = K1 * K2;
    Scalar K12_1_2 = K12 - K1 - K2;
    Scalar K2_1_K1_cab = K2 * (1 - K1) * cos_ab;

    p[4] = (square(K12_1_2) - 4 * K12 * cos_bc * cos_bc);
    p[3] = (4 * (K12_1_2)*K2_1_K1_cab +
            4 * K1 * cos_bc *
                ((K12 + K2 - K1) * cos_ac + 2 * K2 * cos_ab * cos_bc));
    p[2] = (square(2 * K2_1_K1_cab) + 2 * (K12 + K1 - K2) * K12_1_2 +
            4 * K1 *
                ((K1 - K2) * cos_bc * cos_bc + (1 - K2) * K1 * cos_ac * cos_ac -
                 2 * K2 * (1 + K1) * cos_ab * cos_ac * cos_bc));
    p[1] = (4 * (K12 + K1 - K2) * K2_1_K1_cab +
            4 * K1 * ((K12 - K1 + K2) * cos_ac * cos_bc +
                      2 * K12 * cos_ab * cos_ac * cos_ac));
    p[0] = square(K12 + K1 - K2) - 4 * K12 * K1 * cos_ac * cos_ac;
  }
  Eigen::Matrix<Scalar, 3, 1> xi1[3];

  Scalar roots[4];
  Scalar inv_p4 = static_cast<Scalar>(1.0 / p[4]);
  for (int i = 0; i < 4; ++i)
    p[i] *= inv_p4;
  int nr = FindQuarticRealRoots(p[3], p[2], p[1], p[0], roots);

  try {
    int count = 0;
    for (int i = 0; i < nr; ++i) {
      Scalar x = roots[i];
      if (x <= 0)
        continue;
      for (int j = 0; j < 3; ++j)
        x = eval_newton_quartic(p[3], p[2], p[1], p[0], x);
      Scalar xx = x * x;

      Scalar a_den = xx - 2 * x * cos_ab + 1;
      Scalar a = sqrt(ab_sq / a_den);
      Scalar b = a * x;

      Scalar M = 1 - K1;
      Scalar P = 2 * (K1 * cos_ac - x * cos_bc);
      Scalar Q = xx - K1;

      Scalar P1 = -2 * x * cos_bc;
      Scalar Q1 = xx - K2 * a_den;

      Scalar den = M * Q1 - Q;
      if (den == 0) {
        std::cerr << "skipped" << std::endl;
        continue;
      }

      Scalar y = (P1 * Q - P * Q1) / den;
      Scalar c = a * y;

      xi1[0] = a * za;
      xi1[1] = b * zb;
      xi1[2] = c * zc;

      poses.push_back(GetAbsoluteOrientationFrom3Points(xi, xi1));
      ++count;
    }
    return count;
  } catch (std::exception &e) {
    std::cerr << "Threepointpose.h:: Oops: " << e.what() << std::endl;
    return 0;
  }
}

template <typename Scalar>
inline std::pair<bool, slick::SE3Group<Scalar>>
    ComputeRobustAbsolutPoseRANSACMLE(std::vector<
        std::pair<Eigen::Matrix<Scalar, 4, 1>, Eigen::Matrix<Scalar, 2, 1>>> &
                                          observations) {
  using Observations = std::vector<
      std::pair<Eigen::Matrix<Scalar, 4, 1>, Eigen::Matrix<Scalar, 2, 1>>>;
  using PotentialPoses = std::vector<slick::SE3Group<Scalar>>;
  PotentialPoses potential_poses;
  Eigen::Matrix<Scalar, 3, 1> world_points[3];
  Eigen::Matrix<Scalar, 2, 1> observed_points[3];
  int nrand = 0;
  for (size_t i = 0; i < 100; i++) {
    std::random_shuffle(observations.begin(), observations.end());
    nrand = std::rand() % observations.size();
    if (nrand - 3 > 0) {
      world_points[0] = slick::project(observations[nrand].first);
      world_points[1] = slick::project(observations[nrand - 1].first);
      world_points[2] = slick::project(observations[nrand - 2].first);
      observed_points[0] = observations[nrand].second;
      observed_points[1] = observations[nrand - 1].second;
      observed_points[2] = observations[nrand - 2].second;
#if 0
      for (int i1 = 0; i1 < 3; i1++) {
        std::cout << "3dpoints=[" << vWPoint3s[0].transpose()
                  << ";" << vWPoint3s[1].transpose()
                  << ";" << vWPoint3s[2].transpose() << "]\n";
        std::cout << "2dpoints=[" << vObservedPoint2s[0].transpose()
                  << ";" << vObservedPoint2s[1].transpose()
                  << ";" << vObservedPoint2s[2].transpose() << "]\n";
      }
#endif
      try {
#if 0
        std::vector<Eigen::Matrix<SlickScalar, 3, 3>> Rs;
        std::vector<Eigen::Matrix<SlickScalar, 3, 1>> Ts;
        slick::computeAbsolutePose3Point(
              vObservedPoint2s[0], vObservedPoint2s[1], vObservedPoint2s[1],
              vWPoint3s[0], vWPoint3s[1], vWPoint3s[2], Rs, Ts);
#else
        std::vector<slick::SE3Group<Scalar>> vPoses;
        CalcThreePointPoses(world_points, observed_points, vPoses);
#endif
        bool bValid = false;
        for (size_t i = 0; i < vPoses.size(); ++i) {
          slick::SE3Group<Scalar> se3Pose = vPoses[i];
          bValid = false;
          Eigen::Matrix<Scalar, 3, 1> v3Cam =
              se3Pose * slick::project(observations[nrand].first);
          Eigen::Matrix<Scalar, 2, 1> v2Error =
              observations[nrand].second - slick::project(v3Cam);
          if (v3Cam[2] > 0 && v2Error.squaredNorm() < 1.e-9) {
            potential_poses.push_back(se3Pose);
          } else {
            std::cout << "Invalid Pose.........."
                      << "z:" << v3Cam[2] << " error:" << v2Error.squaredNorm()
                      << std::endl;
            return std::make_pair(false, SE3Group<Scalar>());
          }
        }
      } catch (std::exception &e) {
        std::cout << "Error while calculating three_point_pose!\n";
        std::cout << e.what();
      }
    }
  }
  // get the best pose
  Scalar min_errorsq = 1.e+9;
  Scalar sum_errorsq = 0.0;
  slick::SE3Group<Scalar> se3_best;
  slick::Tukey<Scalar> est;
  for (size_t i = 0; i < potential_poses.size(); i++) {
    std::vector<Scalar> vErrorSq;
    vErrorSq.reserve(observations.size());
    for (size_t j = 0; j < observations.size(); j++) {
      Eigen::Matrix<Scalar, 2, 1> v2Error =
          observations[j].second -
          slick::project(potential_poses[i] *
                         slick::project(observations[j].first));
      vErrorSq.push_back(v2Error.squaredNorm());
    }

    est.compute_sigma_squared(vErrorSq);

    sum_errorsq = 0.0;
    for (size_t j = 0; j < observations.size(); j++) {
      Eigen::Matrix<Scalar, 2, 1> v2Error =
          observations[j].second -
          slick::project(potential_poses[i] *
                         slick::project(observations[j].first));
      Scalar dErrorSq = v2Error.squaredNorm();
      const Scalar w = est.weight(dErrorSq);
      sum_errorsq += w * dErrorSq;
    }
    if (min_errorsq > sum_errorsq) {
      min_errorsq = sum_errorsq;
      se3_best = potential_poses[i];
    }
  }
  return std::make_pair(true, se3_best);
}
} // namespace slick