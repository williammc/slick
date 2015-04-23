#pragma once
#include <memory>
#include <random>
#include <gtest/gtest.h>

#define NSAMPLES 1000

namespace slick {

template <typename T>
inline std::pair<std::shared_ptr<std::mt19937>,
          std::shared_ptr<std::uniform_real_distribution<T>>>
sample_uniform_real(const T min, const T max) {
  // static_assert(min < max, "Min value must be less than Max value");
  std::random_device rd;
  std::shared_ptr<std::mt19937> gen(new std::mt19937(rd()));
  std::shared_ptr<std::uniform_real_distribution<T>> dis(
      new std::uniform_real_distribution<T>(min, max));
  return std::make_pair(gen, dis);
}

template <typename Derived1, typename Derived2>
inline void EXPECT_MATRIX_EQUAL(const Eigen::MatrixBase<Derived1>& left,
                         const Eigen::MatrixBase<Derived2>& right) {
  for (int i = 0; i < left.rows(); ++i) {
    for (int j = 0; j < left.cols(); ++j) {
      if (std::is_same<float, typename Derived1::Scalar>::value) {
        EXPECT_FLOAT_EQ(left(i, j), right(i, j));
      } else if (std::is_same<double, typename Derived1::Scalar>::value) {
        EXPECT_DOUBLE_EQ(left(i, j), right(i, j));
      }
    }
  }
}

template <typename Derived1, typename Derived2, typename P>
inline void EXPECT_MATRIX_NEAR(const Eigen::MatrixBase<Derived1>& left,
                        const Eigen::MatrixBase<Derived2>& right,
                        const P& gap) {
  for (int i = 0; i < left.rows(); ++i) {
    for (int j = 0; j < left.cols(); ++j) {
      EXPECT_NEAR(left(i, j), right(i, j), gap);
    }
  }
}

template <typename T>
inline T Gap() {
  if (typeid(T) == typeid(double)) return T(1.e-9);
  if (typeid(T) == typeid(float)) return T(1.e-6);
  return T(0);
}

template <typename T> inline T GenRandNumber(const int N) {
  return (T(std::rand()) / T(RAND_MAX) - T(0.5)) * T(N);
}

template <typename T> inline Eigen::Matrix<T, 2, 1> GenRandPoint2D(const int N) {
  return Eigen::Matrix<T, 2, 1>(GenRandNumber<T>(N), GenRandNumber<T>(N));
}

template <typename T> inline Eigen::Matrix<T, 3, 1> GenRandPoint(const int N) {
  return Eigen::Matrix<T, 3, 1>(GenRandNumber<T>(N), GenRandNumber<T>(N), GenRandNumber<T>(N));
}
}  // namespace slick
