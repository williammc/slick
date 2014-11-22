#ifndef LOOK3D_MATH_UNITTEST_UTIL_H_
#define LOOK3D_MATH_UNITTEST_UTIL_H_
#include <memory>
#include <random>
#include <gtest/gtest.h>

#define NSAMPLES 1000

namespace look3d {

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
  if (typeid(T) == typeid(double)) return double(1.e-9);
  if (typeid(T) == typeid(float)) return float(1.e-6);
}
}  // namespace look3d
#endif  // LOOK3D_MATH_UNITTEST_UTIL_H_
