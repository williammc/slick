#include <ctime>
#include <gtest/gtest.h>
#include <Eigen/Eigen>
#include <Eigen/Geometry>

#include "math/so2.h"
#include "math/test/unittest/util.h"

using namespace std;
using namespace Eigen;
using namespace slick;

// constructors tests ==========================================================
template <typename T>
void DefaultConstructor_Test() {
  // for double precision
  SO2Group<T> rot;
  EXPECT_MATRIX_NEAR(Eigen::Matrix<T, 2, 2>::Identity(), rot.get_matrix(),
                     Gap<T>());
}

TEST(SO2Test, DefaultConstructor) {
  DefaultConstructor_Test<double>();
  DefaultConstructor_Test<float>();
}

template <typename T>
void FromMatrixConstructor_Test() {
  std::srand(std::time(0));  // use current time as seed for random generator
  auto angle = T(std::rand());

  Eigen::Matrix<T, 2, 2> mrot;
  mrot = Eigen::Rotation2D<T>(angle);
  SO2Group<T> rot(mrot);
  EXPECT_MATRIX_NEAR(mrot, rot.get_matrix(), Gap<T>());
}

TEST(SO2Test, FromMatrixConstructor) {
  FromMatrixConstructor_Test<double>();
  FromMatrixConstructor_Test<float>();
}

template <typename T>
void FromAngleConstructor_Test() {
  std::srand(std::time(0));  // use current time as seed for random generator
  auto angle = T(std::rand());

  Eigen::Matrix<T, 2, 2> mrot;
  mrot = Eigen::Rotation2D<T>(angle);
  SO2Group<T> rot(angle);
  EXPECT_MATRIX_NEAR(mrot, rot.get_matrix(), Gap<T>());
}

TEST(SO2Test, FromAngleConstructor) {
  FromAngleConstructor_Test<double>();
  FromAngleConstructor_Test<float>();
}

template <typename T>
void FromListInitializerConstructor_Test() {
  SO2Group<T> rot{1, 0, 0, 1};
  EXPECT_MATRIX_NEAR(Eigen::Matrix<T, 2, 2>::Identity(), rot.get_matrix(),
                     Gap<T>());

  EXPECT_NEAR(0.0, rot.ln(), Gap<T>());
}

TEST(SO2Test, FromListInitializerConstructor) {
  FromListInitializerConstructor_Test<double>();
  FromListInitializerConstructor_Test<float>();
}

// so2 specific functions ======================================================
template <typename T>
void Ln_Test() {
  std::srand(std::time(0));  // use current time as seed for random generator
  auto angle = T(std::rand()) / T(RAND_MAX) * 2 * M_PI;
  SO2Group<T> rot(angle);
  auto est = rot.ln();
  est = (est > 0) ? est : 2 * M_PI + est;
  EXPECT_NEAR(angle, est, Gap<T>());
}

TEST(SO2Test, Ln) {
  Ln_Test<double>();
  Ln_Test<float>();
}

template <typename T>
void Inverse_Test() {
  std::srand(std::time(0));  // use current time as seed for random generator
  auto angle = T(std::rand());
  SO2Group<T> rot(angle);
  Eigen::Matrix<T, 2, 2> m_inv = rot.get_matrix().inverse();
  EXPECT_MATRIX_NEAR(m_inv, rot.inverse().get_matrix(), Gap<T>());
}

TEST(SO2Test, Inverse) {
  Inverse_Test<double>();
  Inverse_Test<float>();
}

template <typename T>
void SO2RightHandMulOperator_Test() {
  std::srand(std::time(0));  // use current time as seed for random generator
  auto angle = T(std::rand());

  Eigen::Matrix<T, 2, 2> mrot1;
  mrot1 = Eigen::Rotation2D<T>(angle);
  SO2Group<T> rot1(angle);

  std::srand(std::time(0));  // use current time as seed for random generator
  angle = T(std::rand());
  Eigen::Matrix<T, 2, 2> mrot2;
  mrot2 = Eigen::Rotation2D<T>(angle);
  SO2Group<T> rot2(angle);

  SO2Group<T> rot(2.0);
  EXPECT_MATRIX_NEAR(mrot1 * mrot2, (rot1 * rot2).get_matrix(), Gap<T>());
}

TEST(SO2Test, SO2RightHandMulOperator) {
  SO2RightHandMulOperator_Test<double>();
  SO2RightHandMulOperator_Test<float>();
}

template <typename T>
void MatrixRightHandMulOperator_Test() {
  std::srand(std::time(0));  // use current time as seed for random generator
  auto angle = T(std::rand());

  Eigen::Matrix<T, 2, 2> mrot;
  mrot = Eigen::Rotation2D<T>(angle);

  SO2Group<T> rot(angle);

  Eigen::Matrix<T, 2, 2> mrand = Eigen::Matrix<T, 2, 2>::Random();
  EXPECT_MATRIX_NEAR(mrot * mrand, rot * mrand, Gap<T>());
}

TEST(SO2Test, MatrixRightHandMulOperator) {
  MatrixRightHandMulOperator_Test<double>();
  MatrixRightHandMulOperator_Test<float>();
}

template <typename T>
void MatrixLeftHandMulOperator_Test() {
  std::srand(std::time(0));  // use current time as seed for random generator
  auto angle = T(std::rand());

  Eigen::Matrix<T, 2, 2> mrot;
  mrot = Eigen::Rotation2D<T>(angle);

  SO2Group<T> rot(angle);

  Eigen::Matrix<T, 2, 2> mrand = Eigen::Matrix<T, 2, 2>::Random();
  EXPECT_MATRIX_NEAR(mrand * mrot, mrand * rot, Gap<T>());
}

TEST(SO2Test, MatrixLeftHandMulOperator) {
  MatrixLeftHandMulOperator_Test<double>();
  MatrixLeftHandMulOperator_Test<float>();
}

int main(int argc, char** argv) {
  std::vector<char*> vars(argc+1);
  for (int i = 0; i < argc; ++i)
    vars[i] = argv[i];
  char ca[50];
  sprintf(ca, "--gtest_repeat=%d", NSAMPLES);
  vars[argc] = ca;
  argc++;
  ::testing::InitGoogleTest(&argc, &vars.front());
  return RUN_ALL_TESTS();
}
