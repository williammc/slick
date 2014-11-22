#include <ctime>
#include <gtest/gtest.h>
#include <Eigen/Eigen>
#include <Eigen/Geometry>

#include "slick/math/se2.h"
#include "slick/test/unittest/util.h"

using namespace std;
using namespace Eigen;
using namespace slick;

// constructors tests ==========================================================
template <typename T>
void DefaultConstructor_Test() {
  // for double precision
  SE2Group<T> trans;
  EXPECT_MATRIX_NEAR(Eigen::Matrix<T, 2, 2>::Identity(), trans.get_rotation().get_matrix(),
                     Gap<T>());
  EXPECT_MATRIX_NEAR(Eigen::Matrix<T, 2, 1>::Zero(), trans.get_translation(),
                     Gap<T>());
}

TEST(SE2Test, DefaultConstructor) {
  DefaultConstructor_Test<double>();
  DefaultConstructor_Test<float>();
}

template <typename T>
void FromSO2VectorConstructor_Test() {
  std::srand(std::time(0));  // use current time as seed for random generator
  auto angle = T(std::rand());

  Eigen::Matrix<T, 2, 2> mrot;
  mrot = Eigen::Rotation2D<T>(angle);

  Eigen::Matrix<T, 2, 1> t = Eigen::Matrix<T, 2, 1>::Random();

  SO2Group<T> rot(angle);
  SE2Group<T> trans(rot, t);

  EXPECT_MATRIX_NEAR(mrot, trans.get_rotation().get_matrix(),
                     Gap<T>());
  EXPECT_MATRIX_NEAR(t, trans.get_translation(),
                     Gap<T>());
}

TEST(SE2Test, FromSO2VectorConstructor) {
  FromSO2VectorConstructor_Test<double>();
  FromSO2VectorConstructor_Test<float>();
}


template <typename T>
void FromVectorConstructor_Test() {
  std::srand(std::time(0));  // use current time as seed for random generator
  auto angle = T(std::rand());

  Eigen::Matrix<T, 2, 2> mrot;
  mrot = Eigen::Rotation2D<T>(angle);

  Eigen::Matrix<T, 2, 1> t = Eigen::Matrix<T, 2, 1>::Random();

  SE2Group<T> trans(Eigen::Matrix<T, 3, 1>(t[0], t[1], angle));

  EXPECT_MATRIX_NEAR(mrot, trans.get_rotation().get_matrix(),
                     Gap<T>());
}

TEST(SE2Test, FromVectorConstructor) {
  FromVectorConstructor_Test<double>();
  FromVectorConstructor_Test<float>();
}


// template <typename T>
// void FromListInitializerConstructor_Test() {
//   SE2Group<T> rot{1, 0, 0, 1};
//   EXPECT_MATRIX_NEAR(Eigen::Matrix<T, 2, 2>::Identity(), rot.get_matrix(),
//                      Gap<T>());

//   EXPECT_NEAR(0.0, rot.ln(), Gap<T>());
// }

// TEST(SE2Test, FromListInitializerConstructor) {
//   FromListInitializerConstructor_Test<double>();
//   FromListInitializerConstructor_Test<float>();
// }

// se2 specific functions ======================================================
template <typename T>
void Ln_Test() {
  std::srand(std::time(0));  // use current time as seed for random generator
  auto angle = T(std::rand()) / T(RAND_MAX) * 2 * M_PI;
  Eigen::Matrix<T, 2, 1> t = Eigen::Matrix<T, 2, 1>::Random();
  SE2Group<T> trans(Eigen::Matrix<T, 3, 1>(t[0], t[1], angle));

  auto est = trans.ln();
  est[2] = (est[2] > 0) ? est[2] : 2 * M_PI + est[2];
  EXPECT_NEAR(angle, est[2], Gap<T>());
}

TEST(SE2Test, Ln) {
  Ln_Test<double>();
  Ln_Test<float>();
}

template <typename T>
void Inverse_Test() {
  std::srand(std::time(0));  // use current time as seed for random generator
  auto angle = T(std::rand());
  Eigen::Matrix<T, 2, 1> t = Eigen::Matrix<T, 2, 1>::Random();
  SE2Group<T> trans(Eigen::Matrix<T, 3, 1>(t[0], t[1], angle));
  Eigen::Matrix<T, 2, 2> m_inv = trans.get_rotation().get_matrix().inverse();
  EXPECT_MATRIX_NEAR(m_inv, trans.inverse().get_rotation().get_matrix(), Gap<T>());
  // EXPECT_MATRIX_NEAR(-t, trans.inverse().get_translation(), Gap<T>());
}

TEST(SE2Test, Inverse) {
  Inverse_Test<double>();
  Inverse_Test<float>();
}

template <typename T>
void SE2RightHandMulOperator_Test() {
  std::srand(std::time(0));  // use current time as seed for random generator
  auto angle = T(std::rand());

  Eigen::Matrix<T, 3, 3> mtrans1 = Eigen::Matrix<T, 3, 3>::Identity();
  Eigen::Matrix<T, 2, 2> m;
  m = Eigen::Rotation2D<T>(angle);
  mtrans1.block(0, 0, 2, 2) = m;
  Eigen::Matrix<T, 2, 1> t = Eigen::Matrix<T, 2, 1>::Random();
  SE2Group<T> trans1(Eigen::Matrix<T, 3, 1>(t[0], t[1], angle));
  mtrans1.block(0, 2, 2, 1) = trans1.get_translation();

  std::srand(std::time(0));  // use current time as seed for random generator
  angle = T(std::rand());
  Eigen::Matrix<T, 3, 3> mtrans2 = Eigen::Matrix<T, 3, 3>::Identity();
  m = Eigen::Rotation2D<T>(angle);
  mtrans2.block(0, 0, 2, 2) = m;
  t = Eigen::Matrix<T, 2, 1>::Random();
  SE2Group<T> trans2(Eigen::Matrix<T, 3, 1>(t[0], t[1], angle));
  mtrans2.block(0, 2, 2, 1) = trans2.get_translation();

  EXPECT_MATRIX_NEAR((mtrans1 * mtrans2).block(0, 0, 2, 3), 
    (trans1 * trans2).get_matrix(), Gap<T>());
}

TEST(SE2Test, SE2RightHandMulOperator) {
  SE2RightHandMulOperator_Test<double>();
  SE2RightHandMulOperator_Test<float>();
}

template <typename T>
void MatrixRightHandMulOperator_Test() {
  std::srand(std::time(0));  // use current time as seed for random generator
  auto angle = T(std::rand());

  Eigen::Matrix<T, 3, 3> mtrans = Eigen::Matrix<T, 3, 3>::Identity();
  Eigen::Matrix<T, 2, 2> m;
  m = Eigen::Rotation2D<T>(angle);
  mtrans.block(0, 0, 2, 2) = m;
  Eigen::Matrix<T, 2, 1> t = Eigen::Matrix<T, 2, 1>::Random();
  SE2Group<T> trans(Eigen::Matrix<T, 3, 1>(t[0], t[1], angle));
  mtrans.block(0, 2, 2, 1) = trans.get_translation();

  Eigen::Matrix<T, 3, 3> mrand = Eigen::Matrix<T, 3, 3>::Random();
  EXPECT_MATRIX_NEAR(mtrans * mrand, trans * mrand, Gap<T>());
}

TEST(SE2Test, MatrixRightHandMulOperator) {
  MatrixRightHandMulOperator_Test<double>();
  MatrixRightHandMulOperator_Test<float>();
}

template <typename T>
void MatrixLeftHandMulOperator_Test() {
  std::srand(std::time(0));  // use current time as seed for random generator
  auto angle = T(std::rand());

  Eigen::Matrix<T, 3, 3> mtrans = Eigen::Matrix<T, 3, 3>::Identity();
  Eigen::Matrix<T, 2, 2> m;
  m = Eigen::Rotation2D<T>(angle);
  mtrans.block(0, 0, 2, 2) = m;
  Eigen::Matrix<T, 2, 1> t = Eigen::Matrix<T, 2, 1>::Random();
  SE2Group<T> trans(Eigen::Matrix<T, 3, 1>(t[0], t[1], angle));
  mtrans.block(0, 2, 2, 1) = trans.get_translation();

  Eigen::Matrix<T, 3, 3> mrand = Eigen::Matrix<T, 3, 3>::Random();
  EXPECT_MATRIX_NEAR(mrand*mtrans , mrand*trans, Gap<T>());
}

TEST(SE2Test, MatrixLeftHandMulOperator) {
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
