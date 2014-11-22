#include <ctime>
#include <gtest/gtest.h>
#include <Eigen/Eigen>
#include <Eigen/Geometry>

#include "math/se3.h"
#include "math/test/unittest/util.h"

using namespace std;
using namespace Eigen;
using namespace slick;

// constructors tests ==========================================================
template <typename T>
void DefaultConstructor_Test() {
  // for double precision
  SE3Group<T> trans;
  EXPECT_MATRIX_NEAR(Eigen::Matrix<T, 2, 2>::Identity(), trans.get_rotation().get_matrix(),
                     Gap<T>());
  EXPECT_MATRIX_NEAR(Eigen::Matrix<T, 2, 1>::Zero(), trans.get_translation(),
                     Gap<T>());
}

TEST(SE3Test, DefaultConstructor) {
  DefaultConstructor_Test<double>();
  DefaultConstructor_Test<float>();
}

template <typename T>
void FromSO3VectorConstructor_Test() {
  Eigen::Matrix<T, 3, 1> aaxis = Eigen::Matrix<T, 3, 1>::Random();
  aaxis.normalize();
  std::srand(std::time(0));  // use current time as seed for random generator
  auto angle = T(std::rand());

  aaxis *= angle;
  SO3Group<T> rot(aaxis);
  Eigen::Matrix<T, 3, 1> t = Eigen::Matrix<T, 3, 1>::Random();
  SE3Group<T> trans(rot, t);

  EXPECT_MATRIX_NEAR(rot.get_matrix(), trans.get_rotation().get_matrix(),
                     Gap<T>());
  EXPECT_MATRIX_NEAR(t, trans.get_translation(),
                     Gap<T>());
}

TEST(SE3Test, FromSO3VectorConstructor) {
  FromSO3VectorConstructor_Test<double>();
  FromSO3VectorConstructor_Test<float>();
}


template <typename T>
void FromVectorConstructor_Test() {
  Eigen::Matrix<T, 3, 1> aaxis = Eigen::Matrix<T, 3, 1>::Random();
  aaxis.normalize();
  std::srand(std::time(0));  // use current time as seed for random generator
  auto angle = T(std::rand());

  aaxis *= angle;
  SO3Group<T> rot(aaxis);

  Eigen::Matrix<T, 6, 1> v;
  v.segment(0, 3) = Eigen::Matrix<T, 3, 1>::Random();
  v.segment(3, 3) = aaxis;
  SE3Group<T> trans(v);

  EXPECT_MATRIX_NEAR(rot.get_matrix(), trans.get_rotation().get_matrix(),
                     Gap<T>());
}

TEST(SE3Test, FromVectorConstructor) {
  FromVectorConstructor_Test<double>();
  FromVectorConstructor_Test<float>();
}


// template <typename T>
// void FromListInitializerConstructor_Test() {
//   SE3Group<T> rot{1, 0, 0, 1};
//   EXPECT_MATRIX_NEAR(Eigen::Matrix<T, 2, 2>::Identity(), rot.get_matrix(),
//                      Gap<T>());

//   EXPECT_NEAR(0.0, rot.ln(), Gap<T>());
// }

// TEST(SE3Test, FromListInitializerConstructor) {
//   FromListInitializerConstructor_Test<double>();
//   FromListInitializerConstructor_Test<float>();
// }

// se2 specific functions ======================================================
template <typename T>
void Ln_Test() {
  Eigen::Matrix<T, 3, 1> aaxis = Eigen::Matrix<T, 3, 1>::Random();
  aaxis.normalize();
  std::srand(std::time(0));  // use current time as seed for random generator
  auto angle = T(std::rand()) / T(RAND_MAX) * 2 * M_PI;
  aaxis *= angle;
  SO3Group<T> rot(aaxis);
  Eigen::Matrix<T, 3, 1> t = Eigen::Matrix<T, 3, 1>::Random();
  SE3Group<T> trans(rot, t);

  Eigen::Matrix<T, 3, 1> est_axis = trans.ln().segment(3, 3);
  T d = est_axis.dot(aaxis);
  if (d < 0) {
    auto a = -est_axis.norm() + 2 * M_PI;
    est_axis = -est_axis.normalized() * a;
  }
  EXPECT_MATRIX_NEAR(aaxis, est_axis, Gap<T>());
}

TEST(SE3Test, Ln) {
  Ln_Test<double>();
  Ln_Test<float>();
}

template <typename T>
void Inverse_Test() {
  Eigen::Matrix<T, 6, 1> v = Eigen::Matrix<T, 6, 1>::Random();
  SE3Group<T> trans(v);
  Eigen::Matrix<T, 3, 3> m_inv = trans.get_rotation().get_matrix().inverse();
  EXPECT_MATRIX_NEAR(m_inv, trans.inverse().get_rotation().get_matrix(), Gap<T>());
  // EXPECT_MATRIX_NEAR(-t, trans.inverse().get_translation(), Gap<T>());
}

TEST(SE3Test, Inverse) {
  Inverse_Test<double>();
  Inverse_Test<float>();
}

template <typename T>
void SE3RightHandMulOperator_Test() {
  Eigen::Matrix<T, 6, 1> v = Eigen::Matrix<T, 6, 1>::Random();
  SE3Group<T> trans1(v);
  Eigen::Matrix<T, 4, 4> mtrans1 = Eigen::Matrix<T, 4, 4>::Identity();
  mtrans1.block(0, 0, 3, 4) = trans1.get_matrix();

  v = Eigen::Matrix<T, 6, 1>::Random();
  SE3Group<T> trans2(v);
  Eigen::Matrix<T, 4, 4> mtrans2 = Eigen::Matrix<T, 4, 4>::Identity();
  mtrans2.block(0, 0, 3, 4) = trans2.get_matrix();

  EXPECT_MATRIX_NEAR((mtrans1 * mtrans2).block(0, 0, 3, 4), 
    (trans1 * trans2).get_matrix(), Gap<T>());
}

TEST(SE3Test, SE3RightHandMulOperator) {
  SE3RightHandMulOperator_Test<double>();
  SE3RightHandMulOperator_Test<float>();
}

template <typename T>
void MatrixRightHandMulOperator_Test() {
  Eigen::Matrix<T, 6, 1> v = Eigen::Matrix<T, 6, 1>::Random();
  SE3Group<T> trans(v);
  Eigen::Matrix<T, 4, 4> mtrans = Eigen::Matrix<T, 4, 4>::Identity();
  mtrans.block(0, 0, 3, 4) = trans.get_matrix();

  Eigen::Matrix<T, 4, 4> mrand = Eigen::Matrix<T, 4, 4>::Random();
  EXPECT_MATRIX_NEAR(mtrans * mrand, trans * mrand, Gap<T>());
}

TEST(SE3Test, MatrixRightHandMulOperator) {
  MatrixRightHandMulOperator_Test<double>();
  MatrixRightHandMulOperator_Test<float>();
}

template <typename T>
void MatrixLeftHandMulOperator_Test() {
  Eigen::Matrix<T, 6, 1> v = Eigen::Matrix<T, 6, 1>::Random();
  SE3Group<T> trans(v);
  Eigen::Matrix<T, 4, 4> mtrans = Eigen::Matrix<T, 4, 4>::Identity();
  mtrans.block(0, 0, 3, 4) = trans.get_matrix();

  Eigen::Matrix<T, 4, 4> mrand = Eigen::Matrix<T, 4, 4>::Random();
  EXPECT_MATRIX_NEAR(mrand*mtrans , mrand*trans, Gap<T>());
}

TEST(SE3Test, MatrixLeftHandMulOperator) {
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
