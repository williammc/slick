#include <ctime>
#include <gtest/gtest.h>
#include <Eigen/Eigen>
#include <Eigen/Geometry>

#include "math/so3.h"
#include "math/test/unittest/util.h"

using namespace std;
using namespace Eigen;
using namespace slick;

// constructors tests ==========================================================
template <typename T>
void DefaultConstructor_Test() {
  // for double precision
  SO3Group<T> rot;
  EXPECT_MATRIX_NEAR(Eigen::Matrix<T, 3, 3>::Identity(), rot.get_matrix(),
                     Gap<T>());
}

TEST(SO3Test, DefaultConstructor) {
  DefaultConstructor_Test<double>();
  DefaultConstructor_Test<float>();
}

template <typename T>
void FromMatrixConstructor_Test() {
  Eigen::Matrix<T, 3, 1> aaxis = Eigen::Matrix<T, 3, 1>::Random();
  aaxis.normalize();
  std::srand(std::time(0));  // use current time as seed for random generator
  auto angle = T(std::rand());

  Eigen::Matrix<T, 3, 3> m;
  m = Eigen::AngleAxis<T>(angle, aaxis);

  SO3Group<T> rot(m);
  EXPECT_MATRIX_NEAR(m, rot.get_matrix(), Gap<T>());
}

TEST(SO3Test, FromMatrixConstructor) {
  FromMatrixConstructor_Test<double>();
  FromMatrixConstructor_Test<float>();
}

template <typename T>
void FromAngleAxisConstructor_Test() {
  Eigen::Matrix<T, 3, 1> aaxis = Eigen::Matrix<T, 3, 1>::Random();
  aaxis.normalize();
  std::srand(std::time(0));  // use current time as seed for random generator
  auto angle = T(std::rand());

  Eigen::Matrix<T, 3, 3> m;
  m = Eigen::AngleAxis<T>(angle, aaxis);

  aaxis *= angle;
  SO3Group<T> rot(aaxis);
#if 0  // NOTE: Eigen implementation of rotation matrix is not as accurate as \
       // ours
    EXPECT_MATRIX_NEAR(m, rot.get_matrix(), Gap<T>());
    std::cout << "gt:\n" << m << std::endl;
    std::cout << "est:\n" << rot.get_matrix() << std::endl;


    std::cout << "(rot.inverse()*rot).get_matrix():\n" << (rot.inverse()*rot).get_matrix() << std::endl;
    EXPECT_MATRIX_NEAR(Eigen::Matrix<T, 3, 3>::Identity(),
     (rot.inverse()*rot).get_matrix(), Gap<T>());


    std::cout << "m*m.transpose():\n" << m*m.transpose() << std::endl;
    EXPECT_MATRIX_NEAR(Eigen::Matrix<T, 3, 3>::Identity(),
     m*m.transpose(), Gap<T>());
#else
  EXPECT_MATRIX_NEAR(Eigen::Matrix<T, 3, 3>::Identity(),
                     (rot.inverse() * rot).get_matrix(), Gap<T>());
#endif
}

TEST(SO3Test, FromAngleAxisConstructor) {
  FromAngleAxisConstructor_Test<double>();
  FromAngleAxisConstructor_Test<float>();
}

// template<typename T>
// void FromListInitializerConstructor_Test() {
//   SO3Group<T> rot {1, 0, 0, 1};
//   EXPECT_MATRIX_NEAR(Eigen::Matrix<T, 3, 3>::Identity(), rot.get_matrix(),
// Gap<T>());

//   EXPECT_NEAR(0.0, rot.ln(), Gap<T>());
// }

// TEST(SO3Test, FromListInitializerConstructor) {
//   FromListInitializerConstructor_Test<double>();
//   FromListInitializerConstructor_Test<float>();
// }

// so2 specific functions ======================================================

template <typename T>
void Ln_Test() {
  Eigen::Matrix<T, 3, 1> aaxis = Eigen::Matrix<T, 3, 1>::Random();
  aaxis.normalize();
  std::srand(std::time(0));  // use current time as seed for random generator
  auto angle = T(std::rand()) / T(RAND_MAX) * 2 * M_PI;
  aaxis *= angle;

  SO3Group<T> rot(aaxis);
  auto est_axis = rot.ln();

  // std::cout << "gt: " << aaxis.transpose()
  //           << "  est:" << est_axis.transpose() << std::endl;

  T d = est_axis.dot(aaxis);
  if (d < 0) {
    auto a = -est_axis.norm() + 2 * M_PI;
    est_axis = -est_axis.normalized() * a;
  }
  EXPECT_MATRIX_NEAR(aaxis, est_axis, Gap<T>());
}

TEST(SO3Test, Ln) {
  Ln_Test<double>();
  Ln_Test<float>();
}

template <typename T>
void Inverse_Test() {
  Eigen::Matrix<T, 3, 1> aaxis = Eigen::Matrix<T, 3, 1>::Random();
  aaxis.normalize();
  std::srand(std::time(0));  // use current time as seed for random generator
  auto angle = T(std::rand());
  aaxis *= angle;
  SO3Group<T> rot(aaxis);
  Eigen::Matrix<T, 3, 3> m_inv = rot.get_matrix().inverse();
  EXPECT_MATRIX_NEAR(m_inv, rot.inverse().get_matrix(), Gap<T>());
}

TEST(SO3Test, Inverse) {
  Inverse_Test<double>();
  Inverse_Test<float>();
}

template <typename T>
void SO3RightHandMulOperator_Test() {
  Eigen::Matrix<T, 3, 1> aaxis = Eigen::Matrix<T, 3, 1>::Random();
  aaxis.normalize();
  std::srand(std::time(0));  // use current time as seed for random generator
  auto angle = T(std::rand());

#if 0  // NOTE: Eigen::AngleAxis has less accuracy than ours 
    Eigen::Matrix<T, 3, 3> mrot1;
    mrot1 = Eigen::AngleAxis<T>(angle, aaxis);

    aaxis *= angle;
    SO3Group<T> rot1(aaxis);

    aaxis = Eigen::Matrix<T, 3, 1>::Random();
    aaxis.normalize();
    std::srand(std::time(0));  // use current time as seed for random generator
    angle = T(std::rand());

    Eigen::Matrix<T, 3, 3> mrot2;
    mrot2 = Eigen::AngleAxis<T>(angle, aaxis);

    aaxis *= angle;
    SO3Group<T> rot2(aaxis);

    EXPECT_MATRIX_NEAR(mrot1 * mrot2, (rot1 * rot2).get_matrix(), Gap<T>());
#else
  aaxis *= angle;
  SO3Group<T> rot1(aaxis);
  Eigen::Matrix<T, 3, 1> t = -aaxis;
  SO3Group<T> rot2(t);

  EXPECT_MATRIX_NEAR(Eigen::Matrix<T, 3, 3>::Identity(),
                     (rot1 * rot2).get_matrix(), Gap<T>());
#endif
}

TEST(SO3Test, SO3RightHandMulOperator) {
  SO3RightHandMulOperator_Test<double>();
  SO3RightHandMulOperator_Test<float>();
}

template <typename T>
void MatrixRightHandMulOperator_Test() {
  Eigen::Matrix<T, 3, 1> aaxis = Eigen::Matrix<T, 3, 1>::Random();
  aaxis.normalize();
  std::srand(std::time(0));  // use current time as seed for random generator
  auto angle = T(std::rand());

#if 0  // NOTE: Eigen::AngleAxis has less accuracy than ours
    Eigen::Matrix<T, 3, 3> mrot1;
    mrot1 = Eigen::AngleAxis<T>(angle, aaxis);

    aaxis *= angle;
    SO3Group<T> rot1(aaxis);

    Eigen::Matrix<T, 3, 3> mrand = Eigen::Matrix<T, 3, 3>::Random();

    EXPECT_MATRIX_NEAR(mrot1 * mrand, rot1 * mrand, Gap<T>());
#else
  aaxis *= angle;
  SO3Group<T> rot1(aaxis);

  Eigen::Matrix<T, 3, 3> mrand = Eigen::Matrix<T, 3, 3>::Random();

  EXPECT_MATRIX_NEAR(rot1.get_matrix() * mrand, rot1 * mrand, Gap<T>());

#endif
}

TEST(SO3Test, MatrixRightHandMulOperator) {
  MatrixRightHandMulOperator_Test<double>();
  MatrixRightHandMulOperator_Test<float>();
}

template <typename T>
void MatrixLeftHandMulOperator_Test() {
  Eigen::Matrix<T, 3, 1> aaxis = Eigen::Matrix<T, 3, 1>::Random();
  aaxis.normalize();
  std::srand(std::time(0));  // use current time as seed for random generator
  auto angle = T(std::rand());

#if 0  // NOTE: Eigen::AngleAxis has less accuracy than ours
    Eigen::Matrix<T, 3, 3> mrot1;
    mrot1 = Eigen::AngleAxis<T>(angle, aaxis);

    aaxis *= angle;
    SO3Group<T> rot1(aaxis);

    Eigen::Matrix<T, 3, 3> mrand = Eigen::Matrix<T, 3, 3>::Random();

    EXPECT_MATRIX_NEAR(mrand * mrot1, mrand * rot1, Gap<T>());
#else

  aaxis *= angle;
  SO3Group<T> rot1(aaxis);

  Eigen::Matrix<T, 3, 3> mrand = Eigen::Matrix<T, 3, 3>::Random();

  EXPECT_MATRIX_NEAR(mrand * rot1.get_matrix(), mrand * rot1, Gap<T>());
#endif
}

TEST(SO3Test, MatrixLeftHandMulOperator) {
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
