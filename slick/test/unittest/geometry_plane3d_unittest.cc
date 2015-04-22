#include <gtest/gtest.h>

#include "slick/geometry/plane3d.h"
#include "slick/test/unittest/util.h"

template <typename T> inline Eigen::Matrix<T, 3, 1> GenRandPoint() {
  const T N = 10;
  const T x = T(std::rand()) / T(RAND_MAX) * N;
  const T y = T(std::rand()) / T(RAND_MAX) * N;
  const T z = T(std::rand()) / T(RAND_MAX) * N;
  return Eigen::Matrix<T, 3, 1>(x, y, z);
}

template <typename T>
void check_on_plane(const slick::Plane3DBase<T> &pln,
                    const Eigen::Matrix<T, 3, 1> &pt) {

  const T t = pt.dot(pln.equation().head<3>()) + pln.equation()[3];
  EXPECT_NEAR(std::fabs(t), 0.0, slick::Gap<T>());
};

// plane3d specific functions ==================================================
template <typename T> void Project_Test() {
  std::srand(
      std::time(nullptr)); // use current time as seed for random generator
  slick::Plane3DBase<T> plane(GenRandPoint<T>(), GenRandPoint<T>(),
                              GenRandPoint<T>());
  const auto pr = plane.project(GenRandPoint<T>());
  const auto eq = plane.equation();
  const T t = pr[0] * eq[0] + pr[1] * eq[1] + pr[2] * eq[2] + eq[3];
  EXPECT_NEAR(std::fabs(t), 0.0, slick::Gap<T>());
}

TEST(Plane3DBaseTest, Project) {
  Project_Test<double>();
  Project_Test<float>();
}

template <typename T> void perpendicular_distance_Test() {
  std::srand(
      std::time(nullptr)); // use current time as seed for random generator
  slick::Plane3DBase<T> plane(GenRandPoint<T>(), GenRandPoint<T>(),
                              GenRandPoint<T>());
  const auto pt = GenRandPoint<T>();
  const auto pr = plane.project(pt);
  const T t = plane.perpendicular_distance(pt);
  EXPECT_NEAR(t, (pt - pr).norm(), slick::Gap<T>());
}

TEST(Plane3DBaseTest, perpendicular_distance) {
  perpendicular_distance_Test<double>();
  perpendicular_distance_Test<float>();
}

template <typename T> void intersect_line_Test() {
  std::srand(
      std::time(nullptr)); // use current time as seed for random generator
  slick::Plane3DBase<T> plane(GenRandPoint<T>(), GenRandPoint<T>(),
                              GenRandPoint<T>());
  slick::Line3DBase<T> line(GenRandPoint<T>(), GenRandPoint<T>());
  auto res = plane.intersect(line);
  if (res.first) {
    check_on_plane(plane, res.second);
  }
}

TEST(Plane3DBaseTest, intersect_line) {
  intersect_line_Test<double>();
  intersect_line_Test<float>();
}

template <typename T> void intersect_plane_Test() {
  std::srand(
      std::time(nullptr)); // use current time as seed for random generator
  slick::Plane3DBase<T> plane1(GenRandPoint<T>(), GenRandPoint<T>(),
                               GenRandPoint<T>());
  slick::Plane3DBase<T> plane2(GenRandPoint<T>(), GenRandPoint<T>(),
                               GenRandPoint<T>());
  auto res = plane1.intersect(plane2);
  if (res.first) {
    check_on_plane(plane1, res.second.point1());
    check_on_plane(plane1, res.second.point2());
    check_on_plane(plane2, res.second.point1());
    check_on_plane(plane2, res.second.point2());
  }
}

TEST(Plane3DBaseTest, intersect_plane) {
  intersect_plane_Test<double>();
  intersect_plane_Test<float>();
}

int main(int argc, char **argv) {
  std::vector<char *> vars(argc + 1);
  for (int i = 0; i < argc; ++i)
    vars[i] = argv[i];
  char ca[50];
  sprintf(ca, "--gtest_repeat=%d", NSAMPLES);
  vars[argc] = ca;
  argc++;
  ::testing::InitGoogleTest(&argc, &vars.front());
  return RUN_ALL_TESTS();
}
