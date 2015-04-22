#include <gtest/gtest.h>

#include "slick/geometry/line2d.h"
#include "slick/test/unittest/util.h"

template <typename Scalar> inline Eigen::Matrix<Scalar, 2, 1> GenRandPoint() {
  const Scalar x = Scalar(std::rand()) / Scalar(RAND_MAX);
  const Scalar y = Scalar(std::rand()) / Scalar(RAND_MAX);
  return Eigen::Matrix<Scalar, 2, 1>(x, y);
}

// line2d specific functions ===================================================
template <typename T> void Project_Test() {
  std::srand(std::time(0)); // use current time as seed for random generator
  slick::Line2DBase<T> line(GenRandPoint<T>(), GenRandPoint<T>());
  const auto pr = line.project(GenRandPoint<T>());
  const auto lv = (line.point2() - line.point1()).normalized();
  const T t = (pr - line.point1()).normalized().dot(lv);
  EXPECT_NEAR(std::fabs(t), 1.0, slick::Gap<T>());
}

TEST(Line2DBaseTest, Project) {
  Project_Test<double>();
  Project_Test<float>();
}

template <typename T> void perpendicular_distance_Test() {
  std::srand(std::time(0)); // use current time as seed for random generator
  slick::Line2DBase<T> line(GenRandPoint<T>(), GenRandPoint<T>());
  const auto pt = GenRandPoint<T>();
  const auto pr = line.project(pt);
  const T t = line.perpendicular_distance(pt);
  EXPECT_NEAR(t, (pt - pr).norm(), slick::Gap<T>());
}

TEST(Line2DBaseTest, perpendicular_distance) {
  perpendicular_distance_Test<double>();
  perpendicular_distance_Test<float>();
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
