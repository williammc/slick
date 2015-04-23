#include <gtest/gtest.h>

#include "slick/geometry/line3d.h"
#include "slick/test/unittest/util.h"

const int N = 10;

// line3d specific functions ===================================================
template <typename T> void Project_Test() {
  std::srand(
      std::time(nullptr)); // use current time as seed for random generator
  slick::Line3DBase<T> line(slick::GenRandPoint<T>(N),
                            slick::GenRandPoint<T>(N));
  const auto pr = line.project(slick::GenRandPoint<T>(N));
  const auto lv = (line.point2() - line.point1()).normalized();
  const T t = (pr - line.point1()).normalized().dot(lv);
  EXPECT_NEAR(std::fabs(t), 1.0, slick::Gap<T>());
}

TEST(Line3DBaseTest, Project) {
  Project_Test<double>();
  Project_Test<float>();
}

template <typename T> void perpendicular_distance_Test() {
  std::srand(
      std::time(nullptr)); // use current time as seed for random generator
  slick::Line3DBase<T> line(slick::GenRandPoint<T>(N),
                            slick::GenRandPoint<T>(N));
  const auto pt = slick::GenRandPoint<T>(N);
  const auto pr = line.project(pt);
  const T t = line.perpendicular_distance(pt);
  EXPECT_NEAR(t, (pt - pr).norm(), slick::Gap<T>());
}

TEST(Line3DBaseTest, perpendicular_distance) {
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
