#include <gtest/gtest.h>
#include <Eigen/Geometry>
#include "slick/geometry/line2d.h"
#include "slick/test/unittest/util.h"

const int N = 10;

template <typename T>
void check_on_line(const slick::Line2DBase<T> &line,
                   const Eigen::Matrix<T, 2, 1> &pt) {
  const auto lv = (line.point2() - line.point1()).normalized();
  const T t = (pt - line.point1()).normalized().dot(lv);
  EXPECT_NEAR(std::fabs(t), 1.0, slick::Gap<T>());
};

// line2d specific functions ===================================================
template <typename T> void Project_Test() {
  std::srand(
      std::time(nullptr)); // use current time as seed for random generator
  slick::Line2DBase<T> line(slick::GenRandPoint2D<T>(N),
                            slick::GenRandPoint2D<T>(N));
  const auto pr = line.project(slick::GenRandPoint2D<T>(N));
  const auto lv = (line.point2() - line.point1()).normalized();
  const T t = (pr - line.point1()).normalized().dot(lv);
  EXPECT_NEAR(std::fabs(t), 1.0, slick::Gap<T>());
}

TEST(Line2DBaseTest, Project) {
  Project_Test<double>();
  Project_Test<float>();
}

template <typename T> void perpendicular_distance_Test(T gap) {
  std::srand(
      std::time(nullptr)); // use current time as seed for random generator
  slick::Line2DBase<T> line(slick::GenRandPoint2D<T>(N),
                            slick::GenRandPoint2D<T>(N));
  const auto pt = slick::GenRandPoint2D<T>(N);
  const auto pr = line.project(pt);
  const T t = line.perpendicular_distance(pt);
  EXPECT_NEAR(t, (pt - pr).norm(), gap);
}

TEST(Line2DBaseTest, perpendicular_distance) {
  perpendicular_distance_Test<double>(slick::Gap<double>());
  perpendicular_distance_Test<float>(1.e-3);
}

template <typename T> void perpendicular_sqdistance_Test(T gap) {
  std::srand(
      std::time(nullptr)); // use current time as seed for random generator
  slick::Line2DBase<T> line(slick::GenRandPoint2D<T>(N),
                            slick::GenRandPoint2D<T>(N));
  const auto pt = slick::GenRandPoint2D<T>(N);
  const auto pr = line.project(pt);
  const T t = line.perpendicular_sqdistance(pt);
  EXPECT_NEAR(t, (pt - pr).squaredNorm(), gap);
}

TEST(Line2DBaseTest, perpendicular_sqdistance) {
  perpendicular_sqdistance_Test<double>(slick::Gap<double>());
  // perpendicular_sqdistance_Test<float>(1.e-3);
}

template <typename T> void intersect_line_Test() {
  std::srand(
      std::time(nullptr)); // use current time as seed for random generator
  slick::Line2DBase<T> line1(slick::GenRandPoint2D<T>(N),
                               slick::GenRandPoint2D<T>(N));
  slick::Line2DBase<T> line2(slick::GenRandPoint2D<T>(N),
                               slick::GenRandPoint2D<T>(N));
  auto pt = line1.intersect(line2);
  if (pt) {
    check_on_line(line1, *pt);
    check_on_line(line1, *pt);
    check_on_line(line2, *pt);
    check_on_line(line2, *pt);
  }
}

TEST(Line2DBaseTest, intersect_line) {
  intersect_line_Test<double>();
  intersect_line_Test<float>();
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
