#include <gtest/gtest.h>

#include "slick/geometry/plane3d.h"
#include "slick/test/unittest/util.h"

const int N = 10;

template <typename T>
slick::Plane3DBase<T> GenerateAxisAlignedPlane(int index, T value) {
  auto p1 = slick::GenRandPoint<T>(N);
  p1[index] = value;
  auto p2 = slick::GenRandPoint<T>(N);
  p2[index] = value;
  auto p3 = slick::GenRandPoint<T>(N);
  p3[index] = value;
  return slick::Plane3DBase<T>(p1, p2, p3);
}

template <typename T>
void check_on_plane(const slick::Plane3DBase<T> &pln,
                    const Eigen::Matrix<T, 3, 1> &pt) {

  const T t = pt.dot(pln.equation().head(3)) + pln.equation()[3];
  EXPECT_NEAR(std::fabs(t), 0.0, slick::Gap<T>());
};

// plane3d specific functions ==================================================
template <typename T> void Project_Test() {
  std::srand(
      std::time(nullptr)); // use current time as seed for random generator
  slick::Plane3DBase<T> plane(slick::GenRandPoint<T>(N),
                              slick::GenRandPoint<T>(N),
                              slick::GenRandPoint<T>(N));
  const auto pr = plane.project(slick::GenRandPoint<T>(N));
  const auto eq = plane.equation();
  const T t = pr[0] * eq[0] + pr[1] * eq[1] + pr[2] * eq[2] + eq[3];
  EXPECT_NEAR(std::fabs(t), 0.0, slick::Gap<T>());
}

TEST(Plane3DBaseTest, Project) {
  Project_Test<double>();
  Project_Test<float>();
}

/// check for axis-aligned planes
template <typename T> void axisaligned_perpendicular_distance_Test() {
  for (int i = 0; i < 3; ++i) {
    const auto t = slick::GenRandNumber<T>(N);
    const auto pln = GenerateAxisAlignedPlane<T>(i, t);
    const auto pt = slick::GenRandPoint<T>(N);
    const T d = pln.perpendicular_distance(pt);
    EXPECT_NEAR(d, std::fabs(t - pt[i]), slick::Gap<T>());
  }
}

TEST(Plane3DBaseTest, axisaligned_perpendicular_distance) {
  axisaligned_perpendicular_distance_Test<double>();
  axisaligned_perpendicular_distance_Test<float>();
}

template <typename T> void perpendicular_distance_Test() {
  std::srand(
      std::time(nullptr)); // use current time as seed for random generator
  slick::Plane3DBase<T> plane(slick::GenRandPoint<T>(N),
                              slick::GenRandPoint<T>(N),
                              slick::GenRandPoint<T>(N));
  const auto pt = slick::GenRandPoint<T>(N);
  const auto pr = plane.project(pt);
  const T t = plane.perpendicular_distance(pt);
  EXPECT_NEAR(t, (pt - pr).norm(), slick::Gap<T>());
}

TEST(Plane3DBaseTest, perpendicular_distance) {
  perpendicular_distance_Test<double>();
  perpendicular_distance_Test<float>();
}


/// check for axis-aligned planes
template <typename T> void axisaligned_intersect_ray_Test() {
  for (int i = 0; i < 3; ++i) {
    const auto t = slick::GenRandNumber<T>(N);
    const auto pln = GenerateAxisAlignedPlane<T>(i, t);
    const auto ray = slick::GenRandPoint<T>(N);
    auto pt = pln.intersect(ray);
    if (pt) {
      EXPECT_NEAR(t, (*pt)[i], slick::Gap<T>());
      const auto d = ray.normalized().dot(pt->normalized());
      EXPECT_NEAR(std::fabs(d), 1.0, slick::Gap<T>());
    }
  }
}

TEST(Plane3DBaseTest, axisaligned_intersect_ray) {
  axisaligned_intersect_ray_Test<double>();
  axisaligned_intersect_ray_Test<float>();
}

template <typename T> void intersect_ray_Test() {
  std::srand(
    std::time(nullptr)); // use current time as seed for random generator
  slick::Plane3DBase<T> plane(slick::GenRandPoint<T>(N),
                              slick::GenRandPoint<T>(N),
                              slick::GenRandPoint<T>(N));
  const auto ray = slick::GenRandPoint<T>(N);
  auto pt = plane.intersect(ray);
  if (pt) {
    check_on_plane(plane, *pt);
    const auto d = ray.normalized().dot(pt->normalized());
    EXPECT_NEAR(std::fabs(d), 1.0, slick::Gap<T>());
  }
}

TEST(Plane3DBaseTest, intersect_ray) {
  intersect_ray_Test<double>();
  intersect_ray_Test<float>();
}

/// check for axis-aligned planes
template <typename T> void axisaligned_intersect_line_Test() {
  for (int i = 0; i < 3; ++i) {
    const auto t = slick::GenRandNumber<T>(N);
    const auto pln = GenerateAxisAlignedPlane<T>(i, t);
    slick::Line3DBase<T> line(slick::GenRandPoint<T>(N),
                              slick::GenRandPoint<T>(N));
    auto pt = pln.intersect(line);
    if (pt) {
      EXPECT_NEAR(t, (*pt)[i], slick::Gap<T>());
    }
  }
}


TEST(Plane3DBaseTest, axisaligned_intersect_line) {
  axisaligned_intersect_line_Test<double>();
  axisaligned_intersect_line_Test<float>();
}

template <typename T> void intersect_line_Test() {
  std::srand(
      std::time(nullptr)); // use current time as seed for random generator
  slick::Plane3DBase<T> plane(slick::GenRandPoint<T>(N),
                              slick::GenRandPoint<T>(N),
                              slick::GenRandPoint<T>(N));
  slick::Line3DBase<T> line(slick::GenRandPoint<T>(N),
                            slick::GenRandPoint<T>(N));
  auto pt = plane.intersect(line);
  if (pt) {
    check_on_plane(plane, *pt);
  }
}

TEST(Plane3DBaseTest, intersect_line) {
  intersect_line_Test<double>();
  intersect_line_Test<float>();
}

/// check for axis-aligned planes intersection
template <typename T> void axisaligned_intersect_plane_Test() {
  for (int i = 0; i < 3; ++i) {
    const auto t = slick::GenRandNumber<T>(N);
    const auto pln = GenerateAxisAlignedPlane<T>(i, t);
    const auto t1 = slick::GenRandNumber<T>(N);
    int i1 = (i + 1) % 3;
    const auto pln1 = GenerateAxisAlignedPlane<T>(i1, t1);
    auto line = pln.intersect(pln1);
    EXPECT_NEAR(t, (*line).point1()[i], slick::Gap<T>());
    EXPECT_NEAR(t, (*line).point2()[i], slick::Gap<T>());
    EXPECT_NEAR(t1, (*line).point1()[i1], slick::Gap<T>());
    EXPECT_NEAR(t1, (*line).point2()[i1], slick::Gap<T>());
  }
}

TEST(Plane3DBaseTest, axisaligned_intersect_plane) {
  axisaligned_intersect_plane_Test<double>();
  axisaligned_intersect_plane_Test<float>();
}

template <typename T> void intersect_plane_Test() {
  std::srand(
      std::time(nullptr)); // use current time as seed for random generator
  slick::Plane3DBase<T> plane1(slick::GenRandPoint<T>(N),
                               slick::GenRandPoint<T>(N),
                               slick::GenRandPoint<T>(N));
  slick::Plane3DBase<T> plane2(slick::GenRandPoint<T>(N),
                               slick::GenRandPoint<T>(N),
                               slick::GenRandPoint<T>(N));
  auto line = plane1.intersect(plane2);
  if (line) {
    check_on_plane(plane1, (*line).point1());
    check_on_plane(plane1, (*line).point2());
    check_on_plane(plane2, (*line).point1());
    check_on_plane(plane2, (*line).point2());
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
