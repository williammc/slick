#include "slick/scene/triangulation.h"
#include <iostream>
#include <Eigen/Core>

#include "slick/math/se3.h"

int test_triangulate() {
  std::cout << "generting 3D points..." << std::endl;
  std::vector<Eigen::Matrix<slick::SlickScalar, 3, 1>> v3_points(100);
  for (size_t i = 0; i < v3_points.size(); ++i) {
    v3_points[i].setRandom();
    v3_points[i] *= 20.0; // scale up to value range to [-100,100]
    if (v3_points[i][2] < 0)
      v3_points[i][2] = -v3_points[i][2];
    v3_points[i][2] = 17.0; // all 3D points are on the same plane
  }
  // generates a  3degree rotation around axis-y & 1unit translate along axis-a
  slick::SE3 pose1;
  Eigen::Matrix<slick::SlickScalar, 3, 1> v3(0, 3.0 * M_PI / 180.0, 0);
  slick::SO3Group<slick::SlickScalar> so3_r(v3);
  slick::SE3Group<slick::SlickScalar> pose2(
      so3_r, Eigen::Matrix<slick::SlickScalar, 3, 1>(-1.0, 0.0, 0.0));

  std::vector<slick::SlickScalar> norms;
  for (int i = 0; i < v3_points.size(); ++i) {
    Eigen::Matrix<slick::SlickScalar, 2, 1> v2_p1 = slick::project(pose1 * v3_points[i]);
    Eigen::Matrix<slick::SlickScalar, 2, 1> v2_p2 = slick::project(pose2 * v3_points[i]);
    std::pair<Eigen::Matrix<slick::SlickScalar, 3, 1>, slick::SlickScalar> res =
        slick::triangulate(v2_p1, v2_p2, pose1, pose2);

    std::cout << "norm(real-estimate):" << (v3_points[i] - res.first).norm()
              << std::endl;
    norms.push_back((v3_points[i] - res.first).norm());
  }

  slick::SlickScalar sum_norm(0.0);
  for (std::vector<slick::SlickScalar>::iterator j = norms.begin();
       j != norms.end(); ++j)
    sum_norm += *j;
  std::cout << "sum_norm:"
            << "mean_norm:" << sum_norm / norms.size() << std::endl;
  if (sum_norm / norms.size() > 1.e-9)
    return -1;
  return 1;
}

int main(int argc, char **argv) {
  bool ok = true, final_ok = true;
  ok = test_triangulate();
  std::cout << "test_triangulate():" << ok << std::endl;
  final_ok &= ok;

  return (final_ok) ? 0 : -1;
}
