#include "slick/scene/triangulation.h"
#include <iostream>
#include <Eigen/Core>

#include "slick/math/se3.h"

using namespace Eigen;
using namespace std;
using namespace slick;

int test_triangulate() {
  cout << "generting 3D points..." << endl;
  std::vector<Matrix<SlickScalar,3,1> > v3_points(100);
  for (size_t i = 0; i < v3_points.size(); ++i) {
    v3_points[i].setRandom();
    v3_points[i] *= 20.0;  // scale up to value range to [-100,100]
    if(v3_points[i][2] < 0)
      v3_points[i][2] =-v3_points[i][2];
    v3_points[i][2] = 17.0; // all 3D points are on the same plane
  }
  // generates a  3degree rotation around axis-y & 1unit translate along axis-a
  SE3 pose1;
  Matrix<SlickScalar,3,1> v3(0,3.0*M_PI/180.0,0);
  SO3Group<SlickScalar> so3_r(v3);
  SE3Group<SlickScalar> pose2(so3_r,Matrix<SlickScalar,3,1>(-1.0,0.0,0.0));

  std::vector<SlickScalar> norms;
  for (int i=0; i<v3_points.size(); ++i) {
    Matrix<SlickScalar,2,1> v2_p1 = project(pose1*v3_points[i]);
    Matrix<SlickScalar,2,1> v2_p2 = project(pose2*v3_points[i]);
    std::pair<Eigen::Matrix<SlickScalar,3,1>,SlickScalar> res = triangulate(v2_p1, v2_p2, pose1, pose2);

    cout << "norm(real-estimate):" << (v3_points[i]-res.first).norm() << endl;
    norms.push_back((v3_points[i]-res.first).norm());
  }

  SlickScalar sum_norm=0.0;
  for (std::vector<SlickScalar>::iterator j=norms.begin();j!=norms.end();++j)
      sum_norm += *j;
  cout << "sum_norm:" << sum_norm << endl;
  cout << "mean_norm:" << sum_norm/norms.size() << endl;
  if(sum_norm/norms.size() >1.e-9)
    return -1;
  return 1;
}

int main (int argc, char** argv) {
  bool ok=true, final_ok=true;
  ok = test_triangulate();
  cout << "test_triangulate():" << ok << endl;
  final_ok &= ok;

  return (final_ok)?1:-1;
}
