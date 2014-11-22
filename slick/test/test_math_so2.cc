#include <iostream>
#include <iomanip>

#include <Eigen/Core>

#include <math/so2.h>

using namespace std;
using namespace Eigen;
using namespace look3d;

void test_constructors(){
  cout << "SO2():\n" << SO2() << endl;
  cout << "SO2(1.0):\n" << SO2(1.0) << endl;
  cout << "SO2(1.0)*SO2(-1.0):\n" << SO2(1.0)*SO2(-1.0) << endl;
}

void test_operators(){
  Matrix<DefaultScalarType, 2, 2> m2;
  m2.setIdentity();
  SO2 rot = m2;
  cout << "SO2 rot = m2:\n" << rot << endl;
  rot = m2.block<2, 2>(0, 0);
  cout << "SO2 rot = m2.block<2, 2>(0, 0):\n" << rot << endl;
  SO2 rot1 = rot;
  cout << "SO2 rot1 = rot:\n" << rot1 << endl;
#if _MSC_VER >= 1800
  SO2 rot2 {1, 0, 0};
  cout << "SO2 rot2 {1, 0, 0, 1}:\n" << rot2 << endl;
#endif
  rot *= SO2(2.0);
  cout << "rot *= SO2(2.0):\n" << rot << endl;
  cout << "SO2(1.0)*m2:\n" << SO2(1.0)*m2.block<2, 2>(0, 0) << endl;
  cout << "m2*SO2(1.0):\n" << m2*SO2(1.0) << endl;
  cout << "SO2(1.0).inverse():\n" << SO2(1.0).inverse() << endl;
  cout << "SO2(1.0).ln():\n" << SO2(1.0).ln() << endl;
}

int main(int , char ** ){
  cout << "testing constructors ...\n" << endl;
  test_constructors();
  cout << "testing operators ...\n" << endl;
  test_operators();
  return 0;
}
