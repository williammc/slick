#include <iostream>
#include <iomanip>

#include <Eigen/Eigen>

#include "slick/math/so3.h"

using namespace std;
using namespace Eigen;
using namespace slick;

void test_constructors(){
  cout << "SO3():\n" << SO3() << endl;
  Matrix<SlickScalar,3,1> v3(1.7,2.12,0.0);
  cout << "v3:" << v3 << endl;
  cout << "SO3(v3):\n" << SO3(v3) << endl;
  cout << "SO3::exp(v3):\n" << SO3::exp(v3).get_matrix().cast<float>() << endl;
  Matrix<SlickScalar,3,1> v3Minus = -v3;
  cout << "SO3(-v3):\n" << SO3(v3Minus) << endl;
  cout << "SO3(v3)*SO3(-v3):\n" << SO3(v3)*SO3(v3Minus) << endl;
  Matrix<SlickScalar,3,3> m3;
  m3.setIdentity();
  SO3 so3 = m3;
  cout << "SO3 SO3=m2:\n" << so3 << endl;
  cout << "SO3(v3)*m2:\n" << SO3(v3)*m3.block<3,3>(0,0) << endl;
  cout << "m2*SO3(v3):\n" << m3*SO3(v3) << endl;
  cout << "SO3(v3).inverse():\n" << SO3(v3).inverse() << endl;
  cout << "SO3(v3).ln():\n" << SO3(v3).ln() << endl;
  SO3 so3Test(v3);
  cout << "so3Test:" << so3Test << endl;
  cout << "so3Test.ln():" << so3Test.ln() << endl;
  Matrix<float,3,1> v3f=v3.cast<float>();
  cout << "SO3(v3)*SO3(-v3):" << SO3(v3)*SO3(v3Minus) << endl;
}

void test_operators(){
  Matrix<SlickScalar, 3, 1> v3(1.7, 2.12, 0.0);
  cout << "v3:" << v3 << endl;
  SO3 so3_1(v3);
  cout << "so3_1 = SO3(v3):\n" << so3_1 << endl;
  SO3 so3_2;
  so3_2=so3_1;
  cout << "so3_2 = so3_1:\n" << so3_2 << endl;
}

int main(int , char ** ){
  cout << "testing constructors ...\n" << endl;
  test_constructors();
  cout << "testing operators ...\n" << endl;
  test_operators();
  return 0;
}

