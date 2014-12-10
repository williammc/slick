#include <iostream>
#include <iomanip>

#include <Eigen/Eigen>

#include "slick/math/so3.h"
#include "slick/math/se3.h"

using namespace std;
using namespace Eigen;
using namespace slick;

void testConstructors(){
  cout << "SE3():\n" << SE3() << endl;
  SO3 so3;
  Matrix<SlickScalar,3,1> v3=Matrix<SlickScalar,3,1>::Zero();
  v3[0]=2;
  cout << "SE3(so3,v3):\n" << SE3(so3,v3) << endl;
  Matrix<SlickScalar,6,1> v6 = Matrix<SlickScalar,6,1>::Zero();
  v6[0]=7;
  v6[3]=2.11;
  v6[4]=4.21;
  cout << "SE3(v6):\n" << SE3(v6) << endl;
}

void testOperators(){
  SO3 so3;
  Matrix<SlickScalar,6,1> v6 = Matrix<SlickScalar,6,1>::Zero();
  v6[0]=7;
  v6[3]=2.11;
  v6[4]=4.21;
  cout << "v6:" << v6.transpose() << endl;
  Matrix<SlickScalar,6,1> v6Minus = -v6;
  cout << "SE3(v6)*SE3(-v6):\n" << SE3(v6)*SE3(v6Minus) << endl;
  Matrix<SlickScalar,4,4> m4 = Matrix<SlickScalar,4,4>::Identity();
  cout << "SE3(v6)*m4:\n" << SE3(v6)*m4 << endl;
  cout << "m4*SE3(v6):\n" << endl;
  cout << m4*SE3(v6) << endl;
  cout << "SE3(v6).inverse():\n" << SE3(v6).inverse() << endl;
  cout << "SE3(v6).ln():\n" << SE3(v6).ln() << endl;
  SE3 se3_1(v6);
  cout << "se3_1:\n" << se3_1 << endl;
  SE3 se3_2;
  se3_2 = se3_1;
  cout << "se3_2=se3_1:\n" << se3_2 << endl;
}

int main(int , char ** ){
  cout << "testing constructors ...\n" << endl;
  testConstructors();
  cout << "testing operators ...\n" << endl;
  testOperators();
  return 0;
}


