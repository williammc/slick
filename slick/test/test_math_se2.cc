#include <iostream>

#include <Eigen/Eigen>

#include <math/so2.h>
#include <math/se2.h>

using namespace std;
using namespace Eigen;
using namespace look3d;

void test_constructors(){
  cout << "SE2():\n" << SE2() << endl;
  SO2 so2;
  Matrix<DefaultScalarType,2,1> v2=Matrix<DefaultScalarType,2,1>::Zero();
  v2[0]=1;
  cout << "SE2(so2,v2):\n" << SE2(so2,v2) << endl;
  Matrix<DefaultScalarType,3,1> v3 = Matrix<DefaultScalarType,3,1>::Zero();
  v3[0]=1;
  cout << "SE2(v3):\n" << SE2(v3) << endl;
}

void test_operators(){
  SO2 so2;
  Matrix<DefaultScalarType,3,1> v3 = Matrix<DefaultScalarType,3,1>::Zero();
  v3[0]=1;
  Matrix<DefaultScalarType,3,1> v3Minus = -v3;
  cout << "SE2(v3)*SE2(-v3):\n" << SE2(v3)*SE2(v3Minus) << endl;
  Matrix<DefaultScalarType,3,3> m3 = Matrix<DefaultScalarType,3,3>::Identity();
  cout << "SE2(v3)*m3:\n" << SE2(v3)*m3 << endl;
  cout << "m3*SE2(v3):\n" << endl;
  cout << m3*SE2(v3) << endl;
  cout << "SE2(v3).inverse():\n" << SE2(v3).inverse() << endl;
  cout << "SE2(v3).ln():\n" << SE2(v3).ln() << endl;
  cout << "SO2()*SE2(v3):\n" << SO2()*SE2(v3) << endl;
}

int main(int , char ** )
{
  cout << "testing constructors ...\n" << endl;
  test_constructors();
  cout << "testing operators ...\n" << endl;
  test_operators();
  return 0;
}

