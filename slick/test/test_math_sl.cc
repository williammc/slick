
#include <iostream>
#include <iomanip>

#include <math/sl.h>

using namespace std;
using namespace Eigen;
using namespace look3d;
#if 0
void test_sl(){
  Matrix<DefaultScalarType,8,1> m_test;
  m_test << 1,0,-1,0,0,0,0,0;
  SL<DefaultScalarType,3> h(m_test);
  std::cout << h << std::endl;
  std::cout << h.inverse() << std::endl;
  std::cout << SL<DefaultScalarType,3>::exp((Eigen::Matrix<DefaultScalarType,8,1>() << -1,0,1,0,0,0,0,0).finished()) << std::endl;
  std::cout << h * h.inverse() << std::endl;
  h *= h.inverse();
  cout << h << endl;

  for(int i = 0; i <  SL<DefaultScalarType,3>::dim; ++i)
    std::cout << "generator " << i << "\n" <<  SL<DefaultScalarType,3>::generator(i) << std::endl;

  for(int i = 0; i <  SL<DefaultScalarType,2>::dim; ++i)
    std::cout << "generator " << i << "\n" <<  SL<DefaultScalarType,2>::generator(i) << std::endl;

  std::cout <<  SL<DefaultScalarType,2>::exp(Eigen::Matrix<DefaultScalarType,3,1>(1,2,3)) << std::endl;

  h =  SL<DefaultScalarType,3>::exp( (Eigen::Matrix<DefaultScalarType,8,1>() << 1,0,-1,0,0,0,1,0).finished());

  std::cout << h << "\n";
  Eigen::Matrix<DefaultScalarType,3,1> t(0,1,2);
  std::cout << "with vector " << t << "\n";
  std::cout << h * t << "\n";
  std::cout << t.transpose() * h << "\n";

  Eigen::Matrix<DefaultScalarType,3,5> m = Eigen::Matrix<DefaultScalarType,3,5>::Zero();
  m.row(0) = (Eigen::Matrix<DefaultScalarType, 5, 1>() << 0, 1, 2, 3, 4).finished();
  m.row(1) = (Eigen::Matrix<DefaultScalarType, 5, 1>() << 1, 2, 3, 4, -5).finished();
  m.row(2) = (Eigen::Matrix<DefaultScalarType, 5, 1>() << 2, 3, 4, 5, 8).finished();

  std::cout << "with matrix " << m << "\n";
  std::cout << h * m << "\n";
  std::cout << m.transpose() * h << "\n";

  std::cout << std::endl;
}
#endif

void test_operators(){
  return;
}

int main(int , char ** ){
  cout << "testing sl ...\n" << endl;
  //test_sl();
  cout << "testing operators ...\n" << endl;
  test_operators();
  return 0;
}

