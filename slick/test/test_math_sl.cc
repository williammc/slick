
#include <iostream>
#include <iomanip>

#include "slick/math/sl.h"

using namespace std;
using namespace Eigen;
using namespace slick;
#if 0
void test_sl(){
  Matrix<SlickScalar,8,1> m_test;
  m_test << 1,0,-1,0,0,0,0,0;
  SL<SlickScalar,3> h(m_test);
  std::cout << h << std::endl;
  std::cout << h.inverse() << std::endl;
  std::cout << SL<SlickScalar,3>::exp((Eigen::Matrix<SlickScalar,8,1>() << -1,0,1,0,0,0,0,0).finished()) << std::endl;
  std::cout << h * h.inverse() << std::endl;
  h *= h.inverse();
  cout << h << endl;

  for(int i = 0; i <  SL<SlickScalar,3>::dim; ++i)
    std::cout << "generator " << i << "\n" <<  SL<SlickScalar,3>::generator(i) << std::endl;

  for(int i = 0; i <  SL<SlickScalar,2>::dim; ++i)
    std::cout << "generator " << i << "\n" <<  SL<SlickScalar,2>::generator(i) << std::endl;

  std::cout <<  SL<SlickScalar,2>::exp(Eigen::Matrix<SlickScalar,3,1>(1,2,3)) << std::endl;

  h =  SL<SlickScalar,3>::exp( (Eigen::Matrix<SlickScalar,8,1>() << 1,0,-1,0,0,0,1,0).finished());

  std::cout << h << "\n";
  Eigen::Matrix<SlickScalar,3,1> t(0,1,2);
  std::cout << "with vector " << t << "\n";
  std::cout << h * t << "\n";
  std::cout << t.transpose() * h << "\n";

  Eigen::Matrix<SlickScalar,3,5> m = Eigen::Matrix<SlickScalar,3,5>::Zero();
  m.row(0) = (Eigen::Matrix<SlickScalar, 5, 1>() << 0, 1, 2, 3, 4).finished();
  m.row(1) = (Eigen::Matrix<SlickScalar, 5, 1>() << 1, 2, 3, 4, -5).finished();
  m.row(2) = (Eigen::Matrix<SlickScalar, 5, 1>() << 2, 3, 4, 5, 8).finished();

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

