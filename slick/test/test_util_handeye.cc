#include <iostream>
#include <fstream>
#ifndef WIN32
#include <tr1/tuple>
#else
#include <tuple>
#endif
#include "geometry/handeye.h"

#if defined(HAVE_TOON) && defined(HAVE_TAG)
#include <TooN/TooN.h>
#include <tag/handeye.h>
#endif

#include "config.h"

using namespace slick;
using namespace std;

int test_computeHandEye() {
  stringstream ss;
  ss << DATA_DIR << "/" << "handeye-1.txt";
  fstream fs(ss.str().c_str(), std::ios_base::in);
  vector<SO3 > sensorworld2sensor_rots, cam2vision_rots;
  Eigen::Matrix<DefaultScalarType, 3, 1> temp_v;
  for (int i=0; i<6; ++i) {
    fs >> temp_v[0] >> temp_v[1] >> temp_v[2];
    SO3 sw2s_rot(temp_v);

    fs >> temp_v[0] >> temp_v[1] >> temp_v[2];
    SO3 c2v_rot(temp_v);
    sensorworld2sensor_rots.push_back(sw2s_rot);
    cam2vision_rots.push_back(c2v_rot);
  }

  std::pair<SO3 , SO3 > res  = ComputeHandEye( sensorworld2sensor_rots, cam2vision_rots );

  cout << "sensor2cam:\n" << res.first << endl;
  cout << "vision2sensorworld:\n" << res.second << endl;

  return 1;
}


int test_computeHandEye_TooN() {
  
#if defined(HAVE_TOON) && defined(HAVE_TAG)
  stringstream ss;
  ss << DATA_DIR << "/" << "handeye-1.txt";
  fstream fs(ss.str().c_str(), std::ios_base::in);
  vector<TooN::SO3 > sensorworld2sensor_rots, cam2vision_rots;
  TooN::Vector<3> temp_v;
  for (int i=0; i<6; ++i) {
    fs >> temp_v[0] >> temp_v[1] >> temp_v[2];
    TooN::SO3 sw2s_rot(temp_v);

    fs >> temp_v[0] >> temp_v[1] >> temp_v[2];
    TooN::SO3 c2v_rot(temp_v);
    sensorworld2sensor_rots.push_back(sw2s_rot);
    cam2vision_rots.push_back(c2v_rot);
  }

  std::pair<TooN::SO3 , TooN::SO3 > res  = tag::computeHandEye( sensorworld2sensor_rots, cam2vision_rots );

  cout << "sensor2cam:\n" << res.first << endl;
  cout << "vision2sensorworld:\n" << res.second << endl;
#endif
  return 1;
}

int main (int argc, char** argv) {
  bool final_ok=true;
  bool ok = test_computeHandEye();
  cout << "test_computeHandEye():" << ok  << endl;
  test_computeHandEye_TooN();
  final_ok &= ok;
  return (final_ok)?1:-1;
}
