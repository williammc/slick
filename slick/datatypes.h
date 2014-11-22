#ifndef LOOK3D_MATH_COMMON_H_
#define LOOK3D_MATH_COMMON_H_
#include <Eigen/Core>
#include <Eigen/StdVector>
#include <vector>

namespace look3d {
typedef double DefaultScalarType;  // make it easy to change for entire lib.

#ifdef WIN32
typedef unsigned char uint8_t;
#endif

typedef Eigen::Matrix<DefaultScalarType, 2, 1> Vector2;
typedef Eigen::Matrix<DefaultScalarType, 3, 1> Vector3;
typedef Eigen::Matrix<DefaultScalarType, 4, 1> Vector4;
typedef Eigen::Matrix<DefaultScalarType, 6, 1> Vector6;

typedef Eigen::Matrix<DefaultScalarType, 2, 2> Matrix2;
typedef Eigen::Matrix<DefaultScalarType, 2, 3> Matrix2x3;
typedef Eigen::Matrix<DefaultScalarType, 3, 3> Matrix3;
typedef Eigen::Matrix<DefaultScalarType, 4, 4> Matrix4;
typedef Eigen::Matrix<DefaultScalarType, 3, 4> Matrix3x4;
typedef Eigen::Matrix<DefaultScalarType, 6, 6> Matrix6;

}  // namespace look3d

// Defines common vectors with Eigen elements
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(std::pair<int, look3d::Vector2>)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(std::pair<look3d::Vector4, look3d::Vector2>)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(std::pair<look3d::Vector3, look3d::Vector3>)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(look3d::Matrix2)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(look3d::Matrix3)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(look3d::Matrix4)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(look3d::Vector2)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(look3d::Vector3)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(look3d::Vector4)
//sEIGEN_DEFINE_STL_VECTOR_SPECIALIZATION()
#endif  // LOOK3D_MATH_COMMON_H_
