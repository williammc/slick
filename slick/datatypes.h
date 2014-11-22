#pragma once
#include <Eigen/Core>
#include <Eigen/StdVector>
#include <vector>

namespace slick {
typedef double SlickScalar;  // make it easy to change for entire lib.

// single precision types
typedef Eigen::Matrix<float, 2, Eigen::Dynamic, Eigen::ColMajor> ColMat2Xf;
typedef Eigen::Matrix<float, 4, Eigen::Dynamic, Eigen::ColMajor> ColMat4Xf;
typedef Eigen::Matrix<float, 6, Eigen::Dynamic, Eigen::ColMajor> ColMat6Xf;

typedef Eigen::Matrix<float, 2, 1> Vec2f;
typedef Eigen::Matrix<float, 3, 1> Vec3f;
typedef Eigen::Matrix<float, 4, 1> Vec4f;
typedef Eigen::Matrix<float, 6, 1> Vec6f;

typedef Eigen::Matrix<float, 2, 2> Mat2f;
typedef Eigen::Matrix<float, 3, 3> Mat3f;
typedef Eigen::Matrix<float, 4, 4> Mat4f;
typedef Eigen::Matrix<float, 6, 6> Mat6f;

typedef Eigen::Matrix<float, 3, 4> Matrix34f;

typedef std::vector<Vec2f, Eigen::aligned_allocator<Vec2f>> AlignVec2fs;
typedef std::vector<Vec4f, Eigen::aligned_allocator<Vec4f>> AlignVec4fs;
typedef std::vector<Vec6f, Eigen::aligned_allocator<Vec6f>> AlignVec6fs;

// double precision types
typedef Eigen::Matrix<double, 2, 1> Vec2d;
typedef Eigen::Matrix<double, 3, 1> Vec3d;
typedef Eigen::Matrix<double, 4, 1> Vec4d;
typedef Eigen::Matrix<double, 6, 1> Vec6d;

typedef Eigen::Matrix<double, 2, 2> Mat2d;
typedef Eigen::Matrix<double, 3, 3> Mat3d;
typedef Eigen::Matrix<double, 4, 4> Mat4d;
typedef Eigen::Matrix<double, 6, 6> Mat6d;

// default precision types (this can affect entire lib(s))
typedef Eigen::Matrix<SlickScalar, 2, 1> Vec2;
typedef Eigen::Matrix<SlickScalar, 3, 1> Vec3;
typedef Eigen::Matrix<SlickScalar, 4, 1> Vec4;
typedef Eigen::Matrix<SlickScalar, 6, 1> Vec6;

typedef Eigen::Matrix<SlickScalar, 2, 2> Mat2;
typedef Eigen::Matrix<SlickScalar, 2, 3> Mat2x3;
typedef Eigen::Matrix<SlickScalar, 3, 3> Mat3;
typedef Eigen::Matrix<SlickScalar, 4, 4> Mat4;
typedef Eigen::Matrix<SlickScalar, 3, 4> Mat3x4;
typedef Eigen::Matrix<SlickScalar, 6, 6> Mat6;

}  // namespace slick

// Defines common vectors with Eigen elements

EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(slick::Vec2f)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(slick::Vec3f)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(slick::Vec4f)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(slick::Mat2f)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(slick::Mat3f)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(slick::Mat4f)

EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(slick::Vec2d)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(slick::Vec3d)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(slick::Vec4d)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(slick::Mat2d)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(slick::Mat3d)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(slick::Mat4d)

EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(std::pair<int, slick::Vec2>)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(std::pair<slick::Vec4, slick::Vec2>)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(std::pair<slick::Vec3, slick::Vec3>)
