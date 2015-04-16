#pragma once
#include <Eigen/Core>
#include <Eigen/StdVector>
#include <vector>

namespace slick {
using SlickScalar = double;  // make it easy to change for entire lib.

// single precision types
using ColMat2Xf = Eigen::Matrix<float, 2, Eigen::Dynamic, Eigen::ColMajor>;
using ColMat4Xf = Eigen::Matrix<float, 4, Eigen::Dynamic, Eigen::ColMajor>;
using ColMat6Xf = Eigen::Matrix<float, 6, Eigen::Dynamic, Eigen::ColMajor>;

using Vec2f = Eigen::Matrix<float, 2, 1>;
using Vec3f = Eigen::Matrix<float, 3, 1>;
using Vec4f = Eigen::Matrix<float, 4, 1>;
using Vec6f = Eigen::Matrix<float, 6, 1>;

using Mat2f = Eigen::Matrix<float, 2, 2>;
using Mat3f = Eigen::Matrix<float, 3, 3>;
using Mat4f = Eigen::Matrix<float, 4, 4>;
using Mat6f = Eigen::Matrix<float, 6, 6>;

using Matrix34f = Eigen::Matrix<float, 3, 4>;

using AlignVec2fs = std::vector<Vec2f, Eigen::aligned_allocator<Vec2f>>;
using AlignVec4fs = std::vector<Vec4f, Eigen::aligned_allocator<Vec4f>>;
using AlignVec6fs = std::vector<Vec6f, Eigen::aligned_allocator<Vec6f>>;

// double precision types
using Vec2d = Eigen::Matrix<double, 2, 1>;
using Vec3d = Eigen::Matrix<double, 3, 1>;
using Vec4d = Eigen::Matrix<double, 4, 1>;
using Vec6d = Eigen::Matrix<double, 6, 1>;

using Mat2d = Eigen::Matrix<double, 2, 2>;
using Mat3d = Eigen::Matrix<double, 3, 3>;
using Mat4d = Eigen::Matrix<double, 4, 4>;
using Mat6d = Eigen::Matrix<double, 6, 6>;

// default precision types (this can affect entire lib(s))
using Vec2 = Eigen::Matrix<SlickScalar, 2, 1>;
using Vec3 = Eigen::Matrix<SlickScalar, 3, 1>;
using Vec4 = Eigen::Matrix<SlickScalar, 4, 1>;
using Vec6 = Eigen::Matrix<SlickScalar, 6, 1>;

using Mat2 = Eigen::Matrix<SlickScalar, 2, 2>;
using Mat2x3 = Eigen::Matrix<SlickScalar, 2, 3>;
using Mat3 = Eigen::Matrix<SlickScalar, 3, 3>;
using Mat4 = Eigen::Matrix<SlickScalar, 4, 4>;
using Mat3x4 = Eigen::Matrix<SlickScalar, 3, 4>;
using Mat6 = Eigen::Matrix<SlickScalar, 6, 6>;

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

EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(std::pair<int, slick::Vec2d>)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(std::pair<slick::Vec4d, slick::Vec2d>)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(std::pair<slick::Vec3d, slick::Vec3d>)