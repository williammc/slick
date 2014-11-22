// Copyright 2012, 2013, 2014 The Look3D Authors. All rights reserved.
#ifndef LOOK3D_MATH_MESTIMATOR_H_
#define LOOK3D_MATH_MESTIMATOR_H_
#include <vector>
#include <cassert>
#include <algorithm>

#include "math/look3d_math_api.h"

namespace look3d {

template <typename Precision>
struct LOOK3D_MATH_API MEstimator {
  MEstimator() {}
  virtual ~MEstimator() {}
  typedef Precision value_type;
  virtual Precision compute_sigma_squared(std::vector<Precision> values) = 0;
  virtual Precision weight(Precision error_squared) const = 0;
  virtual int get_outliers() const = 0;
  virtual Precision get_sigma_squared() const = 0;
};

template <typename Precision>
struct LOOK3D_MATH_API LeastSquares : public MEstimator<Precision> {
  LeastSquares();
  Precision sigmaSquared;

  Precision compute_sigma_squared(std::vector<Precision> values);
  Precision weight(Precision error_squared) const { return 1; }
  int get_outliers() const { return 0; }
  Precision get_sigma_squared() const { return sigmaSquared; }
};

template <typename Precision>
struct LOOK3D_MATH_API Tukey : public MEstimator<Precision> {
  Precision sigmaSquared;
  mutable int outliers;

  Tukey();

  Precision compute_sigma_squared(std::vector<Precision> values);
  Precision weight(Precision error_squared) const;
  int get_outliers() const { return outliers; }
  Precision get_sigma_squared() const { return sigmaSquared; }
};

template <typename Precision>
struct LOOK3D_MATH_API Huber : public MEstimator<Precision> {
  Precision sigmaSquared;
  mutable int outliers;

  Huber();

  Precision compute_sigma_squared(std::vector<Precision> values);
  Precision weight(Precision error_squared) const;
  int get_outliers() const { return outliers; }
  Precision get_sigma_squared() const { return sigmaSquared; }
};

// Instantiate template functions ==============================================
//#define INSTANTIATE_MESTIMATORS(T)                                             \
//template struct LOOK3D_MATH_API MEstimator<T>;\
//template struct LOOK3D_MATH_API LeastSquares<T>;\
//template struct LOOK3D_MATH_API Tukey<T>;\
//template struct LOOK3D_MATH_API Huber<T>;
////template LOOK3D_MATH_API T LeastSquares<T>::compute_sigma_squared(std::vector<T>);\
////template LOOK3D_MATH_API Tukey<T>::Tukey();\
////template LOOK3D_MATH_API T Tukey<T>::compute_sigma_squared(std::vector<T>);\
////template LOOK3D_MATH_API T Tukey<T>::weight(T) const;\
////template LOOK3D_MATH_API Huber<T>::Huber();\
////template LOOK3D_MATH_API T Huber<T>::compute_sigma_squared(std::vector<T>);\
////template LOOK3D_MATH_API T Huber<T>::weight(T) const;
//
//INSTANTIATE_MESTIMATORS(float)
//INSTANTIATE_MESTIMATORS(double)
}  // namespace look3d
#endif  // LOOK3D_MATH_MESTIMATOR_H_
