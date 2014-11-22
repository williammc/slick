// Copyright 2014 The Slick Authors. All rights reserved.
#pragma once
#include <vector>
#include <cassert>
#include <algorithm>

#include "slick/slick_api.h"

namespace slick {

template <typename Precision>
struct SLICK_API MEstimator {
  MEstimator() {}
  virtual ~MEstimator() {}
  typedef Precision value_type;
  virtual Precision compute_sigma_squared(std::vector<Precision> values) = 0;
  virtual Precision weight(Precision error_squared) const = 0;
  virtual int get_outliers() const = 0;
  virtual Precision get_sigma_squared() const = 0;
};

template <typename Precision>
struct SLICK_API LeastSquares : public MEstimator<Precision> {
  LeastSquares();
  Precision sigmaSquared;

  Precision compute_sigma_squared(std::vector<Precision> values);
  Precision weight(Precision error_squared) const { return 1; }
  int get_outliers() const { return 0; }
  Precision get_sigma_squared() const { return sigmaSquared; }
};

template <typename Precision>
struct SLICK_API Tukey : public MEstimator<Precision> {
  Precision sigmaSquared;
  mutable int outliers;

  Tukey();

  Precision compute_sigma_squared(std::vector<Precision> values);
  Precision weight(Precision error_squared) const;
  int get_outliers() const { return outliers; }
  Precision get_sigma_squared() const { return sigmaSquared; }
};

template <typename Precision>
struct SLICK_API Huber : public MEstimator<Precision> {
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
//template struct SLICK_API MEstimator<T>;\
//template struct SLICK_API LeastSquares<T>;\
//template struct SLICK_API Tukey<T>;\
//template struct SLICK_API Huber<T>;
////template SLICK_API T LeastSquares<T>::compute_sigma_squared(std::vector<T>);\
////template SLICK_API Tukey<T>::Tukey();\
////template SLICK_API T Tukey<T>::compute_sigma_squared(std::vector<T>);\
////template SLICK_API T Tukey<T>::weight(T) const;\
////template SLICK_API Huber<T>::Huber();\
////template SLICK_API T Huber<T>::compute_sigma_squared(std::vector<T>);\
////template SLICK_API T Huber<T>::weight(T) const;
//
//INSTANTIATE_MESTIMATORS(float)
//INSTANTIATE_MESTIMATORS(double)
}  // namespace slick
