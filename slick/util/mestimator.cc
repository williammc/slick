// Copyright 2012, 2013, 2014 The Look3D Authors. All rights reserved.
#include "math/mestimator.h"

namespace look3d {

template <typename Precision>
LeastSquares<Precision>::LeastSquares() : sigmaSquared(0) {}

template <typename Precision>
Precision LeastSquares<Precision>::compute_sigma_squared(std::vector<Precision> values) {
  sigmaSquared = 0;
  if (values.size() <= 0)
    return sigmaSquared;
  for (typename std::vector<Precision>::const_iterator i = values.begin();
      i != values.end(); ++i)
    sigmaSquared += *i;
  sigmaSquared /= values.size();
  return sigmaSquared;
}

template <typename Precision>
Tukey<Precision>::Tukey() : sigmaSquared(0), outliers(0) {}

template <typename Precision>
Precision Tukey<Precision>::compute_sigma_squared(std::vector<Precision> values) {
  if (values.size() == 0)
    return 0;
  std::sort(values.begin(), values.end());
  Precision median_squared = values[values.size()/2];
  if (median_squared == 0)
    median_squared = Precision(1.e-19);
  // from Georg Klein
  Precision sigma = Precision(1.4826) *
      (1 + Precision(5) / (values.size() * 2 - 6)) * sqrt(median_squared);
  sigma = Precision(4.6851) * sigma;
  sigmaSquared = sigma * sigma;
  return sigmaSquared;
}
template <typename Precision>
Precision Tukey<Precision>::weight(Precision error_squared) const {
  if (error_squared > sigmaSquared) {
    ++outliers;
    return 0;
  } else {
    const Precision w = 1 - error_squared / sigmaSquared;
    return w*w;
  }
}

template <typename Precision>
Huber<Precision>::Huber() : sigmaSquared(0), outliers(0) {}


template <typename Precision>
Precision Huber<Precision>::compute_sigma_squared(std::vector<Precision> values) {
  assert(values.size() > 0);
  std::sort(values.begin(), values.end());
  Precision median_squared = values[values.size()/2];
  // from Georg Klein
  Precision sigma = Precision(1.4826) *
      (1 + Precision(5) / (values.size() * 2 - 6)) * sqrt(median_squared);
  sigma = Precision(1.345) * sigma;
  sigmaSquared = sigma * sigma;
  return sigmaSquared;
}

template <typename Precision>
Precision Huber<Precision>::weight(Precision error_squared) const {
  if (error_squared > sigmaSquared) {
    ++outliers;
    return sqrt(sigmaSquared/error_squared);
  } else {
    return 1;
  }
}

#define INSTANTIATE_MESTIMATORS(T)\
template LeastSquares<T>::LeastSquares();\
template T LeastSquares<T>::compute_sigma_squared(std::vector<T>);\
template Tukey<T>::Tukey();\
template T Tukey<T>::compute_sigma_squared(std::vector<T>);\
template T Tukey<T>::weight(T) const;\
template Huber<T>::Huber();\
template T Huber<T>::compute_sigma_squared(std::vector<T>);\
template T Huber<T>::weight(T) const;

INSTANTIATE_MESTIMATORS(float)
INSTANTIATE_MESTIMATORS(double)

}  // namespace look3d
