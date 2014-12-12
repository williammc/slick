// Copyright 2014 The Slick Authors. All rights reserved.
#pragma once
#include <vector>
#include <cassert>
#include <algorithm>

namespace slick {

template <typename Scalar>
struct MEstimator {
  MEstimator() {}
  virtual ~MEstimator() {}
  typedef Scalar value_type;
  virtual Scalar compute_sigma_squared(std::vector<Scalar> values) = 0;
  virtual Scalar weight(Scalar error_squared) const = 0;
  virtual int get_outliers() const = 0;
  virtual Scalar get_sigma_squared() const = 0;
};

template <typename Scalar>
struct LeastSquares : public MEstimator<Scalar> {
  LeastSquares() : sq_sigma(0) {}
  Scalar sq_sigma;

  Scalar compute_sigma_squared(std::vector<Scalar> values) {
    sq_sigma = 0;
    if (values.size() <= 0)
      return sq_sigma;
    for (typename std::vector<Scalar>::const_iterator i = values.begin();
        i != values.end(); ++i)
      sq_sigma += *i;
    sq_sigma /= values.size();
    return sq_sigma;
  }

  Scalar weight(Scalar error_squared) const { return 1; }
  int get_outliers() const { return 0; }
  Scalar get_sigma_squared() const { return sq_sigma; }
};

template <typename Scalar>
struct Tukey : public MEstimator<Scalar> {
  Scalar sq_sigma;
  mutable int outliers;

  Tukey() : sq_sigma(0), outliers(0) {}

  Scalar compute_sigma_squared(std::vector<Scalar> values) {
    if (values.size() == 0)
      return 0;
    std::sort(values.begin(), values.end());
    Scalar median_squared = values[values.size()/2];
    if (median_squared == 0)
      median_squared = Scalar(1.e-19);
    // from Georg Klein
    Scalar sigma = Scalar(1.4826) *
        (1 + Scalar(5) / (values.size() * 2 - 6)) * sqrt(median_squared);
    sigma = Scalar(4.6851) * sigma;
    sq_sigma = sigma * sigma;
    return sq_sigma;
  }

  Scalar weight(Scalar error_squared) const {
    if (error_squared > sq_sigma) {
      ++outliers;
      return 0;
    } else {
      const Scalar w = 1 - error_squared / sq_sigma;
      return w*w;
    }
  }

  int get_outliers() const { return outliers; }
  Scalar get_sigma_squared() const { return sq_sigma; }
};

template <typename Scalar>
struct Huber : public MEstimator<Scalar> {
  Scalar sq_sigma;
  mutable int outliers;

  Huber() : sq_sigma(0), outliers(0) {}

  Scalar compute_sigma_squared(std::vector<Scalar> values) {
    assert(values.size() > 0);
    std::sort(values.begin(), values.end());
    Scalar median_squared = values[values.size()/2];
    // from Georg Klein
    Scalar sigma = Scalar(1.4826) *
        (1 + Scalar(5) / (values.size() * 2 - 6)) * sqrt(median_squared);
    sigma = Scalar(1.345) * sigma;
    sq_sigma = sigma * sigma;
    return sq_sigma;
  }

  Scalar weight(Scalar error_squared) const {
    if (error_squared > sq_sigma) {
      ++outliers;
      return sqrt(sq_sigma/error_squared);
    } else {
      return 1;
    }
  }

  int get_outliers() const { return outliers; }
  Scalar get_sigma_squared() const { return sq_sigma; }
};
}  // namespace slick
