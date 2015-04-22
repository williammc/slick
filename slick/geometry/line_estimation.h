#pragma once

namespace slick {

// Estimate 2D line from list of 2D points.
/// @return line and strength ([0,1] weak to strong)
template <typename T, typename Vec2Type>
inline std::pair<Line2DBase<T>, T>
LeastSquareFitLine(const std::vector<Vec2Type> &points) {
  // calculate eigen value to find lines
  int left = std::numeric_limits<int>::max(), right = 0,
      top = std::numeric_limits<int>::max(), bottom = 0;
  T sum_x = 0, sum_y = 0;
  for (auto pt : points) {
    sum_x += pt.x();
    sum_y += pt.y();
    left = std::min<int>(left, pt.x());
    right = std::max<int>(right, pt.x());
    top = std::min<int>(top, pt.y());
    bottom = std::max<int>(bottom, pt.y());
  }
  const T avg_x = sum_x / points.size();
  const T avg_y = sum_y / points.size();
  // check if line support region is long enough
  Eigen::Matrix<T, 2, 2> m2_cov = Eigen::Matrix<T, 2, 2>::Zero();
  for (auto pt : points) {
    m2_cov(0, 0) += (pt.x() - avg_x) * (pt.x() - avg_x);
    m2_cov(0, 1) += (pt.x() - avg_x) * (pt.y() - avg_y);
    m2_cov(1, 1) += (pt.y() - avg_y) * (pt.y() - avg_y);
  }
  m2_cov(1, 0) = m2_cov(0, 1);
  Eigen::EigenSolver<Eigen::Matrix<T, 2, 2>> eig_solv(m2_cov);
  T confidence =
      eig_solv.eigenvalues()[0].real() / eig_solv.eigenvalues()[1].real();
  T theta = atan2(eig_solv.eigenvectors().col(1)[1].real(),
                  eig_solv.eigenvectors().col(1)[0].real());
  if (eig_solv.eigenvalues()[0].real() < eig_solv.eigenvalues()[1].real()) {
    confidence = 1.0 / confidence;
    theta = atan2(eig_solv.eigenvectors().col(0)[1].real(),
                  eig_solv.eigenvectors().col(0)[0].real());
  }

  // confidence ranges [5, 1000] bad to good
  T strength = (confidence - 5) / (1000 - 5);

  Line2DBase<T> line_seg;
  // range 0.009 to 0.14 from good to bad
  if (theta == 0.0) { // purely vertical line
    theta = 1.e-18;
    line_seg = Line2DBase<T>(Eigen::Vector2d(avg_x, top),
                             Eigen::Vector2d(avg_x, bottom));
  } else {
    const T left_y = (avg_x - left) * (1 / tan(theta)) + avg_y;
    Eigen::Vector2d p1, p2;
    if (left_y < top) {
      const T top_x = (avg_y - top) * tan(theta) + avg_x;
      p1[0] = top_x;
      p1[1] = top;
    } else if (left_y > bottom) {
      const T bottom_x = (avg_y - bottom) * tan(theta) + avg_x;
      p1[0] = bottom_x;
      p1[1] = bottom;
    } else {
      p1[0] = left;
      p1[1] = left_y;
    }
    const T right_y = (avg_x - right) * (1 / tan(theta)) + avg_y;
    if (right_y < top) {
      const T top_x = (avg_y - top) * tan(theta) + avg_x;
      p2[0] = top_x;
      p2[1] = top;
    } else if (right_y > bottom) {
      const T bottom_x = (avg_y - bottom) * tan(theta) + avg_x;
      p2[0] = bottom_x;
      p2[1] = bottom;
    } else {
      p2[0] = right;
      p2[1] = right_y;
    }
    line_seg = Line2DBase<T>(p1, p2);
  }
  return std::make_pair(line_seg, strength);
}

} // namespace slick