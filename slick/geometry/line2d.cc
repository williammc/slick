// Copyright 2012, 2013, 2014 The LooK3D Authors. All rights reserved.
#include "slick/geometry/line2d.h"

#include <algorithm>
#include <Eigen/Core>

namespace slick {
const SlickScalar Line2d::epsilon = 1.e-30;

Line2d::Line2d(const Eigen::Vector2d& pt1, const Eigen::Vector2d& pt2,
               bool linevec) {
  assert((pt1 -pt2).norm() > epsilon);
  if (linevec) {
    is_bounded = false;
    point1 = pt1;
    point2 = pt1 + pt2.normalized();
  } else {
    is_bounded = true;
  }
  point1 = pt1;
  point2 = pt2;
  calc_line_vector();

  score = 0;
  CalcParams();
}

Line2d::Line2d(SlickScalar a, SlickScalar b, SlickScalar c) {
  a = a;
  b = b;
  c = c;
  score = 0;
  // calculate two points on the line
  CalcPointsFromLineEquation();
}

bool Line2d::IntersectLines(const Line2d& line, Eigen::Vector2d& intersection) const {
  // intersection point of line defined by p1 & p2 and line defined by p3 & p4
  // http://local.wasp.uwa.edu.au/~pbourke/geometry/lineline2d/
  SlickScalar x1, x2, x3, x4, y1, y2, y3, y4;
  x1 = point1[0];
  y1 = point1[1];
  x2 = point2[0];
  y2 = point2[1];
  x3 = line.point1[0];
  y3 = line.point1[1];
  x4 = line.point2[0];
  y4 = line.point2[1];

  if (std::fabs((y4-y3)*(x2-x1)-(x4-x3)*(y2-y1)) < epsilon)
      // 	warning('line_intersect.m degenerate --dclee');
    return false;
  intersection[0] = x1 + (x2-x1) * ((x4-x3)*(y1-y3)-(y4-y3)*(x1-x3))/((y4-y3)*(x2-x1)-(x4-x3)*(y2-y1));
  intersection[1] = y1 + (y2-y1) * ((x4-x3)*(y1-y3)-(y4-y3)*(x1-x3))/((y4-y3)*(x2-x1)-(x4-x3)*(y2-y1));
  return true;
}

bool Line2d::CheckIntersectLines(const Line2d& line) const {
  // using namespace std;
  if (!is_bounded  && !line.is_bounded) {
    // cout<<"all unbounded  divisor= "<<a_ * line.b_ - line.a_ * b_<<endl;
    return !is_zero(a * line.b - line.a * b);
  }

  if (std::max(point1.x(), point2.x()) < std::min(line.point1.x(), line.point2.x()) ||
      std::min(point1.x(), point2.x()) > std::max(line.point1.x(), line.point2.x()) ||
      std::max(point1.y(), point2.y()) < std::min(line.point1.y(), line.point2.y()) ||
      std::min(point1.y(), point2.y()) > std::max(line.point1.y(), line.point2.y())) {
    return false;
  }

  if ((line.a*point1.x()+line.b*point1.y()+line.c)*
      (line.a*point2.x()+line.b*point2.y()+line.c) >= 0.0) return false;

  if ((a*line.point1.x()+b*line.point1.y()+c)*
      (a*line.point2.x()+b*line.point2.y()+c) >= 0.0) return false;

  return true;
}

int Line2d::GetQuadrant() const {
  if (line_vector.x() >= 0.0) {
    if (line_vector.y() >= 0.0)
      return 1;
    else
      return 4;
  } else {
    if (line_vector.y() >= 0.0)
      return 2;
    else
      return 3;
  }
}

void Line2d::CalcParams(void) {
  SlickScalar x1, y1, x2, y2;

  x1 = point1.x();
  y1 = point1.y();
  x2 = point2.x();
  y2 = point2.y();

  if (fabs(x1 - x2) < epsilon) {
    a = 1.0;
    b = 0.0;
    c = -x1;
  } else {
    a = -((y1 - y2)/(x1 - x2));
    b = 1.0;
    c = -(y1 + a*x1);
  }

  return;
}

void Line2d::CalcPointsFromLineEquation(int offset) {
  if (is_zero(b)) {  // vertical line
    point1.x() = -c;
    point1.y() = 0.0;
    point2.x() = -c;
    point2.y() = offset;
  } else {
    point1.x() = 0.0;
    point1.y() = -c/b;
    point2.x() = offset;
    point2.y() = -(offset * a + c)/b;
  }

  calc_line_vector();
  is_bounded = false;
}
}  // namespace slick
