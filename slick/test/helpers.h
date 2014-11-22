//! convenient funcs
#ifndef LOOK3D_TEST_HELPERS_H_
#define LOOK3D_TEST_HELPERS_H_
#include <Eigen/Core>
#include <opencv2/core/core.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <math/se3.h>
#include <math/utilities.h>
#include <track/types.h>

namespace look3d {

// holding data for synthesis point
struct SynthesisPoint {
  SynthesisPoint(){}
  SynthesisPoint(Eigen::Matrix<DefaultScalarType, 3, 1> c, float radius, cv::Scalar color, float thick)
    : center(c), color(color), thickness(thick) {
    calcInternalData(radius);
  }

  void calcInternalData(float radius=2.0) {
    top = center; top[1] = top[1]-radius;
    bottom = center; bottom[1] = bottom[1]+radius;
    left = center; left[0] = left[0]-radius;
    right = center; right[0] = right[0]+radius;

    normal = (top-center).cross(left-center).normalized();
  }

  SynthesisPoint transformToFrame(const SE3 &pose) const {
    SynthesisPoint p = *this;
    p.center = pose*center;
    p.top = pose*top;
    p.bottom = pose*bottom;
    p.left = pose*left;
    p.right = pose*right;
    p.thickness = thickness*p.center[2]/center[2];
    return p;
  }

  template<typename Camera>
  void drawOnFrame(const SE3&pose, const Camera& linear_cam, cv::Mat& frame) const {
    SynthesisPoint p_cam = transformToFrame(pose);
    Eigen::Vector2d v2_top = linear_cam.Project(project(p_cam.top));
    Eigen::Vector2d v2_bottom = linear_cam.Project(project(p_cam.bottom));
    Eigen::Vector2d v2_left = linear_cam.Project(project(p_cam.left));
    Eigen::Vector2d v2_right = linear_cam.Project(project(p_cam.right));
    cv::line(frame, cv::Point(v2_top[0],v2_top[1]), cv::Point(v2_bottom[0],v2_bottom[1]),
             p_cam.color, std::max<float>(2.f,p_cam.thickness));
    cv::line(frame, cv::Point(v2_left[0],v2_left[1]), cv::Point(v2_right[0],v2_right[1]),
             p_cam.color, std::max<float>(2.f,p_cam.thickness));
  }

  Eigen::Matrix<DefaultScalarType, 3, 1> center;
  Eigen::Matrix<DefaultScalarType, 3, 1> top;
  Eigen::Matrix<DefaultScalarType, 3, 1> bottom;
  Eigen::Matrix<DefaultScalarType, 3, 1> left;
  Eigen::Matrix<DefaultScalarType, 3, 1> right;

  Eigen::Matrix<DefaultScalarType, 3, 1> normal;

  cv::Scalar color;
  float thickness;

};

// holding 3D synthesis points
struct SynthesisModel {
  std::vector<SynthesisPoint> world_points;
};

// check if a synthesis point (3D point) is visible in given keyframe
bool check3DPointVisible(const Eigen::Matrix<DefaultScalarType, 3, 1>& point, KeyframeType* kf) {
  Eigen::Vector2d v2_implane = kf->project(project(kf->pose*point));
  return(in_image_with_border(kf->scale_images[0], cv::Point(v2_implane[0], v2_implane[1]), 1));
}

//! generate 3D points in square on given z=1 plane then scale them
void generate3DPointsOnPlane(SynthesisModel& model, const int npoints=100,
                             const DefaultScalarType scale_xyz=1.0) {
  Eigen::Matrix<DefaultScalarType, 3, 1> point(0,0,1);
  for (int i = 0; i < npoints; ++i) {
    point[0] = DefaultScalarType(std::rand()%1000-500)/500;
    point[1] = DefaultScalarType(std::rand()%1000-500)/500;
    SynthesisPoint p(point*scale_xyz, DefaultScalarType(std::rand()%100)/1000,
                     cv::Scalar(100+std::rand()%155,100+std::rand()%155,100+std::rand()%155),
                     DefaultScalarType(std::rand()%1000)/200);
    model.world_points.push_back(p);
  }
}

// generate 3D points in square on given z=1 plane then scale them
void generate3DPointsOffPlane(SynthesisModel& model, const int npoints=100,
                              const DefaultScalarType scale_xyz=1.0) {
  Eigen::Matrix<DefaultScalarType, 3, 1> point(0,0,1);
  for (int i = 0; i < npoints; ++i) {
    point[0] = DefaultScalarType(std::rand()%1000-500)/500;
    point[1] = DefaultScalarType(std::rand()%1000-500)/500;
    point[2] = 0.5 + DefaultScalarType(std::rand()%1000)/500;
    SynthesisPoint p(point*scale_xyz, DefaultScalarType(std::rand()%100)/1000,
                     cv::Scalar(100+std::rand()%155,100+std::rand()%155,100+std::rand()%155),
                     DefaultScalarType(std::rand()%1000)/200);
    model.world_points.push_back(p);
  }
}

//! drawing a synthesis frame given pose
template<typename Camera>
void generateSynthesisFrame(cv::Mat &syn_frame, const SynthesisModel& model,
                       const SE3& pose, const Camera& linear_cam) {
  int w, h;
  linear_cam.get_resolution(w, h);
  if (syn_frame.empty())
    syn_frame.create(h,w, CV_8UC1);
  syn_frame.setTo(cv::Scalar(0));
  for (int i = 0; i < model.world_points.size(); ++i) {
    model.world_points[i].drawOnFrame(pose, linear_cam, syn_frame);
  }
}

//! generates keypoints in given frame pose
template<typename Camera>
void generateKeyPointsOnFrame(std::vector<cv::KeyPoint>& keypoints,
                              std::map<int, const SynthesisPoint*>& map_kpid_modelpoint,
                              const SynthesisModel& model, SE3&pose,
                              Camera& linear_cam, cv::Mat& frame) {
  keypoints.clear();
  for (int i = 0; i < model.world_points.size(); ++i) {
    SynthesisPoint p_cam = model.world_points[i].transformToFrame(pose);
    Eigen::Vector2d v2_implane = linear_cam.Project(project(p_cam.center));
    if (in_image_with_border(frame, cv::Point(v2_implane[0],v2_implane[1]),1)) {
      keypoints.push_back(cv::KeyPoint(cv::Point2f(v2_implane[0], v2_implane[1]), 8.f));
      map_kpid_modelpoint[keypoints.size()-1] = &model.world_points[i];
    }
  }
}
typedef std::pair<Eigen::Vector2d, Eigen::Vector2d> Vector2Pair;
// generates matches between source keyframe & target keyframe
void generateSourceTargetMatches(std::vector<Vector2Pair, Eigen::aligned_allocator<Vector2Pair> >& matches,
                     std::map<int, const SynthesisPoint*>& map_matchid_modelpoint,
                     const SynthesisModel& model,
                     KeyframeType* kf_source, KeyframeType* kf_target) {
  matches.clear();
  for (size_t i = 0; i < model.world_points.size(); ++i) {
    if ( check3DPointVisible(model.world_points[i].center, kf_source ) &&
        check3DPointVisible(model.world_points[i].center, kf_target)) {
      Eigen::Vector2d v2_source = kf_source->project(unproject(model.world_points[i].center));
      Eigen::Vector2d v2_target = kf_target->project(unproject(model.world_points[i].center));
      matches.push_back(std::make_pair(v2_source, v2_target));
      map_matchid_modelpoint[matches.size()-1] = &model.world_points[i];
    }
  }
}

// generate tracking model (include features, keyframes) from synthesis data
template<typename Camera>
void generateTrackingModel(ModelType& model,
                           const SynthesisModel& syn_model,
                           const std::vector<SE3 >& poses, Camera& linear_cam) {
  std::map<const SynthesisPoint*, PointFeatureType*> map_synpoint_pointfeature;
  for (size_t i = 0; i < poses.size(); ++i) {
    cv::Mat syn_frame;
    generateSynthesisFrame(syn_frame, syn_model, poses[i], linear_cam);
    KeyframeType::Ptr kf_temp(new KeyframeType(syn_frame, 4, linear_cam));
    kf_temp->pose = poses[i];
    model.frames.push_back(kf_temp);
    for (size_t j = 0; j < syn_model.world_points.size(); ++j) {
      SynthesisPoint p_cam = syn_model.world_points[j].transformToFrame(poses[i]);
      Eigen::Vector2d v2_implane = linear_cam.Project(project(p_cam.center));
      if (in_image_with_border(syn_frame, cv::Point(v2_implane[0],v2_implane[1]),1)) {

        PointFeatureType* pf;
        if (map_synpoint_pointfeature.count(&syn_model.world_points[j]) == 0) {
          model.features.push_back(PointFeatureType::Ptr(
                                     new PointFeatureType(kf_temp.get(), v2_implane, syn_model.world_points[j].center, syn_model.world_points[j].normal)));
          pf = model.features.back().get();
          map_synpoint_pointfeature[&syn_model.world_points[j]] = model.features.back().get();
        } else {
          pf = map_synpoint_pointfeature[&syn_model.world_points[j]];
        }

        Measurement meas((Measurement::FeatureType*)pf,(Measurement::FrameType*)kf_temp.get(), v2_implane);
        kf_temp->observations.push_back(meas);
      }
    }
  }
}

// draw Cross of given point in an image
void drawCross(cv::Mat& image, cv::Point& center, int radius,
               cv::Scalar color, int thick = 1) {
  cv::Point top(center); top.y -=radius;
  cv::Point bottom(center); bottom.y +=radius;
  cv::Point left(center); left.x -=radius;
  cv::Point right(center); right.x +=radius;
  cv::line(image, top, bottom, color, thick);
}
}  // namespace look3d
#endif  // LOOK3D_TEST_HELPERS_H_
