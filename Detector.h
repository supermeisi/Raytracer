#ifndef ROOT_RAY_TRACER_DETECTOR_H
#define ROOT_RAY_TRACER_DETECTOR_H

#include "SceneConfig.h"

#include <array>
#include <optional>
#include <vector>

struct DetectorHit {
  double distance;
  double u;
  double v;
  int binX;
  int binY;
};

class PlaneDetectorHistogram {
public:
  explicit PlaneDetectorHistogram(const PlaneDetectorDescription &detector);

  std::optional<DetectorHit> Intersect(const std::array<double, 3> &origin,
                                       const std::array<double, 3> &direction,
                                       double maximumDistance) const;
  void Record(const DetectorHit &hit, double intensity);

  const PlaneDetectorDescription &Description() const;
  const std::vector<double> &Bins() const;
  double TotalIntensity() const;

private:
  PlaneDetectorDescription detector_;
  std::array<double, 3> normal_;
  std::array<double, 3> uAxis_;
  std::array<double, 3> vAxis_;
  std::vector<double> bins_;
  double totalIntensity_ = 0.0;
};

#endif
