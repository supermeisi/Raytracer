#include "Detector.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {

double Dot(const std::array<double, 3> &left,
           const std::array<double, 3> &right) {
  return left[0] * right[0] + left[1] * right[1] + left[2] * right[2];
}

std::array<double, 3> Subtract(const std::array<double, 3> &left,
                               const std::array<double, 3> &right) {
  return {left[0] - right[0], left[1] - right[1], left[2] - right[2]};
}

std::array<double, 3> Normalize(std::array<double, 3> vector) {
  const double magnitude = std::hypot(vector[0], vector[1], vector[2]);
  if (!std::isfinite(magnitude) || magnitude <= 0.0) {
    throw std::invalid_argument("detector basis vectors must be nonzero");
  }
  for (double &component : vector) {
    component /= magnitude;
  }
  return vector;
}

std::array<double, 3> Cross(const std::array<double, 3> &left,
                            const std::array<double, 3> &right) {
  return {left[1] * right[2] - left[2] * right[1],
          left[2] * right[0] - left[0] * right[2],
          left[0] * right[1] - left[1] * right[0]};
}

} // namespace

PlaneDetectorHistogram::PlaneDetectorHistogram(
    const PlaneDetectorDescription &detector)
    : detector_(detector), normal_(Normalize(detector.normal)),
      bins_(static_cast<std::size_t>(detector.binsX) *
                static_cast<std::size_t>(detector.binsY),
            0.0) {
  const double projection = Dot(detector.uAxis, normal_);
  uAxis_ = {detector.uAxis[0] - projection * normal_[0],
            detector.uAxis[1] - projection * normal_[1],
            detector.uAxis[2] - projection * normal_[2]};
  uAxis_ = Normalize(uAxis_);
  vAxis_ = Normalize(Cross(normal_, uAxis_));
}

std::optional<DetectorHit>
PlaneDetectorHistogram::Intersect(const std::array<double, 3> &origin,
                                  const std::array<double, 3> &direction,
                                  double maximumDistance) const {
  const double directionMagnitude =
      std::hypot(direction[0], direction[1], direction[2]);
  if (!std::isfinite(directionMagnitude) || directionMagnitude <= 0.0) {
    return std::nullopt;
  }
  const std::array<double, 3> unitDirection{direction[0] / directionMagnitude,
                                            direction[1] / directionMagnitude,
                                            direction[2] / directionMagnitude};
  const double denominator = Dot(unitDirection, normal_);
  if (std::abs(denominator) <= 1e-12) {
    return std::nullopt;
  }

  const std::array<double, 3> centerToOrigin =
      Subtract(detector_.center, origin);
  const double distance = Dot(centerToOrigin, normal_) / denominator;
  if (distance <= 1e-9 || distance > maximumDistance) {
    return std::nullopt;
  }

  const std::array<double, 3> point{origin[0] + distance * unitDirection[0],
                                    origin[1] + distance * unitDirection[1],
                                    origin[2] + distance * unitDirection[2]};
  const std::array<double, 3> centerToPoint = Subtract(point, detector_.center);
  const double u = Dot(centerToPoint, uAxis_);
  const double v = Dot(centerToPoint, vAxis_);
  if (std::abs(u) > detector_.width * 0.5 ||
      std::abs(v) > detector_.height * 0.5) {
    return std::nullopt;
  }

  const int binX =
      std::min(detector_.binsX - 1,
               static_cast<int>((u / detector_.width + 0.5) * detector_.binsX));
  const int binY = std::min(
      detector_.binsY - 1,
      static_cast<int>((v / detector_.height + 0.5) * detector_.binsY));
  return DetectorHit{distance, u, v, binX, binY};
}

void PlaneDetectorHistogram::Record(const DetectorHit &hit, double intensity) {
  if (!std::isfinite(intensity) || intensity < 0.0) {
    throw std::invalid_argument(
        "detected intensity must be finite and nonnegative");
  }
  const std::size_t index = static_cast<std::size_t>(hit.binY) *
                                static_cast<std::size_t>(detector_.binsX) +
                            static_cast<std::size_t>(hit.binX);
  bins_.at(index) += intensity;
  totalIntensity_ += intensity;
}

const PlaneDetectorDescription &PlaneDetectorHistogram::Description() const {
  return detector_;
}

const std::vector<double> &PlaneDetectorHistogram::Bins() const {
  return bins_;
}

double PlaneDetectorHistogram::TotalIntensity() const {
  return totalIntensity_;
}
