#ifndef ROOT_RAY_TRACER_FRESNEL_H
#define ROOT_RAY_TRACER_FRESNEL_H

#include <array>

struct FresnelSplit {
  double reflectance;
  double transmittance;
  bool hasTransmission;
  std::array<double, 3> reflectedDirection;
  std::array<double, 3> transmittedDirection;
};

FresnelSplit ComputeFresnelSplit(const std::array<double, 3> &incidentDirection,
                                 const std::array<double, 3> &normalAgainstRay,
                                 double incidentIndex, double transmittedIndex);

#endif
