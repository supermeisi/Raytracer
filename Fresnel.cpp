#include "Fresnel.h"

#include <algorithm>
#include <cmath>

FresnelSplit ComputeFresnelSplit(const std::array<double, 3> &incidentDirection,
                                 const std::array<double, 3> &normalAgainstRay,
                                 double incidentIndex,
                                 double transmittedIndex) {
  const double cosIncident =
      std::clamp(-(incidentDirection[0] * normalAgainstRay[0] +
                   incidentDirection[1] * normalAgainstRay[1] +
                   incidentDirection[2] * normalAgainstRay[2]),
                 0.0, 1.0);
  const double eta = incidentIndex / transmittedIndex;
  const double sinTransmittedSquared =
      eta * eta * (1.0 - cosIncident * cosIncident);

  const std::array<double, 3> reflected{
      incidentDirection[0] + 2.0 * cosIncident * normalAgainstRay[0],
      incidentDirection[1] + 2.0 * cosIncident * normalAgainstRay[1],
      incidentDirection[2] + 2.0 * cosIncident * normalAgainstRay[2]};
  if (sinTransmittedSquared >= 1.0) {
    return {1.0, 0.0, false, reflected, {0.0, 0.0, 0.0}};
  }

  const double cosTransmitted = std::sqrt(1.0 - sinTransmittedSquared);
  const double sDenominator =
      incidentIndex * cosIncident + transmittedIndex * cosTransmitted;
  const double pDenominator =
      incidentIndex * cosTransmitted + transmittedIndex * cosIncident;
  const double sAmplitude =
      (incidentIndex * cosIncident - transmittedIndex * cosTransmitted) /
      sDenominator;
  const double pAmplitude =
      (incidentIndex * cosTransmitted - transmittedIndex * cosIncident) /
      pDenominator;
  const double reflectance = std::clamp(
      0.5 * (sAmplitude * sAmplitude + pAmplitude * pAmplitude), 0.0, 1.0);

  const std::array<double, 3> transmitted{
      eta * incidentDirection[0] +
          (eta * cosIncident - cosTransmitted) * normalAgainstRay[0],
      eta * incidentDirection[1] +
          (eta * cosIncident - cosTransmitted) * normalAgainstRay[1],
      eta * incidentDirection[2] +
          (eta * cosIncident - cosTransmitted) * normalAgainstRay[2]};
  return {reflectance, 1.0 - reflectance, true, reflected, transmitted};
}
