#ifndef ROOT_RAY_TRACER_RAY_TRACING_ENGINE_H
#define ROOT_RAY_TRACER_RAY_TRACING_ENGINE_H

#include "SceneConfig.h"
#include "Detector.h"

#include <array>
#include <atomic>
#include <cstddef>
#include <optional>
#include <vector>

struct TracedRaySegment {
  std::array<double, 3> start;
  std::array<double, 3> end;
  int sphereIndex;
  double intensity;
};

struct TracedRay {
  std::vector<TracedRaySegment> segments;
  std::optional<DetectorHit> detectorHit;
  double detectorIntensity = 0.0;
};

std::vector<TracedRay>
TraceRaysParallel(const SceneConfig &scene,
                  const std::vector<std::array<double, 3>> &directions,
                  int workerCount, std::atomic<std::size_t> &completedRayCount);

#endif
