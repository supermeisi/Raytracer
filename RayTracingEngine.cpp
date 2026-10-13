#include "RayTracingEngine.h"

#include "Detector.h"
#include "Fresnel.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <thread>

namespace {

constexpr double kSurfaceEpsilon = 1e-5;
using Vector = std::array<double, 3>;

Vector Add(const Vector &left, const Vector &right) {
  return {left[0] + right[0], left[1] + right[1], left[2] + right[2]};
}

Vector Subtract(const Vector &left, const Vector &right) {
  return {left[0] - right[0], left[1] - right[1], left[2] - right[2]};
}

Vector Scale(const Vector &value, double scale) {
  return {value[0] * scale, value[1] * scale, value[2] * scale};
}

double Dot(const Vector &left, const Vector &right) {
  return left[0] * right[0] + left[1] * right[1] + left[2] * right[2];
}

Vector Unit(const Vector &value) {
  const double magnitude = std::hypot(value[0], value[1], value[2]);
  return Scale(value, 1.0 / magnitude);
}

struct RayState {
  Vector origin;
  Vector direction;
  double travelled;
  int currentSphere;
  double intensity;
  int boundaryInteractions;
};

bool IntersectSphere(const Vector &origin, const Vector &direction,
                     const Vector &center, double radius, double &distance,
                     Vector &normal) {
  const Vector offset = Subtract(origin, center);
  const double halfB = Dot(offset, direction);
  const double c = Dot(offset, offset) - radius * radius;
  const double discriminant = halfB * halfB - c;
  if (discriminant < 0.0) {
    return false;
  }

  const double root = std::sqrt(discriminant);
  const double nearDistance = -halfB - root;
  const double farDistance = -halfB + root;
  distance = nearDistance > kSurfaceEpsilon ? nearDistance : farDistance;
  if (distance <= kSurfaceEpsilon) {
    return false;
  }

  normal = Scale(Subtract(Add(origin, Scale(direction, distance)), center),
                 1.0 / radius);
  return true;
}

TracedRay TraceRay(const SceneConfig &scene, const Vector &source,
                   const Vector &direction,
                   const PlaneDetectorHistogram *detector) {
  std::vector<RayState> pending{
      {source, Unit(direction), 0.0, -1, 1.0, 0}};
  TracedRay result;

  while (!pending.empty()) {
    RayState state = pending.back();
    pending.pop_back();
    while (state.travelled < scene.pathLength &&
           state.boundaryInteractions < scene.maxBoundaryInteractions) {
      double hitDistance = 0.0;
      Vector outwardNormal{};
      const double remaining = scene.pathLength - state.travelled;
      int hitSphere = -1;
      if (state.currentSphere >= 0) {
        const SphereDescription &sphere = scene.spheres[state.currentSphere];
        if (IntersectSphere(state.origin, state.direction, sphere.center,
                            sphere.radius, hitDistance, outwardNormal)) {
          hitSphere = state.currentSphere;
        }
      } else {
        hitDistance = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < scene.spheres.size(); ++i) {
          double candidateDistance = 0.0;
          Vector candidateNormal{};
          const SphereDescription &sphere = scene.spheres[i];
          if (IntersectSphere(state.origin, state.direction, sphere.center,
                              sphere.radius, candidateDistance,
                              candidateNormal) &&
              candidateDistance < hitDistance) {
            hitDistance = candidateDistance;
            outwardNormal = candidateNormal;
            hitSphere = static_cast<int>(i);
          }
        }
      }

      if (detector != nullptr) {
        const double maxDistance =
            hitSphere >= 0 ? std::min(remaining, hitDistance) : remaining;
        const auto detectorHit =
            detector->Intersect(state.origin, state.direction, maxDistance);
        if (detectorHit) {
          const Vector hit = Add(state.origin,
                                 Scale(state.direction, detectorHit->distance));
          result.segments.push_back(
              {state.origin, hit, state.currentSphere, state.intensity});
          result.detectorHit = *detectorHit;
          result.detectorIntensity = state.intensity;
          break;
        }
      }

      if (hitSphere < 0 || hitDistance >= remaining) {
        result.segments.push_back(
            {state.origin, Add(state.origin, Scale(state.direction, remaining)),
             state.currentSphere, state.intensity});
        break;
      }

      const Vector hit = Add(state.origin, Scale(state.direction, hitDistance));
      result.segments.push_back(
          {state.origin, hit, state.currentSphere, state.intensity});
      state.travelled += hitDistance;
      const bool entering = state.currentSphere < 0;
      const Vector normalAgainstRay =
          entering ? outwardNormal : Scale(outwardNormal, -1.0);
      const int transmittedSphere = entering ? hitSphere : -1;
      const double incidentIndex =
          entering ? 1.0 : scene.spheres[state.currentSphere].refractiveIndex;
      const double transmittedIndex =
          entering ? scene.spheres[transmittedSphere].refractiveIndex : 1.0;
      const FresnelSplit split = ComputeFresnelSplit(
          state.direction, normalAgainstRay, incidentIndex, transmittedIndex);

      ++state.boundaryInteractions;
      const double reflectedIntensity = state.intensity * split.reflectance;
      if (reflectedIntensity >= scene.minimumBranchIntensity) {
        RayState reflected = state;
        reflected.direction = split.reflectedDirection;
        reflected.origin =
            Add(hit, Scale(reflected.direction, kSurfaceEpsilon));
        reflected.travelled += kSurfaceEpsilon;
        reflected.intensity = reflectedIntensity;
        pending.push_back(reflected);
      }

      const double transmittedIntensity = state.intensity * split.transmittance;
      if (split.hasTransmission &&
          transmittedIntensity >= scene.minimumBranchIntensity) {
        RayState transmitted = state;
        transmitted.direction = split.transmittedDirection;
        transmitted.origin =
            Add(hit, Scale(transmitted.direction, kSurfaceEpsilon));
        transmitted.currentSphere = transmittedSphere;
        transmitted.travelled += kSurfaceEpsilon;
        transmitted.intensity = transmittedIntensity;
        pending.push_back(transmitted);
      }
    }
  }
  return result;
}

} // namespace

std::vector<TracedRay>
TraceRaysParallel(const SceneConfig &scene,
                  const std::vector<std::array<double, 3>> &directions,
                  int workerCount, std::atomic<std::size_t> &completedRayCount) {
  if (workerCount < 0) {
    throw std::invalid_argument("worker count cannot be negative");
  }
  std::vector<TracedRay> rays(directions.size());
  if (directions.empty()) {
    completedRayCount.store(0, std::memory_order_relaxed);
    return rays;
  }

  if (workerCount == 0) {
    workerCount = static_cast<int>(std::thread::hardware_concurrency());
    if (workerCount == 0) {
      workerCount = 1;
    }
  }
  const std::size_t threadCount =
      std::min<std::size_t>(directions.size(),
                            static_cast<std::size_t>(workerCount));
  std::atomic<std::size_t> nextRay{0};
  completedRayCount.store(0, std::memory_order_relaxed);
  const Vector source = scene.sourcePosition;
  std::optional<PlaneDetectorHistogram> detector;
  if (scene.detector) {
    detector.emplace(*scene.detector);
  }
  std::exception_ptr workerError;
  std::mutex errorMutex;
  std::vector<std::thread> workers;
  workers.reserve(threadCount);

  try {
    for (std::size_t worker = 0; worker < threadCount; ++worker) {
      workers.emplace_back([&] {
        try {
          for (;;) {
            const std::size_t index =
                nextRay.fetch_add(1, std::memory_order_relaxed);
            if (index >= directions.size()) {
              break;
            }
            rays[index] = TraceRay(scene, source, directions[index],
                                   detector ? &*detector : nullptr);
            completedRayCount.fetch_add(1, std::memory_order_relaxed);
          }
        } catch (...) {
          std::lock_guard<std::mutex> lock(errorMutex);
          if (!workerError) {
            workerError = std::current_exception();
          }
        }
      });
    }
  } catch (...) {
    for (std::thread &worker : workers) {
      worker.join();
    }
    throw;
  }
  for (std::thread &worker : workers) {
    worker.join();
  }
  if (workerError) {
    std::rethrow_exception(workerError);
  }
  return rays;
}
