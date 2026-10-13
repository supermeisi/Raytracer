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

Vector Cross(const Vector &left, const Vector &right) {
  return {left[1] * right[2] - left[2] * right[1],
          left[2] * right[0] - left[0] * right[2],
          left[0] * right[1] - left[1] * right[0]};
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
  int currentObject;
  double intensity;
  int boundaryInteractions;
};

bool IntersectObject(const ObjectDescription &object, const Vector &origin,
                     const Vector &direction, double &distance, Vector &normal) {
  const Vector localOrigin = Subtract(origin, object.center);
  double nearest = std::numeric_limits<double>::infinity();
  Vector nearestNormal{};
  auto consider = [&](double candidate, const Vector &candidateNormal) {
    if (candidate > kSurfaceEpsilon && candidate < nearest) {
      nearest = candidate;
      nearestNormal = Unit(candidateNormal);
    }
  };

  if (object.shape == ObjectShape::Box) {
    double nearDistance = -std::numeric_limits<double>::infinity();
    double farDistance = std::numeric_limits<double>::infinity();
    Vector nearNormal{};
    Vector farNormal{};
    for (std::size_t axis = 0; axis < 3; ++axis) {
      const double halfSize = object.size[axis] * 0.5;
      if (std::abs(direction[axis]) <= 1e-14) {
        if (localOrigin[axis] < -halfSize || localOrigin[axis] > halfSize) {
          return false;
        }
        continue;
      }
      double first = (-halfSize - localOrigin[axis]) / direction[axis];
      double second = (halfSize - localOrigin[axis]) / direction[axis];
      Vector firstNormal{};
      Vector secondNormal{};
      firstNormal[axis] = -1.0;
      secondNormal[axis] = 1.0;
      if (first > second) {
        std::swap(first, second);
        std::swap(firstNormal, secondNormal);
      }
      if (first > nearDistance) {
        nearDistance = first;
        nearNormal = firstNormal;
      }
      if (second < farDistance) {
        farDistance = second;
        farNormal = secondNormal;
      }
      if (nearDistance > farDistance) {
        return false;
      }
    }
    if (nearDistance > kSurfaceEpsilon) {
      consider(nearDistance, nearNormal);
    } else {
      consider(farDistance, farNormal);
    }
  } else if (object.shape == ObjectShape::Sphere) {
    const double halfB = Dot(localOrigin, direction);
    const double c = Dot(localOrigin, localOrigin) -
                     object.radius * object.radius;
    const double discriminant = halfB * halfB - c;
    if (discriminant < 0.0) {
      return false;
    }
    const double root = std::sqrt(discriminant);
    const double nearDistance = -halfB - root;
    const double farDistance = -halfB + root;
    const double candidate =
        nearDistance > kSurfaceEpsilon ? nearDistance : farDistance;
    if (candidate <= kSurfaceEpsilon) {
      return false;
    }
    const Vector hit = Add(localOrigin, Scale(direction, candidate));
    consider(candidate, Scale(hit, 1.0 / object.radius));
  } else if (object.shape == ObjectShape::Step) {
    for (const MeshTriangle &triangle : object.triangles) {
      const Vector vertexA = Add(object.center, triangle.a);
      const Vector edge1 = Subtract(triangle.b, triangle.a);
      const Vector edge2 = Subtract(triangle.c, triangle.a);
      const Vector p = Cross(direction, edge2);
      const double determinant = Dot(edge1, p);
      if (std::abs(determinant) <= 1e-12) {
        continue;
      }
      const double inverseDeterminant = 1.0 / determinant;
      const Vector fromVertex = Subtract(origin, vertexA);
      const double u = Dot(fromVertex, p) * inverseDeterminant;
      if (u < 0.0 || u > 1.0) {
        continue;
      }
      const Vector q = Cross(fromVertex, edge1);
      const double v = Dot(direction, q) * inverseDeterminant;
      if (v < 0.0 || u + v > 1.0) {
        continue;
      }
      const double candidate = Dot(edge2, q) * inverseDeterminant;
      const Vector triangleNormal = Cross(edge1, edge2);
      if (std::hypot(triangleNormal[0], triangleNormal[1],
                     triangleNormal[2]) > 1e-14) {
        consider(candidate, triangleNormal);
      }
    }
  } else {
    const double halfHeight = object.size[2] * 0.5;
    const bool isTube = object.shape == ObjectShape::Tube;
    const double slope =
        isTube ? 0.0 : (object.topRadius - object.radius) / object.size[2];
    const double intercept = isTube ? object.radius
                                    : object.radius + slope * halfHeight;
    const double radialA = direction[0] * direction[0] +
                           direction[1] * direction[1];
    const double radialB = 2.0 * (localOrigin[0] * direction[0] +
                                  localOrigin[1] * direction[1]);
    const double radialC = localOrigin[0] * localOrigin[0] +
                           localOrigin[1] * localOrigin[1];
    const double a = radialA - slope * slope * direction[2] * direction[2];
    const double b = radialB -
                     2.0 * slope * direction[2] *
                         (slope * localOrigin[2] + intercept);
    const double c =
        radialC - std::pow(slope * localOrigin[2] + intercept, 2.0);
    if (std::abs(a) <= 1e-14) {
      if (std::abs(b) > 1e-14) {
        const double t = -c / b;
        const Vector hit = Add(localOrigin, Scale(direction, t));
        if (hit[2] >= -halfHeight - kSurfaceEpsilon &&
            hit[2] <= halfHeight + kSurfaceEpsilon) {
          const double radius = slope * hit[2] + intercept;
          Vector sideNormal{hit[0], hit[1], -radius * slope};
          if (std::hypot(hit[0], hit[1]) <= kSurfaceEpsilon &&
              std::abs(slope) > 1e-14) {
            sideNormal = {0.0, 0.0, -slope};
          }
          consider(t, sideNormal);
        }
      }
    } else {
      const double discriminant = b * b - 4.0 * a * c;
      if (discriminant >= 0.0) {
        const double root = std::sqrt(discriminant);
        for (double t : {(-b - root) / (2.0 * a),
                         (-b + root) / (2.0 * a)}) {
          const Vector hit = Add(localOrigin, Scale(direction, t));
          if (hit[2] >= -halfHeight - kSurfaceEpsilon &&
              hit[2] <= halfHeight + kSurfaceEpsilon) {
            const double radius = slope * hit[2] + intercept;
            Vector sideNormal{hit[0], hit[1], -radius * slope};
            if (std::hypot(hit[0], hit[1]) <= kSurfaceEpsilon &&
                std::abs(slope) > 1e-14) {
              sideNormal = {0.0, 0.0, -slope};
            }
            consider(t, sideNormal);
          }
        }
      }
    }

    for (int side : {-1, 1}) {
      if (std::abs(direction[2]) <= 1e-14) {
        continue;
      }
      const double z = side * halfHeight;
      const double t = (z - localOrigin[2]) / direction[2];
      const double x = localOrigin[0] + t * direction[0];
      const double y = localOrigin[1] + t * direction[1];
      const double radius = side < 0 ? object.radius
                                     : (isTube ? object.radius
                                               : object.topRadius);
      if (x * x + y * y <= radius * radius + kSurfaceEpsilon) {
        consider(t, {0.0, 0.0, static_cast<double>(side)});
      }
    }
  }

  if (!std::isfinite(nearest)) {
    return false;
  }
  distance = nearest;
  normal = nearestNormal;
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
      int hitObject = -1;
      if (state.currentObject >= 0) {
        const ObjectDescription &object = scene.objects[state.currentObject];
        if (IntersectObject(object, state.origin, state.direction, hitDistance,
                            outwardNormal)) {
          hitObject = state.currentObject;
        }
      } else {
        hitDistance = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < scene.objects.size(); ++i) {
          double candidateDistance = 0.0;
          Vector candidateNormal{};
          const ObjectDescription &object = scene.objects[i];
          if (IntersectObject(object, state.origin, state.direction,
                              candidateDistance, candidateNormal) &&
              candidateDistance < hitDistance) {
            hitDistance = candidateDistance;
            outwardNormal = candidateNormal;
            hitObject = static_cast<int>(i);
          }
        }
      }

      if (detector != nullptr) {
        const double maxDistance =
            hitObject >= 0 ? std::min(remaining, hitDistance) : remaining;
        const auto detectorHit =
            detector->Intersect(state.origin, state.direction, maxDistance);
        if (detectorHit) {
          const Vector hit = Add(state.origin,
                                 Scale(state.direction, detectorHit->distance));
          result.segments.push_back(
              {state.origin, hit, state.currentObject, state.intensity});
          result.detectorHit = *detectorHit;
          result.detectorIntensity = state.intensity;
          break;
        }
      }

      if (hitObject < 0 || hitDistance >= remaining) {
        result.segments.push_back(
            {state.origin, Add(state.origin, Scale(state.direction, remaining)),
             state.currentObject, state.intensity});
        break;
      }

      const Vector hit = Add(state.origin, Scale(state.direction, hitDistance));
      result.segments.push_back(
          {state.origin, hit, state.currentObject, state.intensity});
      state.travelled += hitDistance;
      const bool entering = state.currentObject < 0;
      const Vector normalAgainstRay =
          entering ? outwardNormal : Scale(outwardNormal, -1.0);
      const int transmittedObject = entering ? hitObject : -1;
      const double incidentIndex =
          entering ? 1.0
                   : scene.objects[state.currentObject].refractiveIndex;
      const double transmittedIndex =
          entering ? scene.objects[transmittedObject].refractiveIndex : 1.0;
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
        transmitted.currentObject = transmittedObject;
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
