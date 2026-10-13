#include "SceneConfig.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

using Json = nlohmann::json;

double ReadFiniteNumber(const Json &value, const std::string &description) {
  if (!value.is_number()) {
    throw std::runtime_error(description + " must be a number");
  }
  const double number = value.get<double>();
  if (!std::isfinite(number)) {
    throw std::runtime_error(description + " must be finite");
  }
  return number;
}

std::array<double, 3> ReadVector3(const Json &value,
                                  const std::string &description) {
  if (!value.is_array() || value.size() != 3) {
    throw std::runtime_error(description + " must contain exactly 3 numbers");
  }
  return {ReadFiniteNumber(value[0], description + "[0]"),
          ReadFiniteNumber(value[1], description + "[1]"),
          ReadFiniteNumber(value[2], description + "[2]")};
}

double DistanceSquared(const std::array<double, 3> &left,
                       const std::array<double, 3> &right) {
  double distanceSquared = 0.0;
  for (std::size_t axis = 0; axis < 3; ++axis) {
    const double delta = left[axis] - right[axis];
    distanceSquared += delta * delta;
  }
  return distanceSquared;
}

void ValidateScene(const SceneConfig &scene) {
  if (scene.rayCount < 1 || scene.rayCount > 10000) {
    throw std::runtime_error("source.ray_count must be between 1 and 10000");
  }
  if (scene.openingAngleDegrees < 0.0 || scene.openingAngleDegrees > 180.0) {
    throw std::runtime_error(
        "source.opening_angle_degrees must be between 0 and 180");
  }
  const double directionMagnitude =
      std::hypot(scene.sourceDirection[0], scene.sourceDirection[1],
                 scene.sourceDirection[2]);
  if (!std::isfinite(directionMagnitude) || directionMagnitude <= 0.0) {
    throw std::runtime_error("source.direction must not be the zero vector");
  }
  if (scene.pathLength <= 0.0) {
    throw std::runtime_error("path_length must be greater than 0");
  }
  if (scene.minimumBranchIntensity < 0.0 ||
      scene.minimumBranchIntensity >= 1.0) {
    throw std::runtime_error(
        "optics.minimum_branch_intensity must be in [0, 1)");
  }
  if (scene.maxBoundaryInteractions < 1 || scene.maxBoundaryInteractions > 12) {
    throw std::runtime_error(
        "optics.max_boundary_interactions must be between 1 and 12");
  }
  if (scene.workerCount < 0 || scene.workerCount > 256) {
    throw std::runtime_error(
        "optics.worker_count must be between 0 (automatic) and 256");
  }
  if (scene.spheres.empty()) {
    throw std::runtime_error("objects must contain at least one sphere");
  }

  if (scene.detector) {
    const PlaneDetectorDescription &detector = *scene.detector;
    if (detector.width <= 0.0 || detector.height <= 0.0) {
      throw std::runtime_error(
          "detector width and height must be greater than 0");
    }
    if (detector.binsX < 1 || detector.binsX > 1000 || detector.binsY < 1 ||
        detector.binsY > 1000 ||
        static_cast<long long>(detector.binsX) * detector.binsY > 1000000) {
      throw std::runtime_error(
          "detector bins_x and bins_y must be 1-1000 with at most 1000000 "
          "total bins");
    }
    if (detector.color < 0 ||
        detector.color > std::numeric_limits<short>::max()) {
      throw std::runtime_error("detector.color is outside the ROOT color "
                               "index range");
    }
    const double normalMagnitude =
        std::hypot(detector.normal[0], detector.normal[1], detector.normal[2]);
    const double axisMagnitude =
        std::hypot(detector.uAxis[0], detector.uAxis[1], detector.uAxis[2]);
    if (!std::isfinite(normalMagnitude) || !std::isfinite(axisMagnitude) ||
        normalMagnitude <= 0.0 || axisMagnitude <= 0.0) {
      throw std::runtime_error(
          "detector.normal and detector.u_axis must be nonzero vectors");
    }
    const double axisNormalDot = (detector.uAxis[0] * detector.normal[0] +
                                  detector.uAxis[1] * detector.normal[1] +
                                  detector.uAxis[2] * detector.normal[2]) /
                                 axisMagnitude / normalMagnitude;
    if (std::abs(axisNormalDot) >= 1.0 - 1e-10) {
      throw std::runtime_error(
          "detector.u_axis must not be parallel to detector.normal");
    }
  }

  for (std::size_t i = 0; i < scene.spheres.size(); ++i) {
    const SphereDescription &sphere = scene.spheres[i];
    if (sphere.radius <= 0.0) {
      throw std::runtime_error("objects[" + std::to_string(i) +
                               "].radius must be greater than 0");
    }
    if (sphere.refractiveIndex <= 0.0) {
      throw std::runtime_error("objects[" + std::to_string(i) +
                               "].refractive_index must be greater than 0");
    }
    const double sourceDistanceSquared =
        DistanceSquared(scene.sourcePosition, sphere.center);
    if (sourceDistanceSquared <= sphere.radius * sphere.radius) {
      throw std::runtime_error("source.position must be outside every sphere");
    }

    for (std::size_t j = 0; j < i; ++j) {
      const SphereDescription &other = scene.spheres[j];
      const double minDistance = sphere.radius + other.radius;
      if (DistanceSquared(sphere.center, other.center) <=
          minDistance * minDistance) {
        throw std::runtime_error(
            "sphere interiors must not overlap or touch (objects[" +
            std::to_string(j) + "] and objects[" + std::to_string(i) + "])");
      }
    }
  }
}

} // namespace

SceneConfig LoadSceneConfig(const std::string &filePath) {
  std::ifstream input(filePath);
  if (!input) {
    throw std::runtime_error("cannot open scene file: " + filePath);
  }

  Json root;
  try {
    input >> root;
  } catch (const Json::exception &error) {
    throw std::runtime_error("invalid JSON in " + filePath + ": " +
                             error.what());
  }

  if (!root.is_object()) {
    throw std::runtime_error("scene root must be a JSON object");
  }

  try {
    const Json &source = root.at("source");
    const Json &optics = root.at("optics");
    if (!source.is_object()) {
      throw std::runtime_error("source must be a JSON object");
    }
    if (!optics.is_object()) {
      throw std::runtime_error("optics must be a JSON object");
    }
    const Json &rayCount = source.at("ray_count");
    if (!rayCount.is_number_integer()) {
      throw std::runtime_error("source.ray_count must be an integer");
    }
    const Json &maxBoundaryInteractions =
        optics.at("max_boundary_interactions");
    if (!maxBoundaryInteractions.is_number_integer()) {
      throw std::runtime_error(
          "optics.max_boundary_interactions must be an integer");
    }
    int workerCount = 0;
    if (optics.contains("worker_count")) {
      if (!optics.at("worker_count").is_number_integer()) {
        throw std::runtime_error("optics.worker_count must be an integer");
      }
      workerCount = optics.at("worker_count").get<int>();
    }

    SceneConfig scene{ReadVector3(source.at("position"), "source.position"),
                      ReadVector3(source.at("direction"), "source.direction"),
                      rayCount.get<int>(),
                      ReadFiniteNumber(source.at("opening_angle_degrees"),
                                       "source.opening_angle_degrees"),
                      ReadFiniteNumber(root.at("path_length"), "path_length"),
                      ReadFiniteNumber(optics.at("minimum_branch_intensity"),
                                       "optics.minimum_branch_intensity"),
                      maxBoundaryInteractions.get<int>(),
                      workerCount,
                      {},
                      std::nullopt};

    const Json &objects = root.at("objects");
    if (!objects.is_array()) {
      throw std::runtime_error("objects must be a JSON array");
    }
    if (objects.size() > 256) {
      throw std::runtime_error("objects may contain at most 256 spheres");
    }

    scene.spheres.reserve(objects.size());
    for (std::size_t i = 0; i < objects.size(); ++i) {
      const Json &object = objects[i];
      const std::string prefix = "objects[" + std::to_string(i) + "]";
      if (!object.is_object()) {
        throw std::runtime_error(prefix + " must be a JSON object");
      }
      if (object.at("type").get<std::string>() != "sphere") {
        throw std::runtime_error(prefix + ".type must be \"sphere\"");
      }

      const std::string name =
          object.value("name", "Sphere " + std::to_string(i + 1));
      const Json &colorValue = object.at("color");
      if (!colorValue.is_number_integer()) {
        throw std::runtime_error(prefix + ".color must be an integer ROOT "
                                          "color index");
      }
      const int color = colorValue.get<int>();
      if (color < 0 || color > std::numeric_limits<short>::max()) {
        throw std::runtime_error(
            prefix + ".color is outside the ROOT color index range");
      }

      scene.spheres.push_back(
          {name, ReadVector3(object.at("center"), prefix + ".center"),
           ReadFiniteNumber(object.at("radius"), prefix + ".radius"),
           ReadFiniteNumber(object.at("refractive_index"),
                            prefix + ".refractive_index"),
           color});
    }

    if (root.contains("detector") && !root.at("detector").is_null()) {
      const Json &detector = root.at("detector");
      if (!detector.is_object()) {
        throw std::runtime_error("detector must be a JSON object or null");
      }
      const Json &binsX = detector.at("bins_x");
      const Json &binsY = detector.at("bins_y");
      const Json &colorValue = detector.at("color");
      if (!binsX.is_number_integer() || !binsY.is_number_integer()) {
        throw std::runtime_error("detector.bins_x and detector.bins_y must "
                                 "be integers");
      }
      if (!colorValue.is_number_integer()) {
        throw std::runtime_error("detector.color must be an integer ROOT "
                                 "color index");
      }

      scene.detector = PlaneDetectorDescription{
          detector.value("name", std::string("Plane detector")),
          ReadVector3(detector.at("center"), "detector.center"),
          ReadVector3(detector.at("normal"), "detector.normal"),
          ReadVector3(detector.at("u_axis"), "detector.u_axis"),
          ReadFiniteNumber(detector.at("width"), "detector.width"),
          ReadFiniteNumber(detector.at("height"), "detector.height"),
          binsX.get<int>(),
          binsY.get<int>(),
          colorValue.get<int>()};
    }

    ValidateScene(scene);
    return scene;
  } catch (const Json::exception &error) {
    throw std::runtime_error("invalid scene in " + filePath + ": " +
                             error.what());
  }
}

std::vector<std::array<double, 3>>
GenerateRayDirections(const SceneConfig &scene) {
  constexpr double kPi = 3.14159265358979323846;
  constexpr double kGoldenAngle = 2.39996322972865332;

  std::array<double, 3> forward = scene.sourceDirection;
  const double magnitude = std::hypot(forward[0], forward[1], forward[2]);
  const double inverseMagnitude = 1.0 / magnitude;
  for (double &component : forward) {
    component *= inverseMagnitude;
  }

  const std::array<double, 3> reference =
      std::abs(forward[2]) < 0.999 ? std::array<double, 3>{0.0, 0.0, 1.0}
                                   : std::array<double, 3>{0.0, 1.0, 0.0};
  std::array<double, 3> right{
      reference[1] * forward[2] - reference[2] * forward[1],
      reference[2] * forward[0] - reference[0] * forward[2],
      reference[0] * forward[1] - reference[1] * forward[0]};
  const double rightMagnitude = std::sqrt(
      right[0] * right[0] + right[1] * right[1] + right[2] * right[2]);
  for (double &component : right) {
    component /= rightMagnitude;
  }
  const std::array<double, 3> up{forward[1] * right[2] - forward[2] * right[1],
                                 forward[2] * right[0] - forward[0] * right[2],
                                 forward[0] * right[1] - forward[1] * right[0]};

  const double halfAngle = 0.5 * scene.openingAngleDegrees * kPi / 180.0;
  const double minimumCosine = std::cos(halfAngle);
  std::vector<std::array<double, 3>> directions;
  directions.reserve(static_cast<std::size_t>(scene.rayCount));
  directions.push_back(forward);

  for (int i = 1; i < scene.rayCount; ++i) {
    const double radialFraction =
        static_cast<double>(i) / static_cast<double>(scene.rayCount - 1);
    const double cosine = 1.0 - radialFraction * (1.0 - minimumCosine);
    const double sine = std::sqrt(std::max(0.0, 1.0 - cosine * cosine));
    const double azimuth = (i - 1) * kGoldenAngle;
    const std::array<double, 3> direction{
        cosine * forward[0] +
            sine * (std::cos(azimuth) * right[0] + std::sin(azimuth) * up[0]),
        cosine * forward[1] +
            sine * (std::cos(azimuth) * right[1] + std::sin(azimuth) * up[1]),
        cosine * forward[2] +
            sine * (std::cos(azimuth) * right[2] + std::sin(azimuth) * up[2])};
    directions.push_back(direction);
  }

  return directions;
}
