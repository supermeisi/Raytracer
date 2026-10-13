#include "SceneConfig.h"

#include <nlohmann/json.hpp>

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
  if (scene.gridRadius < 0 || scene.gridRadius > 100) {
    throw std::runtime_error("source.grid.radius must be between 0 and 100");
  }
  if (scene.targetSpacing <= 0.0) {
    throw std::runtime_error("source.grid.spacing must be greater than 0");
  }
  if (scene.pathLength <= 0.0) {
    throw std::runtime_error("path_length must be greater than 0");
  }
  if (scene.spheres.empty()) {
    throw std::runtime_error("objects must contain at least one sphere");
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
      throw std::runtime_error(
          "source.position must be outside every sphere");
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
    const Json &grid = source.at("grid");
    if (!source.is_object() || !grid.is_object()) {
      throw std::runtime_error("source and source.grid must be JSON objects");
    }
    const Json &gridRadius = grid.at("radius");
    if (!gridRadius.is_number_integer()) {
      throw std::runtime_error("source.grid.radius must be an integer");
    }

    SceneConfig scene{
        ReadVector3(source.at("position"), "source.position"),
        gridRadius.get<int>(),
        ReadFiniteNumber(grid.at("target_z"), "source.grid.target_z"),
        ReadFiniteNumber(grid.at("spacing"), "source.grid.spacing"),
        ReadFiniteNumber(root.at("path_length"), "path_length"),
        {}};

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
        throw std::runtime_error(prefix +
                                 ".color is outside the ROOT color index range");
      }

      scene.spheres.push_back(
          {name,
           ReadVector3(object.at("center"), prefix + ".center"),
           ReadFiniteNumber(object.at("radius"), prefix + ".radius"),
           ReadFiniteNumber(object.at("refractive_index"),
                            prefix + ".refractive_index"),
           color});
    }

    ValidateScene(scene);
    return scene;
  } catch (const Json::exception &error) {
    throw std::runtime_error("invalid scene in " + filePath + ": " +
                             error.what());
  }
}
