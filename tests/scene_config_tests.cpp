#include "Fresnel.h"
#include "SceneConfig.h"

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void Require(bool condition, const std::string &message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

class TemporaryScene {
public:
  explicit TemporaryScene(const std::string &contents) {
    const auto now =
        std::chrono::steady_clock::now().time_since_epoch().count();
    path_ = std::filesystem::temp_directory_path() /
            ("raytracer-scene-" + std::to_string(now) + ".json");
    std::ofstream output(path_);
    output << contents;
    if (!output) {
      throw std::runtime_error("could not write temporary scene file");
    }
  }

  ~TemporaryScene() {
    std::error_code error;
    std::filesystem::remove(path_, error);
  }

  const std::filesystem::path &Path() const { return path_; }

private:
  std::filesystem::path path_;
};

bool FailsToLoad(const std::string &contents) {
  TemporaryScene file(contents);
  try {
    LoadSceneConfig(file.Path().string());
  } catch (const std::runtime_error &) {
    return true;
  }
  return false;
}

} // namespace

int main() {
  try {
    const FresnelSplit normalIncidence =
        ComputeFresnelSplit({0.0, 0.0, 1.0}, {0.0, 0.0, -1.0}, 1.0, 1.5);
    Require(normalIncidence.hasTransmission,
            "air-to-glass incidence should transmit");
    Require(std::abs(normalIncidence.reflectance - 0.04) < 1e-12,
            "normal-incidence Fresnel reflectance should be 4 percent");
    Require(std::abs(normalIncidence.transmittance +
                     normalIncidence.reflectance - 1.0) < 1e-12,
            "reflected and transmitted power fractions should sum to one");
    Require(std::abs(normalIncidence.reflectedDirection[2] + 1.0) < 1e-12 &&
                std::abs(normalIncidence.transmittedDirection[2] - 1.0) < 1e-12,
            "normal-incidence directions should be correct");

    const FresnelSplit obliqueIncidence =
        ComputeFresnelSplit({0.6, 0.0, 0.8}, {0.0, 0.0, -1.0}, 1.0, 1.5);
    Require(obliqueIncidence.hasTransmission &&
                std::abs(obliqueIncidence.transmittedDirection[0] - 0.4) <
                    1e-12,
            "oblique refraction should obey Snell's law");
    Require(std::abs(obliqueIncidence.reflectance - 0.0439) < 1e-4,
            "oblique reflectance should use the unpolarized Fresnel equations");

    const FresnelSplit totalInternalReflection = ComputeFresnelSplit(
        {std::sqrt(0.75), 0.0, 0.5}, {0.0, 0.0, -1.0}, 1.5, 1.0);
    Require(!totalInternalReflection.hasTransmission &&
                totalInternalReflection.reflectance == 1.0,
            "glass-to-air incidence above critical angle should fully reflect");
    Require(
        std::abs(totalInternalReflection.reflectedDirection[2] + 0.5) < 1e-12,
        "total internal reflection should reflect about the surface normal");

    const SceneConfig example =
        LoadSceneConfig(std::string(RAYTRACER_SOURCE_DIR) + "/scene.json");
    Require(example.spheres.size() == 2,
            "example scene should load both sphere objects");
    Require(example.spheres[1].name == "Second glass sphere",
            "object names should be loaded");
    Require(example.spheres[1].center[0] == 4.5,
            "object centers should be loaded");
    Require(example.spheres[1].refractiveIndex == 1.33,
            "object refractive indices should be loaded");
    Require(example.rayCount == 500, "ray count should be loaded");
    Require(example.openingAngleDegrees == 60.0,
            "opening angle should be loaded");
    Require(example.minimumBranchIntensity == 0.01,
            "Fresnel branch threshold should be loaded");
    Require(example.maxBoundaryInteractions == 8,
            "boundary interaction limit should be loaded");

    const std::vector<std::array<double, 3>> directions =
        GenerateRayDirections(example);
    Require(directions.size() == static_cast<std::size_t>(example.rayCount),
            "source should generate the configured number of rays");
    for (const auto &direction : directions) {
      const double magnitude =
          std::hypot(direction[0], direction[1], direction[2]);
      Require(std::abs(magnitude - 1.0) < 1e-12,
              "generated directions should be normalized");
      Require(direction[2] >=
                  std::cos(30.0 * 3.14159265358979323846 / 180.0) - 1e-12,
              "generated rays should stay inside the configured cone");
    }
    Require(std::abs(directions.front()[0]) < 1e-12 &&
                std::abs(directions.front()[1]) < 1e-12 &&
                std::abs(directions.front()[2] - 1.0) < 1e-12,
            "first generated ray should follow the configured source axis");

    SceneConfig angledSource = example;
    angledSource.sourceDirection = {2.0, -1.0, 3.0};
    angledSource.rayCount = 25;
    angledSource.openingAngleDegrees = 40.0;
    const auto angledDirections = GenerateRayDirections(angledSource);
    const double axisMagnitude = std::hypot(2.0, -1.0, 3.0);
    Require(angledDirections.size() == 25,
            "ray count should remain exact for an angled source");
    for (const auto &direction : angledDirections) {
      const double axisDot =
          (2.0 * direction[0] - direction[1] + 3.0 * direction[2]) /
          axisMagnitude;
      Require(axisDot >=
                  std::cos(20.0 * 3.14159265358979323846 / 180.0) - 1e-12,
              "angled source rays should stay inside its cone");
    }

    Require(FailsToLoad("{"), "malformed JSON should be rejected");
    Require(FailsToLoad(R"json({
      "source":{"position":[0,0,-8],"direction":[0,0,1],"ray_count":1,
                "opening_angle_degrees":0},
      "path_length":10,
      "optics":{"minimum_branch_intensity":0.01,"max_boundary_interactions":8},
      "objects":[{"type":"box","center":[0,0,0],"radius":1,
                  "refractive_index":1.5,"color":800}]
    })json"),
            "unsupported object types should be rejected");
    Require(FailsToLoad(R"json({
      "source":{"position":[0,0,-8],"direction":[0,0,1],"ray_count":1,
                "opening_angle_degrees":0},
      "path_length":10,
      "optics":{"minimum_branch_intensity":0.01,"max_boundary_interactions":8},
      "objects":[
        {"type":"sphere","center":[0,0,0],"radius":2,"refractive_index":1.5,"color":800},
        {"type":"sphere","center":[3,0,0],"radius":2,"refractive_index":1.5,"color":600}
      ]
    })json"),
            "overlapping spheres should be rejected");
    Require(FailsToLoad(R"json({
      "source":{"position":[0,0,-8],"direction":[0,0,1],"ray_count":1,
                "opening_angle_degrees":0},
      "path_length":10,
      "optics":{"minimum_branch_intensity":0.01,"max_boundary_interactions":8},
      "objects":[{"type":"sphere","center":[0,0,0],"radius":1,
                  "refractive_index":1.5,"color":999999}]
    })json"),
            "invalid ROOT color indices should be rejected");
    Require(FailsToLoad(R"json({
      "source":{"position":[0,0,-8],"direction":[0,0,1],"ray_count":0,
                "opening_angle_degrees":0},
      "path_length":10,
      "optics":{"minimum_branch_intensity":0.01,"max_boundary_interactions":8},
      "objects":[{"type":"sphere","center":[0,0,0],"radius":1,
                  "refractive_index":1.5,"color":800}]
    })json"),
            "zero ray count should be rejected");
    Require(FailsToLoad(R"json({
      "source":{"position":[0,0,-8],"direction":[0,0,0],"ray_count":1,
                "opening_angle_degrees":0},
      "path_length":10,
      "optics":{"minimum_branch_intensity":0.01,"max_boundary_interactions":8},
      "objects":[{"type":"sphere","center":[0,0,0],"radius":1,
                  "refractive_index":1.5,"color":800}]
    })json"),
            "zero source direction should be rejected");
    Require(FailsToLoad(R"json({
      "source":{"position":[0,0,-8],"direction":[0,0,1],"ray_count":1,
                "opening_angle_degrees":181},
      "path_length":10,
      "optics":{"minimum_branch_intensity":0.01,"max_boundary_interactions":8},
      "objects":[{"type":"sphere","center":[0,0,0],"radius":1,
                  "refractive_index":1.5,"color":800}]
    })json"),
            "opening angles above 180 degrees should be rejected");
    Require(FailsToLoad(R"json({
      "source":{"position":[0,0,-8],"direction":[0,0,1],"ray_count":1,
                "opening_angle_degrees":0},
      "path_length":10,
      "optics":{"minimum_branch_intensity":1,"max_boundary_interactions":8},
      "objects":[{"type":"sphere","center":[0,0,0],"radius":1,
                  "refractive_index":1.5,"color":800}]
    })json"),
            "branch intensity threshold at 1 should be rejected");

    std::cout << "Scene config tests passed\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "Scene config test failed: " << error.what() << '\n';
    return 1;
  }
}
