#include "Detector.h"
#include "Fresnel.h"
#include "RayTracingEngine.h"
#include "SceneConfig.h"

#include <BRepPrimAPI_MakeBox.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <STEPControl_StepModelType.hxx>
#include <STEPControl_Writer.hxx>

#include <atomic>
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

class TemporaryStep {
public:
  TemporaryStep() {
    const auto now =
        std::chrono::steady_clock::now().time_since_epoch().count();
    path_ = std::filesystem::temp_directory_path() /
            ("raytracer-step-" + std::to_string(now) + ".step");
    STEPControl_Writer writer;
    if (writer.Transfer(BRepPrimAPI_MakeBox(2.0, 2.0, 2.0).Shape(),
                        STEPControl_AsIs) != IFSelect_RetDone ||
        writer.Write(path_.string().c_str()) != IFSelect_RetDone) {
      throw std::runtime_error("could not create temporary STEP test file");
    }
  }

  ~TemporaryStep() {
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
    Require(example.objects.size() == 5,
            "example scene should load all five optical objects");
    Require(example.objects[1].name == "Second glass sphere",
            "object names should be loaded");
    Require(example.objects[1].center[0] == 4.5,
            "object centers should be loaded");
    Require(example.objects[1].refractiveIndex == 1.33,
            "object refractive indices should be loaded");
    Require(example.objects[0].shape == ObjectShape::Sphere &&
                example.objects[0].radius == 3.0,
            "sphere shape should be loaded");
    Require(example.rayCount == 500, "ray count should be loaded");
    Require(example.openingAngleDegrees == 60.0,
            "opening angle should be loaded");
    Require(example.minimumBranchIntensity == 0.01,
            "Fresnel branch threshold should be loaded");
    Require(example.maxBoundaryInteractions == 8,
            "boundary interaction limit should be loaded");
    Require(example.workerCount == 4,
            "example scene worker count should be loaded");
    Require(example.detector.has_value(),
            "example scene should configure a plane detector");
    Require(example.detector->binsX == 100 && example.detector->binsY == 100,
            "detector histogram dimensions should be loaded");

    PlaneDetectorHistogram detector(*example.detector);
    const auto detectorHit =
        detector.Intersect({0.0, 0.0, -8.0}, {0.0, 0.0, 1.0}, 30.0);
    Require(detectorHit.has_value(),
            "ray along the detector normal should register a hit");
    Require(std::abs(detectorHit->distance - 16.0) < 1e-12 &&
                std::abs(detectorHit->u) < 1e-12 &&
                std::abs(detectorHit->v) < 1e-12,
            "detector hit should have expected distance and local coordinates");
    detector.Record(*detectorHit, 0.25);
    detector.Record(*detectorHit, 0.5);
    Require(std::abs(detector.TotalIntensity() - 0.75) < 1e-12,
            "detector should accumulate ray intensity");
    Require(std::abs(detector.Bins()[50 * 100 + 50] - 0.75) < 1e-12,
            "center hit should land in the center histogram bin");
    Require(!detector.Intersect({0.0, 0.0, -8.0}, {1.0, 0.0, 0.0}, 30.0),
            "ray parallel to the detector plane should not hit it");
    Require(!detector.Intersect({30.0, 0.0, -8.0}, {0.0, 0.0, 1.0}, 30.0),
            "ray outside detector bounds should not hit it");
    Require(!detector.Intersect({0.0, 0.0, -8.0}, {0.0, 0.0, 1.0}, 10.0),
            "detector beyond the remaining ray length should not hit");
    const auto nonUnitDirectionHit =
        detector.Intersect({0.0, 0.0, -8.0}, {0.0, 0.0, 2.0}, 16.0);
    Require(nonUnitDirectionHit.has_value() &&
                std::abs(nonUnitDirectionHit->distance - 16.0) < 1e-12,
            "detector intersection should use physical distance");

    TemporaryScene detectorDisabledScene(R"json({
      "source":{"position":[0,0,-8],"direction":[0,0,1],"ray_count":1,
                "opening_angle_degrees":0},
      "path_length":10,
      "optics":{"minimum_branch_intensity":0.01,"max_boundary_interactions":8},
      "objects":[{"type":"sphere","center":[0,0,0],"radius":1,
                  "refractive_index":1.5,"color":800}]
    })json");
    Require(!LoadSceneConfig(detectorDisabledScene.Path().string())
                 .detector.has_value(),
            "scene should allow the detector to be omitted");

    TemporaryScene additionalShapesScene(R"json({
      "source":{"position":[0,0,-20],"direction":[0,0,1],"ray_count":1,
                "opening_angle_degrees":0},
      "path_length":40,
      "optics":{"minimum_branch_intensity":0.01,"max_boundary_interactions":8},
      "objects":[
        {"type":"sphere","center":[0,0,0],"radius":1,
         "refractive_index":1.5,"color":800},
        {"type":"box","center":[10,0,0],"size":[2,4,6],
         "refractive_index":1.4,"color":600},
        {"type":"tube","center":[20,0,0],"radius":2,"height":6,
         "refractive_index":1.3,"color":400},
        {"type":"cone","center":[30,0,0],"bottom_radius":2,"top_radius":0.5,
         "height":6,"refractive_index":1.6,"color":2}
      ]
    })json");
    const SceneConfig additionalShapes =
        LoadSceneConfig(additionalShapesScene.Path().string());
    Require(additionalShapes.objects.size() == 4,
            "sphere, box, tube, and cone should load together");
    Require(additionalShapes.objects[1].shape == ObjectShape::Box &&
                additionalShapes.objects[1].size ==
                    std::array<double, 3>{2.0, 4.0, 6.0},
            "box dimensions should be loaded");
    Require(additionalShapes.objects[2].shape == ObjectShape::Tube &&
                additionalShapes.objects[2].radius == 2.0 &&
                additionalShapes.objects[2].size[2] == 6.0,
            "tube radius and height should be loaded");
    Require(additionalShapes.objects[3].shape == ObjectShape::Cone &&
                additionalShapes.objects[3].radius == 2.0 &&
                additionalShapes.objects[3].topRadius == 0.5,
            "cone bottom and top radii should be loaded");

    TemporaryStep stepFile;
    const std::string stepSceneJson =
        std::string(R"json({
      "source":{"position":[1,1,-5],"direction":[0,0,1],"ray_count":1,
                "opening_angle_degrees":0},
      "path_length":10,
      "optics":{"minimum_branch_intensity":0.01,"max_boundary_interactions":8},
      "objects":[{"type":"step","file":")json") +
        stepFile.Path().filename().string() +
        R"json(","center":[0,0,0],"mesh_deflection":0.1,
                  "refractive_index":1.5,"color":800}]
    })json";
    TemporaryScene stepSceneFile(stepSceneJson);
    const SceneConfig stepScene = LoadSceneConfig(stepSceneFile.Path().string());
    Require(stepScene.objects.size() == 1 &&
                stepScene.objects[0].shape == ObjectShape::Step &&
                !stepScene.objects[0].triangles.empty(),
            "STEP geometry should load and tessellate");
    std::atomic<std::size_t> stepProgress{0};
    const std::vector<std::array<double, 3>> stepRayDirection{{0, 0, 1}};
    const auto stepRays = TraceRaysParallel(
        stepScene, stepRayDirection, 1, stepProgress);
    bool stepObjectTraversed = false;
    for (const TracedRaySegment &segment : stepRays[0].segments) {
      stepObjectTraversed =
          stepObjectTraversed || segment.objectIndex == 0;
    }
    Require(stepObjectTraversed && !stepRays[0].segments.empty() &&
                std::abs(stepRays[0].segments[0].end[2]) < 1e-9,
            "ray should intersect the imported STEP box");
    Require(FailsToLoad(R"json({
      "source":{"position":[0,0,-5],"direction":[0,0,1],"ray_count":1,
                "opening_angle_degrees":0},
      "path_length":10,
      "optics":{"minimum_branch_intensity":0.01,"max_boundary_interactions":8},
      "objects":[{"type":"step","file":"missing.step","center":[0,0,0],
                  "refractive_index":1.5,"color":800}]
    })json"),
            "missing STEP files should report a load failure");

    for (const ObjectDescription &object : additionalShapes.objects) {
      if (object.shape == ObjectShape::Sphere) {
        continue;
      }
      SceneConfig shapeScene = additionalShapes;
      shapeScene.objects = {object};
      shapeScene.objects[0].center = {0.0, 0.0, 0.0};
      shapeScene.sourcePosition = {0.0, 0.0, -5.0};
      shapeScene.pathLength = 10.0;
      std::atomic<std::size_t> shapeProgress{0};
      const std::vector<std::array<double, 3>> straightRay{{0.0, 0.0, 1.0}};
      const auto shapeRays =
          TraceRaysParallel(shapeScene, straightRay, 1, shapeProgress);
      bool hitOpticalObject = false;
      for (const auto &segment : shapeRays.front().segments) {
        hitOpticalObject = hitOpticalObject || segment.objectIndex == 0;
      }
      Require(hitOpticalObject,
              "box, tube, and cone should participate in ray tracing");

      shapeScene.sourcePosition = {-5.0, 0.0, 0.0};
      const std::vector<std::array<double, 3>> sideRay{{1.0, 0.0, 0.0}};
      const auto sideRays =
          TraceRaysParallel(shapeScene, sideRay, 1, shapeProgress);
      bool hitSide = false;
      for (const auto &segment : sideRays.front().segments) {
        hitSide = hitSide || segment.objectIndex == 0;
      }
      Require(hitSide,
              "box, tube, and cone side surfaces should intersect rays");
    }

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

    SceneConfig threadedScene = example;
    threadedScene.rayCount = 64;
    threadedScene.openingAngleDegrees = 30.0;
    const auto threadedDirections = GenerateRayDirections(threadedScene);
    std::atomic<std::size_t> serialCompleted{0};
    const auto serialRays =
        TraceRaysParallel(threadedScene, threadedDirections, 1, serialCompleted);
    Require(serialCompleted.load() == threadedDirections.size(),
            "serial tracing progress should reach the ray count");
    std::atomic<std::size_t> parallelCompleted{0};
    const auto parallelRays = TraceRaysParallel(
        threadedScene, threadedDirections, 4, parallelCompleted);
    Require(parallelCompleted.load() == threadedDirections.size(),
            "parallel tracing progress should reach the ray count");
    Require(serialRays.size() == parallelRays.size(),
            "serial and parallel tracing should produce the same ray count");
    for (std::size_t ray = 0; ray < serialRays.size(); ++ray) {
      Require(serialRays[ray].segments.size() ==
                  parallelRays[ray].segments.size(),
              "serial and parallel tracing should produce the same segments");
      for (std::size_t segment = 0;
           segment < serialRays[ray].segments.size(); ++segment) {
        const auto &serial = serialRays[ray].segments[segment];
        const auto &parallel = parallelRays[ray].segments[segment];
        Require(serial.start == parallel.start && serial.end == parallel.end &&
                    serial.objectIndex == parallel.objectIndex &&
                    serial.intensity == parallel.intensity,
                "thread count should not change traced ray results");
      }
      Require(serialRays[ray].detectorHit.has_value() ==
                  parallelRays[ray].detectorHit.has_value(),
              "thread count should not change detector hit results");
      if (serialRays[ray].detectorHit) {
        Require(serialRays[ray].detectorHit->binX ==
                        parallelRays[ray].detectorHit->binX &&
                    serialRays[ray].detectorHit->binY ==
                        parallelRays[ray].detectorHit->binY &&
                    serialRays[ray].detectorIntensity ==
                        parallelRays[ray].detectorIntensity,
                "thread count should not change detector hit bins");
      }
    }

    PlaneDetectorHistogram serialHistogram(*threadedScene.detector);
    PlaneDetectorHistogram parallelHistogram(*threadedScene.detector);
    for (const TracedRay &ray : serialRays) {
      if (ray.detectorHit) {
        serialHistogram.Record(*ray.detectorHit, ray.detectorIntensity);
      }
    }
    for (const TracedRay &ray : parallelRays) {
      if (ray.detectorHit) {
        parallelHistogram.Record(*ray.detectorHit, ray.detectorIntensity);
      }
    }
    Require(serialHistogram.Bins() == parallelHistogram.Bins() &&
                serialHistogram.TotalIntensity() ==
                    parallelHistogram.TotalIntensity(),
            "serial and parallel detector histograms should match");

    Require(FailsToLoad("{"), "malformed JSON should be rejected");
    Require(FailsToLoad(R"json({
      "source":{"position":[0,0,-8],"direction":[0,0,1],"ray_count":1,
                "opening_angle_degrees":0},
      "path_length":10,
      "optics":{"minimum_branch_intensity":0.01,"max_boundary_interactions":8},
      "objects":[{"type":"box","center":[0,0,0],"size":[2,0,2],
                  "refractive_index":1.5,"color":800}]
    })json"),
            "box dimensions must be positive");
    Require(FailsToLoad(R"json({
      "source":{"position":[0,0,-8],"direction":[0,0,1],"ray_count":1,
                "opening_angle_degrees":0},
      "path_length":10,
      "optics":{"minimum_branch_intensity":0.01,"max_boundary_interactions":8,
                "worker_count":-1},
      "objects":[{"type":"sphere","center":[0,0,0],"radius":1,
                  "refractive_index":1.5,"color":800}]
    })json"),
            "negative worker count should be rejected");
    Require(FailsToLoad(R"json({
      "source":{"position":[0,0,-8],"direction":[0,0,1],"ray_count":1,
                "opening_angle_degrees":0},
      "path_length":10,
      "optics":{"minimum_branch_intensity":0.01,"max_boundary_interactions":8,
                "worker_count":257},
      "objects":[{"type":"sphere","center":[0,0,0],"radius":1,
                  "refractive_index":1.5,"color":800}]
    })json"),
            "worker counts above 256 should be rejected");
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
      "optics":{"minimum_branch_intensity":0.01,"max_boundary_interactions":8},
      "detector":{"center":[0,0,8],"normal":[0,0,1],"u_axis":[0,0,2],
                  "width":10,"height":10,"bins_x":10,"bins_y":10,"color":2},
      "objects":[{"type":"sphere","center":[0,0,0],"radius":1,
                  "refractive_index":1.5,"color":800}]
    })json"),
            "detector axis parallel to its normal should be rejected");
    Require(FailsToLoad(R"json({
      "source":{"position":[0,0,-8],"direction":[0,0,1],"ray_count":1,
                "opening_angle_degrees":0},
      "path_length":10,
      "optics":{"minimum_branch_intensity":0.01,"max_boundary_interactions":8},
      "detector":{"center":[0,0,8],"normal":[0,0,1],"u_axis":[1,0,0],
                  "width":10,"height":10,"bins_x":1001,"bins_y":10,"color":2},
      "objects":[{"type":"sphere","center":[0,0,0],"radius":1,
                  "refractive_index":1.5,"color":800}]
    })json"),
            "excessive detector bin count should be rejected");
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
