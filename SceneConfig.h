#ifndef ROOT_RAY_TRACER_SCENE_CONFIG_H
#define ROOT_RAY_TRACER_SCENE_CONFIG_H

#include <array>
#include <optional>
#include <string>
#include <vector>

enum class ObjectShape { Sphere, Box, Tube, Cone };

struct ObjectDescription {
  std::string name;
  ObjectShape shape;
  std::array<double, 3> center;
  std::array<double, 3> size;
  double radius;
  double topRadius;
  double refractiveIndex;
  int color;
};

struct PlaneDetectorDescription {
  std::string name;
  std::array<double, 3> center;
  std::array<double, 3> normal;
  std::array<double, 3> uAxis;
  double width;
  double height;
  int binsX;
  int binsY;
  int color;
};

struct SceneConfig {
  std::array<double, 3> sourcePosition;
  std::array<double, 3> sourceDirection;
  int rayCount;
  double openingAngleDegrees;
  double pathLength;
  double minimumBranchIntensity;
  int maxBoundaryInteractions;
  int workerCount;
  std::vector<ObjectDescription> objects;
  std::optional<PlaneDetectorDescription> detector;
};

SceneConfig LoadSceneConfig(const std::string &filePath);
std::vector<std::array<double, 3>>
GenerateRayDirections(const SceneConfig &scene);

#endif
