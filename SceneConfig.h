#ifndef ROOT_RAY_TRACER_SCENE_CONFIG_H
#define ROOT_RAY_TRACER_SCENE_CONFIG_H

#include <array>
#include <string>
#include <vector>

struct SphereDescription {
  std::string name;
  std::array<double, 3> center;
  double radius;
  double refractiveIndex;
  int color;
};

struct SceneConfig {
  std::array<double, 3> sourcePosition;
  int gridRadius;
  double targetZ;
  double targetSpacing;
  double pathLength;
  std::vector<SphereDescription> spheres;
};

SceneConfig LoadSceneConfig(const std::string &filePath);

#endif
