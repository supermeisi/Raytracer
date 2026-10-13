#include <TApplication.h>

#include "SceneConfig.h"

#include <exception>
#include <iostream>

void raytracer(const SceneConfig &scene);

int main(int argc, char **argv) {
  if (argc > 2) {
    std::cerr << "Usage: " << argv[0] << " [scene.json]\n";
    return 2;
  }

  try {
    const SceneConfig scene =
        LoadSceneConfig(argc == 2 ? argv[1] : "scene.json");
    int rootArgc = 1;
    char *rootArgv[] = {argv[0], nullptr};
    TApplication application("raytracer", &rootArgc, rootArgv);
    raytracer(scene);
    application.Run();
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "Ray tracer error: " << error.what() << '\n';
    return 1;
  }
}
