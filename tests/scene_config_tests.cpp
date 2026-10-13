#include "SceneConfig.h"

#include <chrono>
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
    const auto now = std::chrono::steady_clock::now().time_since_epoch().count();
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

    Require(FailsToLoad("{"), "malformed JSON should be rejected");
    Require(FailsToLoad(R"json({
      "source":{"position":[0,0,-8],"grid":{"radius":0,"target_z":-3,"spacing":1}},
      "path_length":10,
      "objects":[{"type":"box","center":[0,0,0],"radius":1,
                  "refractive_index":1.5,"color":800}]
    })json"),
            "unsupported object types should be rejected");
    Require(FailsToLoad(R"json({
      "source":{"position":[0,0,-8],"grid":{"radius":0,"target_z":-3,"spacing":1}},
      "path_length":10,
      "objects":[
        {"type":"sphere","center":[0,0,0],"radius":2,"refractive_index":1.5,"color":800},
        {"type":"sphere","center":[3,0,0],"radius":2,"refractive_index":1.5,"color":600}
      ]
    })json"),
            "overlapping spheres should be rejected");
    Require(FailsToLoad(R"json({
      "source":{"position":[0,0,-8],"grid":{"radius":0,"target_z":-3,"spacing":1}},
      "path_length":10,
      "objects":[{"type":"sphere","center":[0,0,0],"radius":1,
                  "refractive_index":1.5,"color":999999}]
    })json"),
            "invalid ROOT color indices should be rejected");

    std::cout << "Scene config tests passed\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "Scene config test failed: " << error.what() << '\n';
    return 1;
  }
}
