cmake_minimum_required(VERSION 3.16)

project(ROOTRayTracer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(ROOT CONFIG REQUIRED COMPONENTS Core Graf Graf3d Gpad Hist Physics)
find_package(nlohmann_json 3.2.0 CONFIG REQUIRED)
find_package(Threads REQUIRED)
find_package(OpenCASCADE REQUIRED)

set(STEP_GEOMETRY_LIBRARIES
  TKDESTEP
  TKXSBase
  TKDE
  TKMesh
  TKPrim
  TKBRep
  TKTopAlgo
  TKGeomAlgo
  TKGeomBase
  TKG3d
  TKG2d
  TKMath
  TKernel
)

add_executable(raytracer
  raytracer.C
  raytracer_main.cpp
  SceneConfig.cpp
  Fresnel.cpp
  Detector.cpp
  RayTracingEngine.cpp
  StepGeometry.cpp
)

target_include_directories(raytracer PRIVATE ${OpenCASCADE_INCLUDE_DIR})
target_link_libraries(raytracer PRIVATE
  ROOT::Core
  ROOT::Graf
  ROOT::Graf3d
  ROOT::Gpad
  ROOT::Hist
  ROOT::Physics
  nlohmann_json::nlohmann_json
  Threads::Threads
  ${STEP_GEOMETRY_LIBRARIES}
)

configure_file(scene.json scene.json COPYONLY)

include(CTest)
if(BUILD_TESTING)
  add_executable(scene_config_tests
    SceneConfig.cpp
    Fresnel.cpp
    Detector.cpp
    RayTracingEngine.cpp
    StepGeometry.cpp
    tests/scene_config_tests.cpp
  )
  target_include_directories(scene_config_tests PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    ${OpenCASCADE_INCLUDE_DIR}
  )
  target_link_libraries(scene_config_tests PRIVATE
    nlohmann_json::nlohmann_json
    Threads::Threads
    ${STEP_GEOMETRY_LIBRARIES}
  )
  target_compile_definitions(scene_config_tests PRIVATE
    RAYTRACER_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
  )
  add_test(NAME scene_config_tests COMMAND scene_config_tests)
endif()
