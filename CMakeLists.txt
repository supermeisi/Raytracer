cmake_minimum_required(VERSION 3.16)

project(ROOTRayTracer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(ROOT CONFIG REQUIRED COMPONENTS Core Graf Graf3d Gpad Physics)

add_executable(raytracer
  raytracer.C
  raytracer_main.cpp
)

target_link_libraries(raytracer PRIVATE
  ROOT::Core
  ROOT::Graf
  ROOT::Graf3d
  ROOT::Gpad
  ROOT::Physics
)
