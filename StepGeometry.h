#ifndef ROOT_RAY_TRACER_STEP_GEOMETRY_H
#define ROOT_RAY_TRACER_STEP_GEOMETRY_H

#include "SceneConfig.h"

#include <string>
#include <vector>

std::vector<MeshTriangle> LoadStepTriangles(const std::string &filePath,
                                            double linearDeflection);

#endif
