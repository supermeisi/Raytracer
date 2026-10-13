#include "StepGeometry.h"

#include <BRep_Tool.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <Poly_Triangulation.hxx>
#include <STEPControl_Reader.hxx>
#include <TopAbs_Orientation.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopLoc_Location.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>

#include <algorithm>
#include <stdexcept>

std::vector<MeshTriangle> LoadStepTriangles(const std::string &filePath,
                                            double linearDeflection) {
  if (linearDeflection <= 0.0) {
    throw std::invalid_argument("STEP mesh deflection must be greater than 0");
  }

  STEPControl_Reader reader;
  if (reader.ReadFile(filePath.c_str()) != IFSelect_RetDone) {
    throw std::runtime_error("cannot read STEP file: " + filePath);
  }
  if (reader.TransferRoots() == 0) {
    throw std::runtime_error("STEP file contains no transferable roots: " +
                             filePath);
  }
  const TopoDS_Shape shape = reader.OneShape();
  if (shape.IsNull()) {
    throw std::runtime_error("STEP file contains no geometry: " + filePath);
  }

  BRepMesh_IncrementalMesh mesher(shape, linearDeflection);
  if (!mesher.IsDone()) {
    throw std::runtime_error("could not triangulate STEP file: " + filePath);
  }

  std::vector<MeshTriangle> triangles;
  for (TopExp_Explorer faces(shape, TopAbs_FACE); faces.More(); faces.Next()) {
    const TopoDS_Face face = TopoDS::Face(faces.Current());
    TopLoc_Location location;
    const Handle(Poly_Triangulation) triangulation =
        BRep_Tool::Triangulation(face, location);
    if (triangulation.IsNull()) {
      continue;
    }
    const gp_Trsf transform = location.Transformation();
    const bool reversed = face.Orientation() == TopAbs_REVERSED;
    for (Standard_Integer index = 1; index <= triangulation->NbTriangles();
         ++index) {
      Standard_Integer first = 0;
      Standard_Integer second = 0;
      Standard_Integer third = 0;
      triangulation->Triangle(index).Get(first, second, third);
      if (reversed) {
        std::swap(second, third);
      }
      const gp_Pnt pointA = triangulation->Node(first).Transformed(transform);
      const gp_Pnt pointB = triangulation->Node(second).Transformed(transform);
      const gp_Pnt pointC = triangulation->Node(third).Transformed(transform);
      triangles.push_back({{pointA.X(), pointA.Y(), pointA.Z()},
                           {pointB.X(), pointB.Y(), pointB.Z()},
                           {pointC.X(), pointC.Y(), pointC.Z()}});
    }
  }
  if (triangles.empty()) {
    throw std::runtime_error("STEP file produced no triangles: " + filePath);
  }
  return triangles;
}
