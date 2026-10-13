#include "SceneConfig.h"

#include <TCanvas.h>
#include <TColor.h>
#include <TLegend.h>
#include <TMath.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TView.h>
#include <TVector3.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace {

constexpr double kSurfaceEpsilon = 1e-5;

struct RaySegment {
  TVector3 start;
  TVector3 end;
  int sphereIndex;
};

using RayPath = std::vector<RaySegment>;

struct SceneBounds {
  TVector3 minimum{std::numeric_limits<double>::infinity(),
                   std::numeric_limits<double>::infinity(),
                   std::numeric_limits<double>::infinity()};
  TVector3 maximum{-std::numeric_limits<double>::infinity(),
                   -std::numeric_limits<double>::infinity(),
                   -std::numeric_limits<double>::infinity()};

  void Include(const TVector3 &point) {
    minimum.SetXYZ(std::min(minimum.X(), point.X()),
                   std::min(minimum.Y(), point.Y()),
                   std::min(minimum.Z(), point.Z()));
    maximum.SetXYZ(std::max(maximum.X(), point.X()),
                   std::max(maximum.Y(), point.Y()),
                   std::max(maximum.Z(), point.Z()));
  }

  void IncludeSphere(const TVector3 &center, double radius) {
    Include(center - TVector3(radius, radius, radius));
    Include(center + TVector3(radius, radius, radius));
  }
};

bool IntersectSphere(const TVector3 &origin, const TVector3 &direction,
                     const TVector3 &center, double radius, double &distance,
                     TVector3 &normal) {
  const TVector3 offset = origin - center;
  const double halfB = offset.Dot(direction);
  const double c = offset.Mag2() - radius * radius;
  const double discriminant = halfB * halfB - c;
  if (discriminant < 0.0) {
    return false;
  }

  const double root = std::sqrt(discriminant);
  double nearDistance = -halfB - root;
  double farDistance = -halfB + root;
  distance = nearDistance > kSurfaceEpsilon ? nearDistance : farDistance;
  if (distance <= kSurfaceEpsilon) {
    return false;
  }

  normal = (origin + direction * distance - center) * (1.0 / radius);
  return true;
}

TVector3 ToVector3(const std::array<double, 3> &vector) {
  return TVector3(vector[0], vector[1], vector[2]);
}

RayPath TraceRay(const TVector3 &source, TVector3 direction,
                 const SceneConfig &scene) {
  direction = direction.Unit();
  TVector3 origin = source;
  RayPath path;
  double travelled = 0.0;
  int currentSphere = -1;

  while (travelled < scene.pathLength) {
    double hitDistance = 0.0;
    TVector3 outwardNormal;
    const double remaining = scene.pathLength - travelled;
    int hitSphere = -1;
    if (currentSphere >= 0) {
      const SphereDescription &sphere = scene.spheres[currentSphere];
      if (IntersectSphere(origin, direction, ToVector3(sphere.center),
                          sphere.radius, hitDistance, outwardNormal)) {
        hitSphere = currentSphere;
      }
    } else {
      hitDistance = std::numeric_limits<double>::infinity();
      for (std::size_t i = 0; i < scene.spheres.size(); ++i) {
        double candidateDistance = 0.0;
        TVector3 candidateNormal;
        const SphereDescription &sphere = scene.spheres[i];
        if (IntersectSphere(origin, direction, ToVector3(sphere.center),
                            sphere.radius, candidateDistance,
                            candidateNormal) &&
            candidateDistance < hitDistance) {
          hitDistance = candidateDistance;
          outwardNormal = candidateNormal;
          hitSphere = static_cast<int>(i);
        }
      }
    }

    if (hitSphere < 0 || hitDistance >= remaining) {
      path.push_back({origin, origin + direction * remaining, currentSphere});
      break;
    }

    const TVector3 hit = origin + direction * hitDistance;
    path.push_back({origin, hit, currentSphere});
    travelled += hitDistance;

    const bool entering = currentSphere < 0;
    const TVector3 normalAgainstRay =
        entering ? outwardNormal : -outwardNormal;
    const int nextSphere = entering ? hitSphere : -1;
    const double incidentIndex =
        entering ? 1.0 : scene.spheres[currentSphere].refractiveIndex;
    const double transmittedIndex =
        entering ? scene.spheres[nextSphere].refractiveIndex : 1.0;
    const double eta = incidentIndex / transmittedIndex;
    const double cosIncident =
        std::max(0.0, std::min(1.0, -direction.Dot(normalAgainstRay)));
    const double sinTransmittedSquared =
        eta * eta * (1.0 - cosIncident * cosIncident);

    if (sinTransmittedSquared > 1.0) {
      direction =
          (direction + normalAgainstRay * (2.0 * cosIncident)).Unit();
    } else {
      const double cosTransmitted = std::sqrt(1.0 - sinTransmittedSquared);
      direction =
          (direction * eta +
           normalAgainstRay * (eta * cosIncident - cosTransmitted))
              .Unit();
      currentSphere = nextSphere;
    }

    origin = hit + direction * kSurfaceEpsilon;
    travelled += kSurfaceEpsilon;
  }

  return path;
}

void DrawSphereWireframe(const TVector3 &center, double radius,
                         Color_t color) {
  constexpr int kSamples = 80;
  constexpr int kLatitudeCount = 12;
  constexpr int kLongitudeCount = 16;

  for (int latitude = 1; latitude < kLatitudeCount; ++latitude) {
    const double polar =
        TMath::Pi() * latitude / static_cast<double>(kLatitudeCount);
    TPolyLine3D *line = new TPolyLine3D(kSamples + 1);
    for (int i = 0; i <= kSamples; ++i) {
      const double azimuth = 2.0 * TMath::Pi() * i / kSamples;
      line->SetPoint(
          i, center.X() + radius * std::sin(polar) * std::cos(azimuth),
          center.Y() + radius * std::sin(polar) * std::sin(azimuth),
          center.Z() + radius * std::cos(polar));
    }
    line->SetLineColor(color);
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->Draw("same");
  }

  for (int longitude = 0; longitude < kLongitudeCount; ++longitude) {
    const double azimuth =
        2.0 * TMath::Pi() * longitude / kLongitudeCount;
    TPolyLine3D *line = new TPolyLine3D(kSamples + 1);
    for (int i = 0; i <= kSamples; ++i) {
      const double polar = TMath::Pi() * i / kSamples;
      line->SetPoint(
          i, center.X() + radius * std::sin(polar) * std::cos(azimuth),
          center.Y() + radius * std::sin(polar) * std::sin(azimuth),
          center.Z() + radius * std::cos(polar));
    }
    line->SetLineColor(color);
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->Draw("same");
  }
}

void DrawSegment(const RaySegment &segment, Color_t outsideColor,
                 const SceneConfig &scene, SceneBounds &bounds) {
  TPolyLine3D *line = new TPolyLine3D(2);
  line->SetPoint(0, segment.start.X(), segment.start.Y(), segment.start.Z());
  line->SetPoint(1, segment.end.X(), segment.end.Y(), segment.end.Z());
  const Color_t color =
      segment.sphereIndex >= 0 ? scene.spheres[segment.sphereIndex].color
                               : outsideColor;
  line->SetLineColor(color);
  line->SetLineWidth(2);
  line->Draw("same");
  bounds.Include(segment.start);
  bounds.Include(segment.end);
}

TPolyLine3D *MakeLegendLine(Color_t color) {
  TPolyLine3D *line = new TPolyLine3D(2);
  line->SetPoint(0, 0.0, 0.0, 0.0);
  line->SetPoint(1, 1.0, 1.0, 1.0);
  line->SetLineColor(color);
  line->SetLineWidth(2);
  return line;
}

} // namespace

void raytracer(const SceneConfig &scene) {
  gROOT->SetStyle("Plain");
  gStyle->SetCanvasColor(kWhite);

  TCanvas *canvas =
      new TCanvas("raytracer_canvas", "ROOT 3D Ray Tracer", 1000, 800);
  TView *view = TView::CreateView(1);
  const TVector3 source = ToVector3(scene.sourcePosition);
  SceneBounds bounds;
  bounds.Include(source);

  for (const SphereDescription &sphere : scene.spheres) {
    const TVector3 center = ToVector3(sphere.center);
    DrawSphereWireframe(center, sphere.radius, sphere.color);
    bounds.IncludeSphere(center, sphere.radius);
  }

  const int rayColors[] = {kAzure + 2, kBlue + 1, kCyan + 2, kMagenta + 1,
                           kGreen + 2, kViolet + 1, kRed + 1};

  for (int y = -scene.gridRadius; y <= scene.gridRadius; ++y) {
    for (int x = -scene.gridRadius; x <= scene.gridRadius; ++x) {
      const TVector3 target(x * scene.targetSpacing, y * scene.targetSpacing,
                            scene.targetZ);
      const RayPath path = TraceRay(source, target - source, scene);
      const Color_t color =
          rayColors[(x + y + 2 * scene.gridRadius) % 7];
      for (const RaySegment &segment : path) {
        DrawSegment(segment, color, scene, bounds);
      }
    }
  }

  TPolyMarker3D *sourceMarker = new TPolyMarker3D(1);
  sourceMarker->SetPoint(0, source.X(), source.Y(), source.Z());
  sourceMarker->SetMarkerStyle(29);
  sourceMarker->SetMarkerSize(2.0);
  sourceMarker->SetMarkerColor(kBlack);
  sourceMarker->Draw("same");

  const std::size_t displayedSpheres =
      std::min<std::size_t>(scene.spheres.size(), 8);
  TLegend *legend = new TLegend(
      0.72, std::max(0.05, 0.95 - 0.04 * (displayedSpheres + 2)), 0.97, 0.95);
  legend->AddEntry(MakeLegendLine(rayColors[0]), "Ray in air", "l");
  for (std::size_t i = 0; i < displayedSpheres; ++i) {
    const SphereDescription &sphere = scene.spheres[i];
    legend->AddEntry(MakeLegendLine(sphere.color), sphere.name.c_str(), "l");
  }
  if (displayedSpheres < scene.spheres.size()) {
    const std::string moreSpheres =
        std::to_string(scene.spheres.size() - displayedSpheres) +
        " more spheres";
    legend->AddEntry(MakeLegendLine(kGray), moreSpheres.c_str(), "l");
  }
  legend->AddEntry(sourceMarker, "Point source", "p");
  legend->Draw();

  const double padding =
      0.08 * std::max({bounds.maximum.X() - bounds.minimum.X(),
                       bounds.maximum.Y() - bounds.minimum.Y(),
                       bounds.maximum.Z() - bounds.minimum.Z(), 1.0});
  view->SetRange(bounds.minimum.X() - padding, bounds.minimum.Y() - padding,
                 bounds.minimum.Z() - padding, bounds.maximum.X() + padding,
                 bounds.maximum.Y() + padding,
                 bounds.maximum.Z() + padding);

  canvas->Modified();
  canvas->Update();
}
