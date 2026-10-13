#include "SceneConfig.h"

#include "Detector.h"
#include "RayTracingEngine.h"

#include <TCanvas.h>
#include <TColor.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TMath.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TVector3.h>
#include <TView.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <future>
#include <iostream>
#include <limits>
#include <memory>
#include <thread>
#include <vector>

namespace {

using RaySegment = TracedRaySegment;

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

  void IncludeObject(const ObjectDescription &object) {
    const TVector3 center(object.center[0], object.center[1], object.center[2]);
    if (object.shape == ObjectShape::Step) {
      for (const MeshTriangle &triangle : object.triangles) {
        for (const auto *vertex : {&triangle.a, &triangle.b, &triangle.c}) {
          Include(center + TVector3((*vertex)[0], (*vertex)[1], (*vertex)[2]));
        }
      }
      return;
    }
    if (object.shape == ObjectShape::Sphere) {
      const double radius = object.radius;
      Include(center - TVector3(radius, radius, radius));
      Include(center + TVector3(radius, radius, radius));
      return;
    }
    const double radius =
        object.shape == ObjectShape::Box
            ? 0.0
            : std::max(object.radius, object.topRadius);
    const TVector3 halfSize =
        object.shape == ObjectShape::Box
            ? TVector3(object.size[0] * 0.5, object.size[1] * 0.5,
                       object.size[2] * 0.5)
            : TVector3(radius, radius, object.size[2] * 0.5);
    Include(center - halfSize);
    Include(center + halfSize);
  }
};

TVector3 ToVector3(const std::array<double, 3> &vector) {
  return TVector3(vector[0], vector[1], vector[2]);
}

void DrawDetectorPlane(const PlaneDetectorDescription &detector,
                       SceneBounds &bounds) {
  const TVector3 center = ToVector3(detector.center);
  const TVector3 normal = ToVector3(detector.normal).Unit();
  TVector3 uAxis = ToVector3(detector.uAxis);
  uAxis = (uAxis - normal * uAxis.Dot(normal)).Unit();
  const TVector3 vAxis = normal.Cross(uAxis).Unit();
  const TVector3 halfU = uAxis * (detector.width * 0.5);
  const TVector3 halfV = vAxis * (detector.height * 0.5);
  const TVector3 corners[] = {center - halfU - halfV, center + halfU - halfV,
                              center + halfU + halfV, center - halfU + halfV,
                              center - halfU - halfV};

  TPolyLine3D *outline = new TPolyLine3D(5);
  for (int i = 0; i < 5; ++i) {
    outline->SetPoint(i, corners[i].X(), corners[i].Y(), corners[i].Z());
    bounds.Include(corners[i]);
  }
  outline->SetLineColor(detector.color);
  outline->SetLineWidth(2);
  outline->Draw("same");
}

void DrawDetectorHistogram(const PlaneDetectorHistogram &detectorHistogram) {
  const PlaneDetectorDescription &detector = detectorHistogram.Description();
  TCanvas *canvas =
      new TCanvas("detector_histogram_canvas", detector.name.c_str(), 900, 700);
  TH2D *histogram = new TH2D(
      "detected_light", (detector.name + ";Detector u;Detector v").c_str(),
      detector.binsX, -detector.width * 0.5, detector.width * 0.5,
      detector.binsY, -detector.height * 0.5, detector.height * 0.5);
  histogram->SetDirectory(nullptr);
  const std::vector<double> &bins = detectorHistogram.Bins();
  for (int y = 0; y < detector.binsY; ++y) {
    for (int x = 0; x < detector.binsX; ++x) {
      const std::size_t index = static_cast<std::size_t>(y) *
                                    static_cast<std::size_t>(detector.binsX) +
                                static_cast<std::size_t>(x);
      histogram->SetBinContent(x + 1, y + 1, bins[index]);
    }
  }
  histogram->SetStats(false);
  histogram->Draw("COLZ");
  canvas->Modified();
  canvas->Update();
}

void DrawCircle(const TVector3 &center, double radius, Color_t color,
                bool closed) {
  if (radius <= 0.0) {
    return;
  }
  constexpr int kSamples = 48;
  TPolyLine3D *line = new TPolyLine3D(kSamples + (closed ? 1 : 0));
  const int pointCount = kSamples + (closed ? 1 : 0);
  for (int i = 0; i < pointCount; ++i) {
    const double angle = 2.0 * TMath::Pi() * i / kSamples;
    line->SetPoint(i, center.X() + radius * std::cos(angle),
                   center.Y() + radius * std::sin(angle), center.Z());
  }
  line->SetLineColor(color);
  line->SetLineStyle(2);
  line->Draw("same");
}

void DrawSphereWireframe(const TVector3 &center, double radius, Color_t color) {
  constexpr int kSamples = 80;
  constexpr int kLatitudeCount = 12;
  constexpr int kLongitudeCount = 16;

  for (int latitude = 1; latitude < kLatitudeCount; ++latitude) {
    const double polar =
        TMath::Pi() * latitude / static_cast<double>(kLatitudeCount);
    TPolyLine3D *line = new TPolyLine3D(kSamples + 1);
    for (int i = 0; i <= kSamples; ++i) {
      const double azimuth = 2.0 * TMath::Pi() * i / kSamples;
      line->SetPoint(i,
                     center.X() + radius * std::sin(polar) * std::cos(azimuth),
                     center.Y() + radius * std::sin(polar) * std::sin(azimuth),
                     center.Z() + radius * std::cos(polar));
    }
    line->SetLineColor(color);
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->Draw("same");
  }

  for (int longitude = 0; longitude < kLongitudeCount; ++longitude) {
    const double azimuth = 2.0 * TMath::Pi() * longitude / kLongitudeCount;
    TPolyLine3D *line = new TPolyLine3D(kSamples + 1);
    for (int i = 0; i <= kSamples; ++i) {
      const double polar = TMath::Pi() * i / kSamples;
      line->SetPoint(i,
                     center.X() + radius * std::sin(polar) * std::cos(azimuth),
                     center.Y() + radius * std::sin(polar) * std::sin(azimuth),
                     center.Z() + radius * std::cos(polar));
    }
    line->SetLineColor(color);
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->Draw("same");
  }
}

void DrawOpticalObject(const ObjectDescription &object) {
  const TVector3 center = ToVector3(object.center);
  if (object.shape == ObjectShape::Sphere) {
    DrawSphereWireframe(center, object.radius, object.color);
  } else if (object.shape == ObjectShape::Step) {
    for (const MeshTriangle &triangle : object.triangles) {
      TPolyLine3D *line = new TPolyLine3D(4);
      const std::array<std::array<double, 3>, 4> vertices{
          triangle.a, triangle.b, triangle.c, triangle.a};
      for (int vertex = 0; vertex < 4; ++vertex) {
        line->SetPoint(vertex,
                       center.X() + vertices[vertex][0],
                       center.Y() + vertices[vertex][1],
                       center.Z() + vertices[vertex][2]);
      }
      line->SetLineColor(object.color);
      line->SetLineStyle(2);
      line->Draw("same");
    }
  } else if (object.shape == ObjectShape::Box) {
    const TVector3 half(object.size[0] * 0.5, object.size[1] * 0.5,
                        object.size[2] * 0.5);
    TVector3 corners[8];
    for (int corner = 0; corner < 8; ++corner) {
      corners[corner] =
          center + TVector3((corner & 1) ? half.X() : -half.X(),
                            (corner & 2) ? half.Y() : -half.Y(),
                            (corner & 4) ? half.Z() : -half.Z());
    }
    constexpr int edges[][2] = {{0, 1}, {2, 3}, {4, 5}, {6, 7},
                                {0, 2}, {1, 3}, {4, 6}, {5, 7},
                                {0, 4}, {1, 5}, {2, 6}, {3, 7}};
    for (const auto &edge : edges) {
      TPolyLine3D *line = new TPolyLine3D(2);
      line->SetPoint(0, corners[edge[0]].X(), corners[edge[0]].Y(),
                     corners[edge[0]].Z());
      line->SetPoint(1, corners[edge[1]].X(), corners[edge[1]].Y(),
                     corners[edge[1]].Z());
      line->SetLineColor(object.color);
      line->SetLineStyle(2);
      line->Draw("same");
    }
  } else {
    const double halfHeight = object.size[2] * 0.5;
    const TVector3 bottom(center.X(), center.Y(), center.Z() - halfHeight);
    const TVector3 top(center.X(), center.Y(), center.Z() + halfHeight);
    const double topRadius =
        object.shape == ObjectShape::Tube ? object.radius : object.topRadius;
    DrawCircle(bottom, object.radius, object.color, true);
    DrawCircle(top, topRadius, object.color, true);
    constexpr int kMeridians = 12;
    for (int i = 0; i < kMeridians; ++i) {
      const double angle = 2.0 * TMath::Pi() * i / kMeridians;
      TPolyLine3D *line = new TPolyLine3D(2);
      line->SetPoint(0, bottom.X() + object.radius * std::cos(angle),
                     bottom.Y() + object.radius * std::sin(angle), bottom.Z());
      line->SetPoint(1, top.X() + topRadius * std::cos(angle),
                     top.Y() + topRadius * std::sin(angle), top.Z());
      line->SetLineColor(object.color);
      line->SetLineStyle(2);
      line->Draw("same");
    }
  }
}

void DrawSegment(const RaySegment &segment, Color_t outsideColor,
                 const SceneConfig &scene, SceneBounds &bounds) {
  TPolyLine3D *line = new TPolyLine3D(2);
  line->SetPoint(0, segment.start[0], segment.start[1], segment.start[2]);
  line->SetPoint(1, segment.end[0], segment.end[1], segment.end[2]);
  const Color_t color = segment.objectIndex >= 0
                            ? scene.objects[segment.objectIndex].color
                            : outsideColor;
  line->SetLineColor(color);
  line->SetLineWidth(std::max(
      1,
      static_cast<int>(std::lround(1.0 + 4.0 * std::sqrt(segment.intensity)))));
  line->Draw("same");
  bounds.Include(ToVector3(segment.start));
  bounds.Include(ToVector3(segment.end));
}

void PrintProgress(std::size_t completed, std::size_t total) {
  const int percentage =
      total == 0 ? 100 : static_cast<int>(100 * completed / total);
  std::cout << "\rTracing rays: " << percentage << "% (" << completed << "/"
            << total << ")" << std::flush;
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
  std::unique_ptr<PlaneDetectorHistogram> detectorHistogram;

  for (const ObjectDescription &object : scene.objects) {
    DrawOpticalObject(object);
    bounds.IncludeObject(object);
  }
  if (scene.detector) {
    detectorHistogram =
        std::make_unique<PlaneDetectorHistogram>(*scene.detector);
    DrawDetectorPlane(*scene.detector, bounds);
  }

  const int rayColors[] = {kAzure + 2, kBlue + 1,   kCyan + 2, kMagenta + 1,
                           kGreen + 2, kViolet + 1, kRed + 1};

  const std::vector<std::array<double, 3>> rayDirections =
      GenerateRayDirections(scene);
  std::atomic<std::size_t> completedRayCount{0};
  auto tracing = std::async(std::launch::async, [&] {
    return TraceRaysParallel(scene, rayDirections, scene.workerCount,
                             completedRayCount);
  });
  PrintProgress(0, rayDirections.size());
  int lastDisplayedPercentage = 0;
  auto lastProgressUpdate = std::chrono::steady_clock::now();
  while (tracing.wait_for(std::chrono::milliseconds(0)) !=
         std::future_status::ready) {
    const auto now = std::chrono::steady_clock::now();
    if (now - lastProgressUpdate >= std::chrono::milliseconds(100)) {
      const std::size_t completed =
          completedRayCount.load(std::memory_order_relaxed);
      const int percentage = static_cast<int>(
          100 * completed / rayDirections.size());
      if (percentage != lastDisplayedPercentage) {
        PrintProgress(completed, rayDirections.size());
        lastDisplayedPercentage = percentage;
      }
      lastProgressUpdate = now;
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
  }
  const std::vector<TracedRay> tracedRays = tracing.get();
  PrintProgress(rayDirections.size(), rayDirections.size());
  std::cout << '\n';

  for (std::size_t rayIndex = 0; rayIndex < tracedRays.size(); ++rayIndex) {
    const TracedRay &ray = tracedRays[rayIndex];
    const Color_t color = rayColors[rayIndex % 7];
    for (const RaySegment &segment : ray.segments) {
      DrawSegment(segment, color, scene, bounds);
    }
    if (detectorHistogram && ray.detectorHit) {
      detectorHistogram->Record(*ray.detectorHit, ray.detectorIntensity);
    }
  }

  TPolyMarker3D *sourceMarker = new TPolyMarker3D(1);
  sourceMarker->SetPoint(0, source.X(), source.Y(), source.Z());
  sourceMarker->SetMarkerStyle(29);
  sourceMarker->SetMarkerSize(2.0);
  sourceMarker->SetMarkerColor(kBlack);
  sourceMarker->Draw("same");

  const std::size_t displayedObjects =
      std::min<std::size_t>(scene.objects.size(), 8);
  TLegend *legend = new TLegend(
      0.72, std::max(0.05, 0.95 - 0.04 * (displayedObjects + 2)), 0.97, 0.95);
  legend->AddEntry(MakeLegendLine(rayColors[0]), "Ray in air", "l");
  for (std::size_t i = 0; i < displayedObjects; ++i) {
    const ObjectDescription &object = scene.objects[i];
    legend->AddEntry(MakeLegendLine(object.color), object.name.c_str(), "l");
  }
  if (displayedObjects < scene.objects.size()) {
    const std::string moreObjects =
        std::to_string(scene.objects.size() - displayedObjects) +
        " more objects";
    legend->AddEntry(MakeLegendLine(kGray), moreObjects.c_str(), "l");
  }
  legend->AddEntry(sourceMarker, "Point source", "p");
  legend->AddEntry(MakeLegendLine(rayColors[0]), "Line width = intensity", "l");
  legend->Draw();

  const double padding =
      0.08 * std::max({bounds.maximum.X() - bounds.minimum.X(),
                       bounds.maximum.Y() - bounds.minimum.Y(),
                       bounds.maximum.Z() - bounds.minimum.Z(), 1.0});
  view->SetRange(bounds.minimum.X() - padding, bounds.minimum.Y() - padding,
                 bounds.minimum.Z() - padding, bounds.maximum.X() + padding,
                 bounds.maximum.Y() + padding, bounds.maximum.Z() + padding);

  canvas->Modified();
  canvas->Update();
  if (detectorHistogram) {
    DrawDetectorHistogram(*detectorHistogram);
  }
}
