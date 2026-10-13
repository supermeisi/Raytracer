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

  void IncludeSphere(const TVector3 &center, double radius) {
    Include(center - TVector3(radius, radius, radius));
    Include(center + TVector3(radius, radius, radius));
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

void DrawSegment(const RaySegment &segment, Color_t outsideColor,
                 const SceneConfig &scene, SceneBounds &bounds) {
  TPolyLine3D *line = new TPolyLine3D(2);
  line->SetPoint(0, segment.start[0], segment.start[1], segment.start[2]);
  line->SetPoint(1, segment.end[0], segment.end[1], segment.end[2]);
  const Color_t color = segment.sphereIndex >= 0
                            ? scene.spheres[segment.sphereIndex].color
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

  for (const SphereDescription &sphere : scene.spheres) {
    const TVector3 center = ToVector3(sphere.center);
    DrawSphereWireframe(center, sphere.radius, sphere.color);
    bounds.IncludeSphere(center, sphere.radius);
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
