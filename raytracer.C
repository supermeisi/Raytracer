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

#include <cmath>
#include <vector>

namespace {

constexpr double kSphereRadius = 3.0;
constexpr double kGlassIndex = 1.5;
constexpr double kPathLength = 20.0;
constexpr double kSurfaceEpsilon = 1e-5;

struct RaySegment {
  TVector3 start;
  TVector3 end;
  bool inGlass;
};

using RayPath = std::vector<RaySegment>;

bool IntersectSphere(const TVector3 &origin, const TVector3 &direction,
                     double radius, double &distance, TVector3 &normal) {
  const TVector3 offset = origin;
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

  normal = (origin + direction * distance) * (1.0 / radius);
  return true;
}

RayPath TraceRay(const TVector3 &source, TVector3 direction) {
  direction = direction.Unit();
  TVector3 origin = source;
  RayPath path;
  double travelled = 0.0;
  bool inGlass = false;

  while (travelled < kPathLength) {
    double hitDistance = 0.0;
    TVector3 outwardNormal;
    const double remaining = kPathLength - travelled;
    if (!IntersectSphere(origin, direction, kSphereRadius, hitDistance,
                         outwardNormal) ||
        hitDistance >= remaining) {
      path.push_back({origin, origin + direction * remaining, inGlass});
      break;
    }

    const TVector3 hit = origin + direction * hitDistance;
    path.push_back({origin, hit, inGlass});
    travelled += hitDistance;

    const bool entering = !inGlass;
    const TVector3 normalAgainstRay =
        entering ? outwardNormal : -outwardNormal;
    const double incidentIndex = entering ? 1.0 : kGlassIndex;
    const double transmittedIndex = entering ? kGlassIndex : 1.0;
    const double eta = incidentIndex / transmittedIndex;
    const double cosIncident = -direction.Dot(normalAgainstRay);
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
      inGlass = !inGlass;
    }

    origin = hit + direction * kSurfaceEpsilon;
    travelled += kSurfaceEpsilon;
  }

  return path;
}

void DrawSphereWireframe(double radius) {
  constexpr int kSamples = 80;
  constexpr int kLatitudeCount = 12;
  constexpr int kLongitudeCount = 16;

  for (int latitude = 1; latitude < kLatitudeCount; ++latitude) {
    const double polar =
        TMath::Pi() * latitude / static_cast<double>(kLatitudeCount);
    TPolyLine3D *line = new TPolyLine3D(kSamples + 1);
    for (int i = 0; i <= kSamples; ++i) {
      const double azimuth = 2.0 * TMath::Pi() * i / kSamples;
      line->SetPoint(i, radius * std::sin(polar) * std::cos(azimuth),
                     radius * std::sin(polar) * std::sin(azimuth),
                     radius * std::cos(polar));
    }
    line->SetLineColor(kGray + 1);
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
      line->SetPoint(i, radius * std::sin(polar) * std::cos(azimuth),
                     radius * std::sin(polar) * std::sin(azimuth),
                     radius * std::cos(polar));
    }
    line->SetLineColor(kGray + 1);
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->Draw("same");
  }
}

void DrawSegment(const RaySegment &segment, Color_t outsideColor) {
  TPolyLine3D *line = new TPolyLine3D(2);
  line->SetPoint(0, segment.start.X(), segment.start.Y(), segment.start.Z());
  line->SetPoint(1, segment.end.X(), segment.end.Y(), segment.end.Z());
  line->SetLineColor(segment.inGlass ? kOrange + 7 : outsideColor);
  line->SetLineWidth(2);
  line->Draw("same");
}

} // namespace

void raytracer() {
  gROOT->SetStyle("Plain");
  gStyle->SetCanvasColor(kWhite);

  TCanvas *canvas =
      new TCanvas("raytracer_canvas", "ROOT 3D Ray Tracer", 1000, 800);
  TView *view = TView::CreateView(1);
  view->SetRange(-5.0, -5.0, -9.0, 5.0, 5.0, 10.0);

  DrawSphereWireframe(kSphereRadius);

  const TVector3 source(0.0, 0.0, -8.0);
  const int rayColors[] = {kAzure + 2, kBlue + 1, kCyan + 2, kMagenta + 1,
                           kGreen + 2, kViolet + 1, kRed + 1};
  constexpr int kGridRadius = 3;
  constexpr double kTargetZ = -3.5;
  constexpr double kTargetSpacing = 0.6;

  for (int y = -kGridRadius; y <= kGridRadius; ++y) {
    for (int x = -kGridRadius; x <= kGridRadius; ++x) {
      const TVector3 target(x * kTargetSpacing, y * kTargetSpacing, kTargetZ);
      const RayPath path = TraceRay(source, target - source);
      const Color_t color = rayColors[(x + y + 2 * kGridRadius) % 7];
      for (const RaySegment &segment : path) {
        DrawSegment(segment, color);
      }
    }
  }

  TPolyMarker3D *sourceMarker = new TPolyMarker3D(1);
  sourceMarker->SetPoint(0, source.X(), source.Y(), source.Z());
  sourceMarker->SetMarkerStyle(29);
  sourceMarker->SetMarkerSize(2.0);
  sourceMarker->SetMarkerColor(kBlack);
  sourceMarker->Draw("same");

  TPolyLine3D *outsideSample = new TPolyLine3D(2);
  outsideSample->SetPoint(0, -4.5, -4.5, -8.0);
  outsideSample->SetPoint(1, -4.0, -4.0, -7.5);
  outsideSample->SetLineColor(kAzure + 2);
  outsideSample->SetLineWidth(2);

  TPolyLine3D *insideSample = new TPolyLine3D(2);
  insideSample->SetPoint(0, -4.5, -4.5, -8.0);
  insideSample->SetPoint(1, -4.0, -4.0, -7.5);
  insideSample->SetLineColor(kOrange + 7);
  insideSample->SetLineWidth(2);

  TLegend *legend = new TLegend(0.72, 0.82, 0.97, 0.95);
  legend->AddEntry(outsideSample, "Ray in air", "l");
  legend->AddEntry(insideSample, "Ray in glass", "l");
  legend->AddEntry(sourceMarker, "Point source", "p");
  legend->Draw();

  canvas->Modified();
  canvas->Update();
}
