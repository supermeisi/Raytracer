# ROOT 3D Ray Tracer

This ROOT application traces a configurable beam of light rays from a point
source through a transparent glass scene. It uses Snell's law at sphere
boundaries, reflects rays that undergo total internal reflection, and displays
the ray paths and wireframe objects in an interactive 3D ROOT canvas.

## Build and run

With ROOT installed, configure and build the project:

```sh
cmake -S . -B build -DCMAKE_PREFIX_PATH=/path/to/ROOT
cmake --build build
./build/raytracer
```

If ROOT's `ROOTConfig.cmake` is already in CMake's package search path, omit
`-DCMAKE_PREFIX_PATH`. Alternatively, set `ROOT_DIR` to the directory
containing `ROOTConfig.cmake`. The JSON parser dependency is nlohmann-json
(for example, install `nlohmann-json3-dev` on Debian or Ubuntu).

The default [scene.json](./scene.json) is loaded when no argument is given.
Pass another scene description to load it instead:

```sh
./build/raytracer path/to/scene.json
```

The JSON scene format configures the source position and beam, path length,
Fresnel branching controls, an optional finite plane detector, and an array of
non-overlapping spheres. The source's `ray_count` sets the exact number of
rays, `direction` aims the center ray, and `opening_angle_degrees` sets the full
cone angle (from 0 to 180 degrees). Rays are distributed over the cone using a
deterministic spiral pattern, with the first ray following the configured
direction. The `optics.worker_count` setting controls the number of tracing
threads: `0` (the default) uses the available hardware concurrency, and a
positive value selects an explicit worker count from 1 to 256. The terminal
shows the percentage and count of source rays completed; all ROOT drawing
remains on the main thread. Each sphere has a center, radius, refractive index,
and ROOT color index.

The detector has a center, surface normal, in-plane `u_axis`, physical width
and height, histogram bin counts, and outline color. The `u_axis` is projected
onto the detector plane and the perpendicular `v` axis is derived from the
normal. Rays hitting the detector are absorbed and their branch intensities
are accumulated in the corresponding 2D histogram bin. Its outline appears in
the 3D scene, with the ROOT `COLZ` histogram in a second canvas.

```json
{
  "source": {
    "position": [0, 0, -8],
    "direction": [0, 0, 1],
    "ray_count": 49,
    "opening_angle_degrees": 60
  },
  "path_length": 20,
  "optics": {
    "minimum_branch_intensity": 0.01,
    "max_boundary_interactions": 8,
    "worker_count": 0
  },
  "detector": {
    "name": "Exit plane detector",
    "center": [0, 0, 8],
    "normal": [0, 0, -1],
    "u_axis": [1, 0, 0],
    "width": 20,
    "height": 20,
    "bins_x": 100,
    "bins_y": 100,
    "color": 2
  },
  "objects": [
    {
      "type": "sphere",
      "name": "Glass sphere",
      "center": [0, 0, 0],
      "radius": 3,
      "refractive_index": 1.5,
      "color": 800
    }
  ]
}
```

Set `"detector": null` or omit `detector` to disable detection. Detector bin
counts must be between 1 and 1000 per axis, with at most 1,000,000 bins total.
Only spheres are currently supported as optical objects. Spheres may not
overlap or touch, and the source must be outside every sphere. The direction
vector is normalized automatically. Ray count is limited to 1–10000. The 3D
view can be rotated with the mouse.

Run the JSON loader tests with `ctest --test-dir build --output-on-failure`.

At each boundary, the unpolarized Fresnel equations split ray power into
reflected and transmitted branches; total internal reflection sends all power
into the reflected branch. Branch line width gives a visual indication of
relative intensity. To keep the render bounded, branches below
`optics.minimum_branch_intensity` are omitted, and each path ends after
`optics.max_boundary_interactions` (1–12) boundary encounters. Detector bins
sum the intensity of every Fresnel branch that reaches the plane.
