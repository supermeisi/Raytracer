# ROOT 3D Ray Tracer

This ROOT macro traces a grid of light rays from a point source through a
transparent glass scene. It uses Snell's law at sphere boundaries, reflects
rays that undergo total internal reflection, and displays the ray paths and
wireframe objects in an interactive 3D ROOT canvas.

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

The JSON scene format has a source, ray grid, path length, and an array of
non-overlapping spheres. Each sphere has a center, radius, refractive index,
and ROOT color index:

```json
{
  "source": {
    "position": [0, 0, -8],
    "grid": { "radius": 3, "target_z": -3.5, "spacing": 0.6 }
  },
  "path_length": 20,
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

Only spheres are currently supported. Spheres may not overlap or touch, and
the source must be outside every sphere. `grid.radius` controls the number of
rays in each direction around the center ray; the scene values are validated
when loaded. The 3D view can be rotated with the mouse.

Run the JSON loader tests with `ctest --test-dir build --output-on-failure`.

The tracer follows one path per emitted ray: it handles refraction and total
internal reflection, but does not split rays into separate Fresnel-reflected
and transmitted branches.
