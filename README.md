# ROOT 3D Ray Tracer

This ROOT macro traces a grid of light rays from a point source through a
transparent glass sphere. It uses Snell's law at the sphere boundary, reflects
rays that undergo total internal reflection, and displays the ray paths and a
wireframe sphere in an interactive 3D ROOT canvas.

## Build and run

With ROOT installed, configure and build the project:

```sh
cmake -S . -B build -DCMAKE_PREFIX_PATH=/path/to/ROOT
cmake --build build
./build/raytracer
```

If ROOT's `ROOTConfig.cmake` is already in CMake's package search path, omit
`-DCMAKE_PREFIX_PATH`. Alternatively, set `ROOT_DIR` to the directory
containing `ROOTConfig.cmake`.

The 3D view can be rotated with the mouse. Edit the constants near the top of
`raytracer.C` to adjust the sphere radius, glass refractive index, or traced
path length; the source position and ray-grid spacing are set in `raytracer()`.
The original macro can also be run directly from ROOT with `root -l raytracer.C`.

The tracer follows one path per emitted ray: it handles refraction and total
internal reflection, but does not split rays into separate Fresnel-reflected
and transmitted branches.
