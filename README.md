# medial-axis-3D

`medial-axis-3D` is an experimental C++17 tool for extracting and exploring
an approximate medial-axis sheet complex inside a closed triangle mesh.

The pipeline:

1. reads a TetGen `.node` point set and an optional companion `.face` surface;
2. validates and consistently orients the closed surface;
3. creates deterministic surface samples;
4. builds a 3D Delaunay tetrahedralization and its interior Voronoi dual;
5. validates inward Voronoi poles and medial balls;
6. constructs pole-supported medial sheets;
7. filters them using local feature size, radius continuity, component
   support, and 40/80/160-sample cross-resolution stability; and
8. displays the result interactively with Polyscope.

This is a practical geometric approximation, not an exact symbolic medial-axis
solver. A 3D medial axis is a stratified collection of sheets and curves, so
open sheet boundaries, junctions, and topological holes may be valid results.

## Requirements

- CMake 3.16 or newer
- A C++17 compiler
- Git and internet access during the first CMake configure
- A graphics driver supporting OpenGL 3.3 for the interactive viewer

Polyscope is found as an installed CMake package when available. Otherwise,
CMake downloads and builds it automatically.

## Build

### Windows

Run from PowerShell or a Visual Studio Developer Command Prompt:

```powershell
cmake -S . -B build
cmake --build build --config Release
```

Executables are written to `build\Release\` with the Visual Studio generator.

### Linux or macOS

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

Executables are normally written directly to `build/` with a single-config
generator.

## Quick start

Open the included concave example on Windows:

```powershell
.\build\Release\polyscope_viewer.exe `
  .\examples\dented_cube.node `
  .\build\dented_cube_160.ele `
  --samples 160
```

On Linux or macOS:

```bash
./build/polyscope_viewer \
  ./examples/dented_cube.node \
  ./build/dented_cube_160.ele \
  --samples 160
```

The requested sample resolution is the displayed result. With a valid
companion `.face` file, the tool also evaluates 40, 80, and 160 samples to
estimate cross-resolution stability.

For a faster non-interactive run:

```powershell
.\build\Release\polyscope_viewer.exe `
  .\examples\dented_cube.node `
  --samples 160 `
  --no-gui
```

To save a rendering and exit:

```powershell
.\build\Release\polyscope_viewer.exe `
  .\examples\dented_cube.node `
  --samples 160 `
  --screenshot .\build\dented_cube.png
```

## Input format

The required `.node` file uses the
[TetGen node format](https://www.wias-berlin.de/software/tetgen/fformats.node.html):

```text
<# points> 3 <# attributes> <0 or 1 boundary markers>
<point id> <x> <y> <z> [attributes] [boundary marker]
```

For a non-convex object, provide a triangle surface with the same base name and
a `.face` extension:

```text
<# faces> <0 or 1 boundary markers>
<face id> <node id> <node id> <node id> [boundary marker]
```

For example:

```text
model.node
model.face
```

The surface must be:

- closed and watertight;
- triangular;
- free of non-manifold edges;
- indexed using node IDs from the companion `.node` file.

Face orientation does not need to be consistent on input; the loader repairs
consistent orientation and orients closed components outward.

Without a companion `.face`, the program uses the convex Delaunay boundary.
Surface resampling and cross-resolution analysis require a valid `.face` file.

## Using downloaded 3D models

Most public datasets provide OBJ, OFF, PLY, or STL files. ASCII OFF files can
be converted directly with the included dependency-free Python utility:

```powershell
python .\tools\off_to_tetgen.py .\data\model.off
```

This writes `model.node` and `model.face`. Run the viewer using the generated
node file:

```powershell
.\build\Release\polyscope_viewer.exe `
  .\data\model.node `
  .\build\model_medial.ele `
  --samples 160
```

Alternatively, convert supported surface formats with
[TetGen](https://www.wias-berlin.de/software/tetgen/):

```powershell
tetgen -p .\data\model.off
```

TetGen normally produces:

```text
model.1.node
model.1.face
model.1.ele
```

Use the generated `.node`; the viewer automatically finds the matching `.face`:

```powershell
.\build\Release\polyscope_viewer.exe `
  .\data\model.1.node `
  .\build\model_medial.ele `
  --samples 160
```

Good sources for test meshes include:

- [Princeton Mesh Segmentation Benchmark](https://segeval.cs.princeton.edu/)
  for a downloadable collection of watertight OFF meshes;
- [Stanford 3D Scanning Repository](https://graphics.stanford.edu/data/3Dscanrep/)
  for scanned physical objects;
- [Thingi10K](https://github.com/Thingi10K/Thingi10K) for a wide range of
  3D-printing meshes; and
- [ABC Dataset](https://deep-geometry.github.io/abc-dataset/) for CAD parts.

Repair holes and non-manifold edges before conversion. Large scanned models
should be decimated to a few thousand or tens of thousands of faces because
the current exact surface queries prioritize clarity over large-mesh
performance.

## Command-line options

```text
polyscope_viewer <input.node> [output.ele] [options]
```

| Option | Meaning |
| --- | --- |
| `--samples N` | Generate at least `N` deterministic area-weighted surface samples. Original mesh vertices are always preserved. |
| `--no-gui` | Run the complete pipeline, print statistics, and exit. |
| `--screenshot FILE` | Render one image to `FILE` and exit. |
| `--no-cross-resolution` | Disable the additional 40/80/160 stability runs. |
| `--support-rings N` | Override pole-support propagation rings. `0` uses automatic resolution scaling. |
| `--max-relative-radius-jump X` | Set the maximum allowed relative medial-radius jump. |
| `--fixed-radius-jump` | Use one global radius-jump threshold instead of adaptive LFS thresholds. |
| `--min-component-triangles N` | Set the minimum component triangle count. |
| `--min-component-area-fraction X` | Set the minimum component area fraction from `0` to `1`. |
| `--min-component-confidence X` | Set the minimum mean confidence from `0` to `1`. |
| `--min-component-poles N` | Require at least `N` directly supporting validated poles. |

To make the radius-discontinuity weight maximally permissive:

```text
--fixed-radius-jump --max-relative-radius-jump 1
```

When `output.ele` is supplied, it contains the generated Delaunay tetrahedra,
not the medial-sheet triangles. When resampling is enabled, a matching
`output.node` is written beside it.

## Polyscope controls and layers

The main result is `pole-supported medial sheets`. Useful quantities and
diagnostic structures include:

- `sheet confidence`
- `contact angle`
- `medial radius`
- `relative radius jump`
- `surface resolution`
- `estimated LFS`
- `sampling density h/LFS`
- `cross-resolution stability`
- `pole support weight`
- `medial sheet boundaries`
- `validated sheet terminations`
- `unresolved artificial sheet boundaries`
- `medial sheet seams`
- `medial sheet junctions`
- `restored rejected Voronoi faces`
- `radius-discontinuity weights`
- `removed unstable sheet components`
- `removed resolution-unstable sheets`
- `validated medial balls`
- `contact points` and `contact spokes`
- `rejected sheet candidates`

Most diagnostic layers start disabled. Enable them in Polyscope's structure
list when investigating a result.

The adaptive radius threshold uses local sampling density `d = h / LFS`:

```text
threshold = clamp(0.35 + 0.60 * d / (1 + d), 0.35, 0.85)
```

Cross-resolution stability matches medial vertices by LFS-normalized position
and radius, then matches sheet components using overlap, orientation, pole
support, and topology. Components supported at multiple resolutions retain
their complete triangulation. The 40, 80, and 160 comparison runs use exact
downsample counts; the displayed run may contain more points because it
preserves every original mesh vertex.

Pole support, confidence, and radius continuity are weights rather than
per-triangle deletion tests. Enclosed gaps are restored as connected patches
when their radius, orientation, and topology agree. Boundary loops are
classified as real terminations, seams/junctions, or unresolved artificial
boundaries, and cross-resolution instability removes complete sheet
components instead of isolated triangles.

## Tests

Build and run the regression suite:

```powershell
cmake --build build --config Release --target test_runner
.\build\Release\test_runner.exe
```

Or use CTest:

```powershell
ctest --test-dir build -C Release --output-on-failure
```

The regression suite covers Delaunay construction, explicit surface handling,
sampling, pole validation, medial-sheet construction, adaptive radius
filtering, cross-resolution stability, and topology-preserving gap handling.

## Repository layout

```text
examples/       Small ready-to-run .node + .face examples
src/core/       Geometry, Delaunay, Voronoi, medial-sheet, and filtering code
src/io/         TetGen input/output support
src/viewer/     Polyscope command-line viewer
tests/          Regression test runner
```

The original project description is available in
[`Medial Axis for 3D Project.pdf`](./Medial%20Axis%20for%203D%20Project.pdf).
