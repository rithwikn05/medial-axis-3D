# medial-axis-3D

This repository now contains a first-pass scaffold for a 3D medial-axis pipeline:

- basic 3D vector utilities in src/core/vec3.h
- an initial tetrahedralization scaffold in src/core/delaunay3.h
- placeholder structures for surface samples, Voronoi duals, medial candidates, and filters
- a small inspection executable in src/app/inspect.cpp
- a simple test runner in tests/test_runner.cpp
- a simple test runner in tests/test_runner.cpp

## What we've done so far

- Implemented `Vec3` utilities and helpers in `src/core/vec3.h`.
- Added an initial Delaunay scaffold in `src/core/delaunay3.h` and a simple incremental `insert` API.
- Implemented a circumsphere predicate (`compute_circumsphere`) and `circumsphere_contains(...)` test.
- Implemented cavity detection (`cavity_indices(...)`) for identifying tetrahedra whose circumspheres contain a point.
- Implemented tetrahedron replacement: remove cavity tetrahedra and create a local patch that connects boundary faces to the new point.
- Added basic tetrahedron adjacency (`neighbors`) and a small history DAG (`parents`) per tetrahedron; added `rebuild_adjacency()`.
- Switched to a pivoted linear solver for the circumsphere to improve numerical stability in the scaffold.
- Added `tests/test_runner.cpp` and `src/app/inspect.cpp` and verified the scaffold compiles and runs on the local build.

Build and run locally with the Visual Studio Developer Command Prompt:

```bat
cl /nologo /EHsc /I src tests\test_runner.cpp /Fe:test_runner.exe
.test_runner.exe

cl /nologo /EHsc /I src src\app\inspect.cpp /Fe:inspect.exe
inspect.exe
```

The current Delaunay implementation is intentionally simple and testable. The next step is to replace it with a true local-update Delaunay routine that removes tetrahedra whose circumspheres contain the newly inserted point, matching the conceptual step from the PDF.

## Remaining work / Future improvements

- Implement tetrahedron adjacency and a history DAG (parents) for efficient local updates and rollback. (Basic adjacency and parent links have been added to `src/core/delaunay3.h`.)
- Swap the numeric circumsphere solver for robust predicates (`orient3d`, `insphere`) — consider Shewchuk's `predicates.c`.
- Implement incremental flips and use adjacency to perform localized Delaunay maintenance instead of global rebuilds.
- Build the Voronoi dual from the Delaunay tetrahedralization and expose Voronoi cells for medial analysis.
- Implement pole selection and medial sheet extraction per the project PDF.
- Add filtering and simplification to produce a clean medial complex suitable for visualization.
- Improve test coverage with larger geometry cases and add performance benchmarks.
- Performance improvements: spatial acceleration structures, incremental caches, and memory layout optimizations.

If you'd like, I can expand any bullet into a concrete implementation plan and start the next item.
