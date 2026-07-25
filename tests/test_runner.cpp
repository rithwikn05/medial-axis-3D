#include "../src/core/delaunay3.h"

#include <array>
#include <cassert>
#include <iostream>
#include <vector>

int main() {
    using namespace medial_axis_3d;

    Delaunay3 delaunay;
    assert(delaunay.point_count() == 0);
    assert(delaunay.tetrahedron_count() == 0);

    const std::vector<Vec3> samples = {
        Vec3(0.0, 0.0, 0.0),
        Vec3(1.0, 0.0, 0.0),
        Vec3(0.0, 1.0, 0.0),
        Vec3(0.0, 0.0, 1.0),
        Vec3(0.2, 0.2, 0.2)
    };

    for (const auto& sample : samples) {
        assert(delaunay.insert(sample));
    }

    Delaunay3 initial_tetra;
    for (const auto& sample : std::vector<Vec3>{
        Vec3(0.0, 0.0, 0.0),
        Vec3(1.0, 0.0, 0.0),
        Vec3(0.0, 1.0, 0.0),
        Vec3(0.0, 0.0, 1.0)
    }) {
        initial_tetra.insert(sample);
    }

    const std::array<int, 4> expected_initial_vertices{0, 1, 2, 3};
    assert(initial_tetra.tetrahedra()[0].vertices == expected_initial_vertices);

    assert(delaunay.point_count() == 5);
    assert(delaunay.tetrahedron_count() == 4);
    assert(delaunay.tetrahedra()[0].vertices[0] == 0);
    assert(delaunay.tetrahedra()[3].vertices[3] == 4);

    Delaunay3 cavity_test;
    for (const auto& sample : std::vector<Vec3>{
        Vec3(0.0, 0.0, 0.0),
        Vec3(1.0, 0.0, 0.0),
        Vec3(0.0, 1.0, 0.0),
        Vec3(0.0, 0.0, 1.0)
    }) {
        cavity_test.insert(sample);
    }

    const auto cavity = cavity_test.cavity_indices(Vec3(0.2, 0.2, 0.2));
    assert(!cavity.empty());

    // Check adjacency: each tetrahedron should have neighbors recorded (or -1 for boundary)
    bool any_neighbor_present = false;
    for (std::size_t t = 0; t < delaunay.tetrahedron_count(); ++t) {
        const auto& tet = delaunay.tetrahedra()[t];
        for (int n = 0; n < 4; ++n) {
            const int neigh = tet.neighbors[n];
            if (neigh != -1) {
                any_neighbor_present = true;
                assert(neigh >= 0 && neigh < static_cast<int>(delaunay.tetrahedron_count()));
            }
        }
    }
    assert(any_neighbor_present);

    std::cout << "Delaunay scaffold tests passed.\n";
    return 0;
}
