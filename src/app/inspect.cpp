#include "../core/delaunay3.h"

#include <iostream>
#include <vector>

int main() {
    using namespace medial_axis_3d;

    Delaunay3 delaunay;
    const std::vector<Vec3> samples = {
        Vec3(0.0, 0.0, 0.0),
        Vec3(1.0, 0.0, 0.0),
        Vec3(0.0, 1.0, 0.0),
        Vec3(0.0, 0.0, 1.0),
        Vec3(0.2, 0.2, 0.2)
    };

    for (const auto& sample : samples) {
        delaunay.insert(sample);
    }

    std::cout << "Inserted points: " << delaunay.point_count() << '\n';
    std::cout << "Current tetrahedra: " << delaunay.tetrahedron_count() << '\n';

    for (std::size_t i = 0; i < delaunay.tetrahedra().size(); ++i) {
        const auto& tetrahedron = delaunay.tetrahedra()[i];
        std::cout << "tetrahedron " << i << ": "
                  << tetrahedron.vertices[0] << ", "
                  << tetrahedron.vertices[1] << ", "
                  << tetrahedron.vertices[2] << ", "
                  << tetrahedron.vertices[3] << '\n';
    }

    return 0;
}
