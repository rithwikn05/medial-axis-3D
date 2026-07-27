#pragma once

#include <cstddef>
#include <vector>

#include "delaunay3.h"
#include "vec3.h"

namespace medial_axis_3d {

struct VoronoiCell {
    int owner_point = -1;
    Vec3 center{};
    std::vector<Vec3> vertices;
    std::vector<int> incident_tetrahedra;
};

struct VoronoiDual {
    std::vector<VoronoiCell> cells;
};

inline VoronoiDual build_voronoi_dual(const Delaunay3& delaunay) {
    VoronoiDual dual;
    dual.cells.reserve(delaunay.point_count());

    for (std::size_t point_index = 0; point_index < delaunay.point_count(); ++point_index) {
        VoronoiCell cell;
        cell.owner_point = static_cast<int>(point_index);
        cell.center = delaunay.points()[point_index];

        for (std::size_t tetrahedron_index = 0; tetrahedron_index < delaunay.tetrahedron_count(); ++tetrahedron_index) {
            const auto& tetrahedron = delaunay.tetrahedra()[tetrahedron_index];
            const bool incident_to_point = tetrahedron.vertices[0] == static_cast<int>(point_index) ||
                tetrahedron.vertices[1] == static_cast<int>(point_index) ||
                tetrahedron.vertices[2] == static_cast<int>(point_index) ||
                tetrahedron.vertices[3] == static_cast<int>(point_index);

            if (!incident_to_point) {
                continue;
            }

            cell.incident_tetrahedra.push_back(static_cast<int>(tetrahedron_index));

            Vec3 circumcenter{};
            if (delaunay.circumcenter(tetrahedron_index, circumcenter)) {
                cell.vertices.push_back(circumcenter);
            }
        }

        if (cell.vertices.empty()) {
            cell.vertices.push_back(cell.center);
        }

        dual.cells.push_back(std::move(cell));
    }

    return dual;
}

}  // namespace medial_axis_3d
