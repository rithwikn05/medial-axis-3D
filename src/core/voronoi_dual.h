#pragma once

#include <vector>

#include "vec3.h"

namespace medial_axis_3d {

struct VoronoiCell {
    std::vector<Vec3> vertices;
};

struct VoronoiDual {
    std::vector<VoronoiCell> cells;
};

}  // namespace medial_axis_3d
