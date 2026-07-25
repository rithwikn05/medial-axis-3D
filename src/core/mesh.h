#pragma once

#include <array>
#include <vector>

#include "vec3.h"

namespace medial_axis_3d {

struct Triangle {
    std::array<int, 3> vertices{};
};

struct Mesh {
    std::vector<Vec3> vertices;
    std::vector<Triangle> faces;
};

}  // namespace medial_axis_3d
