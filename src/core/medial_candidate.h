#pragma once

#include "vec3.h"

namespace medial_axis_3d {

struct MedialCandidate {
    Vec3 center{};
    double radius{0.0};
};

}  // namespace medial_axis_3d
