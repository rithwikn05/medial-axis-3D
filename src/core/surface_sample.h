#pragma once

#include "vec3.h"

namespace medial_axis_3d {

struct SurfaceSample {
    Vec3 position{};
    Vec3 normal{};
};

}  // namespace medial_axis_3d
