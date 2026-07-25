#pragma once

#include <limits>

namespace medial_axis_3d {

struct FilterSettings {
    double min_radius{0.0};
    double max_radius{std::numeric_limits<double>::infinity()};
};

}  // namespace medial_axis_3d
