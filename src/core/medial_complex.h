#pragma once

#include <vector>

#include "medial_candidate.h"

namespace medial_axis_3d {

struct MedialComplex {
    std::vector<MedialCandidate> candidates;
};

}  // namespace medial_axis_3d
