#pragma once

#include <cstdint>
#include <vector>

#include "medial_candidate.h"
#include "vec3.h"

namespace medial_axis_3d {

// A segment tested against the surface. Intersections at the two segment
// endpoints follow the same exclusion rules as
// segment_intersects_mesh_surface().
struct SurfaceSegmentQuery {
    Vec3 start{};
    Vec3 end{};
};

// A nearest-surface query. A positive relative tolerance may return multiple
// contacts near the minimum distance, matching nearest_surface_contacts().
struct NearestSurfaceQuery {
    Vec3 point{};
    double relative_tolerance{0.0};
};

using SurfaceQueryClassification = std::uint8_t;
constexpr SurfaceQueryClassification surface_query_false = 0;
constexpr SurfaceQueryClassification surface_query_true = 1;

// Backend-neutral batches for the surface geometry operations used by medial
// complex validation. Implementations may execute on the CPU, a CPU
// acceleration structure, or CUDA. Each output contains one entry per input
// query and preserves input order.
class SurfaceQueryBackend {
public:
    virtual ~SurfaceQueryBackend() = default;

    virtual const char* name() const noexcept = 0;

    virtual void points_inside(
        const std::vector<Vec3>& points,
        std::vector<SurfaceQueryClassification>& results
    ) const = 0;

    virtual void segments_intersect(
        const std::vector<SurfaceSegmentQuery>& segments,
        std::vector<SurfaceQueryClassification>& results
    ) const = 0;

    virtual void nearest_surface_contacts_batch(
        const std::vector<NearestSurfaceQuery>& queries,
        std::vector<std::vector<SurfaceContact>>& results
    ) const = 0;
};

}  // namespace medial_axis_3d
