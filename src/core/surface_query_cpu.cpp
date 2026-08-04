#include "surface_query_cpu.h"

#include <cstddef>

namespace medial_axis_3d {

CpuSurfaceQueryBackend::CpuSurfaceQueryBackend(const Mesh& mesh) noexcept
    : mesh_(mesh) {}

const char* CpuSurfaceQueryBackend::name() const noexcept {
    return "cpu-reference";
}

void CpuSurfaceQueryBackend::points_inside(
    const std::vector<Vec3>& points,
    std::vector<SurfaceQueryClassification>& results
) const {
    results.resize(points.size());
    for (std::size_t index = 0; index < points.size(); ++index) {
        results[index] = point_inside_mesh(mesh_, points[index])
            ? surface_query_true
            : surface_query_false;
    }
}

void CpuSurfaceQueryBackend::segments_intersect(
    const std::vector<SurfaceSegmentQuery>& segments,
    std::vector<SurfaceQueryClassification>& results
) const {
    results.resize(segments.size());
    for (std::size_t index = 0; index < segments.size(); ++index) {
        const SurfaceSegmentQuery& segment = segments[index];
        results[index] = segment_intersects_mesh_surface(
            mesh_,
            segment.start,
            segment.end
        )
            ? surface_query_true
            : surface_query_false;
    }
}

void CpuSurfaceQueryBackend::nearest_surface_contacts_batch(
    const std::vector<NearestSurfaceQuery>& queries,
    std::vector<std::vector<SurfaceContact>>& results
) const {
    results.resize(queries.size());
    for (std::size_t index = 0; index < queries.size(); ++index) {
        const NearestSurfaceQuery& query = queries[index];
        results[index] = medial_axis_3d::nearest_surface_contacts(
            mesh_,
            query.point,
            query.relative_tolerance
        );
    }
}

}  // namespace medial_axis_3d
