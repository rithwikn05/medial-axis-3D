#pragma once

#include "mesh.h"
#include "surface_query_backend.h"

namespace medial_axis_3d {

// Reference implementation that deliberately reuses the existing scalar
// surface-query functions. The referenced mesh must outlive this backend.
class CpuSurfaceQueryBackend final : public SurfaceQueryBackend {
public:
    explicit CpuSurfaceQueryBackend(const Mesh& mesh) noexcept;

    const char* name() const noexcept override;

    void points_inside(
        const std::vector<Vec3>& points,
        std::vector<SurfaceQueryClassification>& results
    ) const override;

    void segments_intersect(
        const std::vector<SurfaceSegmentQuery>& segments,
        std::vector<SurfaceQueryClassification>& results
    ) const override;

    void nearest_surface_contacts_batch(
        const std::vector<NearestSurfaceQuery>& queries,
        std::vector<std::vector<SurfaceContact>>& results
    ) const override;

private:
    const Mesh& mesh_;
};

}  // namespace medial_axis_3d
