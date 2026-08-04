#pragma once

#include <memory>
#include <string>

#include "core/mesh.h"
#include "core/surface_query_backend.h"

namespace medial_axis_3d {

// Reports whether the CUDA runtime can see at least one usable device. This
// does not allocate surface data or launch a kernel.
bool cuda_surface_queries_available(std::string* reason = nullptr) noexcept;

// CUDA implementation of the expensive classification queries. Nearest
// contacts currently use the CPU reference path so their multiple-contact and
// tolerance semantics remain unchanged.
class CudaSurfaceQueryBackend final : public SurfaceQueryBackend {
public:
    explicit CudaSurfaceQueryBackend(const Mesh& mesh);
    ~CudaSurfaceQueryBackend() override;

    CudaSurfaceQueryBackend(const CudaSurfaceQueryBackend&) = delete;
    CudaSurfaceQueryBackend& operator=(const CudaSurfaceQueryBackend&) = delete;
    CudaSurfaceQueryBackend(CudaSurfaceQueryBackend&&) noexcept;
    CudaSurfaceQueryBackend& operator=(CudaSurfaceQueryBackend&&) noexcept;

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
    class Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace medial_axis_3d
