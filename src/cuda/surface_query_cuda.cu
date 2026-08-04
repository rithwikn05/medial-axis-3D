#include "surface_query_cuda.h"

#include <cuda_runtime_api.h>

#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "core/surface_query_cpu.h"
#include "cuda_buffer.cuh"
#include "surface_query_kernels.cuh"

namespace medial_axis_3d {
namespace {

constexpr unsigned int query_threads_per_block = 256;

cuda_detail::DeviceVec3 to_device_vec3(const Vec3& value) {
    return {value.x, value.y, value.z};
}

unsigned int query_block_count(std::size_t query_count) {
    const std::size_t blocks =
        (query_count + query_threads_per_block - 1) /
        query_threads_per_block;
    if (blocks > static_cast<std::size_t>(
                     std::numeric_limits<unsigned int>::max())) {
        throw std::length_error("CUDA surface-query batch is too large");
    }
    return static_cast<unsigned int>(blocks);
}

}  // namespace

bool cuda_surface_queries_available(std::string* reason) noexcept {
    int device_count = 0;
    const cudaError_t error = cudaGetDeviceCount(&device_count);
    if (error != cudaSuccess) {
        if (reason != nullptr) {
            *reason = cudaGetErrorString(error);
        }
        return false;
    }
    if (device_count <= 0) {
        if (reason != nullptr) {
            *reason = "No CUDA-capable device was found";
        }
        return false;
    }
    if (reason != nullptr) {
        reason->clear();
    }
    return true;
}

class CudaSurfaceQueryBackend::Impl {
public:
    explicit Impl(const Mesh& mesh)
        : cpu_fallback(mesh) {
        std::string unavailable_reason;
        if (!cuda_surface_queries_available(&unavailable_reason)) {
            throw std::runtime_error(
                "CUDA surface queries are unavailable: " +
                unavailable_reason
            );
        }

        std::vector<cuda_detail::DeviceTriangle> host_triangles;
        host_triangles.reserve(mesh.faces.size());
        for (const Triangle& face : mesh.faces) {
            for (int vertex : face.vertices) {
                if (vertex < 0 ||
                    vertex >= static_cast<int>(mesh.vertices.size())) {
                    throw std::invalid_argument(
                        "Surface mesh contains an invalid triangle index"
                    );
                }
            }
            host_triangles.push_back({
                to_device_vec3(
                    mesh.vertices[static_cast<std::size_t>(face.vertices[0])]
                ),
                to_device_vec3(
                    mesh.vertices[static_cast<std::size_t>(face.vertices[1])]
                ),
                to_device_vec3(
                    mesh.vertices[static_cast<std::size_t>(face.vertices[2])]
                )
            });
        }

        triangle_count = host_triangles.size();
        if (triangle_count > 0) {
            device_triangles.reserve(triangle_count);
            cuda_detail::check_cuda(
                cudaMemcpy(
                    device_triangles.data(),
                    host_triangles.data(),
                    triangle_count * sizeof(cuda_detail::DeviceTriangle),
                    cudaMemcpyHostToDevice
                ),
                "Copy surface triangles to CUDA"
            );
        }
    }

    CpuSurfaceQueryBackend cpu_fallback;
    std::size_t triangle_count{0};
    cuda_detail::DeviceBuffer<cuda_detail::DeviceTriangle> device_triangles;
    mutable cuda_detail::DeviceBuffer<cuda_detail::DeviceVec3> device_points;
    mutable cuda_detail::DeviceBuffer<cuda_detail::DeviceSegment> device_segments;
    mutable cuda_detail::DeviceBuffer<std::uint8_t> device_results;
};

CudaSurfaceQueryBackend::CudaSurfaceQueryBackend(const Mesh& mesh)
    : impl_(std::make_unique<Impl>(mesh)) {}

CudaSurfaceQueryBackend::~CudaSurfaceQueryBackend() = default;

CudaSurfaceQueryBackend::CudaSurfaceQueryBackend(
    CudaSurfaceQueryBackend&&
) noexcept = default;

CudaSurfaceQueryBackend& CudaSurfaceQueryBackend::operator=(
    CudaSurfaceQueryBackend&&
) noexcept = default;

const char* CudaSurfaceQueryBackend::name() const noexcept {
    return "cuda-hybrid";
}

void CudaSurfaceQueryBackend::points_inside(
    const std::vector<Vec3>& points,
    std::vector<SurfaceQueryClassification>& results
) const {
    results.resize(points.size());
    if (points.empty()) {
        return;
    }

    std::vector<cuda_detail::DeviceVec3> host_points;
    host_points.reserve(points.size());
    for (const Vec3& point : points) {
        host_points.push_back(to_device_vec3(point));
    }

    impl_->device_points.reserve(points.size());
    impl_->device_results.reserve(points.size());
    cuda_detail::check_cuda(
        cudaMemcpy(
            impl_->device_points.data(),
            host_points.data(),
            points.size() * sizeof(cuda_detail::DeviceVec3),
            cudaMemcpyHostToDevice
        ),
        "Copy point-containment queries to CUDA"
    );

    cuda_detail::points_inside_kernel<<<
        query_block_count(points.size()),
        query_threads_per_block
    >>>(
        impl_->device_triangles.data(),
        impl_->triangle_count,
        impl_->device_points.data(),
        points.size(),
        impl_->device_results.data()
    );
    cuda_detail::check_cuda(
        cudaGetLastError(),
        "Launch CUDA point-containment kernel"
    );
    cuda_detail::check_cuda(
        cudaMemcpy(
            results.data(),
            impl_->device_results.data(),
            results.size() * sizeof(SurfaceQueryClassification),
            cudaMemcpyDeviceToHost
        ),
        "Copy point-containment results from CUDA"
    );
}

void CudaSurfaceQueryBackend::segments_intersect(
    const std::vector<SurfaceSegmentQuery>& segments,
    std::vector<SurfaceQueryClassification>& results
) const {
    results.resize(segments.size());
    if (segments.empty()) {
        return;
    }

    std::vector<cuda_detail::DeviceSegment> host_segments;
    host_segments.reserve(segments.size());
    for (const SurfaceSegmentQuery& segment : segments) {
        host_segments.push_back({
            to_device_vec3(segment.start),
            to_device_vec3(segment.end)
        });
    }

    impl_->device_segments.reserve(segments.size());
    impl_->device_results.reserve(segments.size());
    cuda_detail::check_cuda(
        cudaMemcpy(
            impl_->device_segments.data(),
            host_segments.data(),
            segments.size() * sizeof(cuda_detail::DeviceSegment),
            cudaMemcpyHostToDevice
        ),
        "Copy segment-intersection queries to CUDA"
    );

    cuda_detail::segments_intersect_kernel<<<
        query_block_count(segments.size()),
        query_threads_per_block
    >>>(
        impl_->device_triangles.data(),
        impl_->triangle_count,
        impl_->device_segments.data(),
        segments.size(),
        impl_->device_results.data()
    );
    cuda_detail::check_cuda(
        cudaGetLastError(),
        "Launch CUDA segment-intersection kernel"
    );
    cuda_detail::check_cuda(
        cudaMemcpy(
            results.data(),
            impl_->device_results.data(),
            results.size() * sizeof(SurfaceQueryClassification),
            cudaMemcpyDeviceToHost
        ),
        "Copy segment-intersection results from CUDA"
    );
}

void CudaSurfaceQueryBackend::nearest_surface_contacts_batch(
    const std::vector<NearestSurfaceQuery>& queries,
    std::vector<std::vector<SurfaceContact>>& results
) const {
    impl_->cpu_fallback.nearest_surface_contacts_batch(queries, results);
}

}  // namespace medial_axis_3d
