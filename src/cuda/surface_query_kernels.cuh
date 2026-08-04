#pragma once

#include <cuda_runtime.h>

#include <cstddef>
#include <cstdint>

namespace medial_axis_3d::cuda_detail {

struct DeviceVec3 {
    double x;
    double y;
    double z;
};

struct DeviceTriangle {
    DeviceVec3 a;
    DeviceVec3 b;
    DeviceVec3 c;
};

struct DeviceSegment {
    DeviceVec3 start;
    DeviceVec3 end;
};

__device__ inline DeviceVec3 subtract(
    const DeviceVec3& a,
    const DeviceVec3& b
) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

__device__ inline double dot_product(
    const DeviceVec3& a,
    const DeviceVec3& b
) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

__device__ inline DeviceVec3 cross_product(
    const DeviceVec3& a,
    const DeviceVec3& b
) {
    return {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}

__device__ inline double squared_length(const DeviceVec3& value) {
    return dot_product(value, value);
}

__device__ inline double vector_length(const DeviceVec3& value) {
    return sqrt(squared_length(value));
}

__global__ void points_inside_kernel(
    const DeviceTriangle* triangles,
    std::size_t triangle_count,
    const DeviceVec3* points,
    std::size_t point_count,
    std::uint8_t* results
) {
    const std::size_t point_index =
        static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (point_index >= point_count) {
        return;
    }

    const DeviceVec3 point = points[point_index];
    double solid_angle_sum = 0.0;
    for (std::size_t triangle_index = 0;
         triangle_index < triangle_count;
         ++triangle_index) {
        const DeviceTriangle triangle = triangles[triangle_index];
        const DeviceVec3 a = subtract(triangle.a, point);
        const DeviceVec3 b = subtract(triangle.b, point);
        const DeviceVec3 c = subtract(triangle.c, point);

        const double length_a = vector_length(a);
        const double length_b = vector_length(b);
        const double length_c = vector_length(c);
        if (length_a == 0.0 || length_b == 0.0 || length_c == 0.0) {
            results[point_index] = 1;
            return;
        }

        const double numerator = dot_product(a, cross_product(b, c));
        const double denominator =
            length_a * length_b * length_c +
            dot_product(a, b) * length_c +
            dot_product(b, c) * length_a +
            dot_product(c, a) * length_b;
        solid_angle_sum += 2.0 * atan2(numerator, denominator);
    }

    constexpr double four_pi = 12.5663706143591729538;
    results[point_index] =
        fabs(solid_angle_sum / four_pi) > 0.5 ? 1 : 0;
}

__device__ inline bool segment_intersects_triangle(
    const DeviceSegment& segment,
    const DeviceTriangle& triangle
) {
    const DeviceVec3 direction = subtract(segment.end, segment.start);
    const DeviceVec3 edge_ab = subtract(triangle.b, triangle.a);
    const DeviceVec3 edge_ac = subtract(triangle.c, triangle.a);
    const DeviceVec3 p = cross_product(direction, edge_ac);
    const double determinant = dot_product(edge_ab, p);

    constexpr double determinant_epsilon = 1e-12;
    constexpr double barycentric_epsilon = 1e-10;
    constexpr double endpoint_epsilon = 1e-10;
    const double determinant_product =
        vector_length(direction) *
        vector_length(edge_ab) *
        vector_length(edge_ac);
    const double determinant_scale =
        determinant_product > 1e-30 ? determinant_product : 1e-30;
    if (fabs(determinant) <= determinant_epsilon * determinant_scale) {
        return false;
    }

    const double inverse_determinant = 1.0 / determinant;
    const DeviceVec3 from_a = subtract(segment.start, triangle.a);
    const double u = dot_product(from_a, p) * inverse_determinant;
    if (u < -barycentric_epsilon || u > 1.0 + barycentric_epsilon) {
        return false;
    }

    const DeviceVec3 q = cross_product(from_a, edge_ab);
    const double v = dot_product(direction, q) * inverse_determinant;
    if (v < -barycentric_epsilon ||
        u + v > 1.0 + barycentric_epsilon) {
        return false;
    }

    const double t = dot_product(edge_ac, q) * inverse_determinant;
    return t > endpoint_epsilon && t < 1.0 - endpoint_epsilon;
}

__global__ void segments_intersect_kernel(
    const DeviceTriangle* triangles,
    std::size_t triangle_count,
    const DeviceSegment* segments,
    std::size_t segment_count,
    std::uint8_t* results
) {
    const std::size_t segment_index =
        static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (segment_index >= segment_count) {
        return;
    }

    const DeviceSegment segment = segments[segment_index];
    const DeviceVec3 direction = subtract(segment.end, segment.start);
    if (squared_length(direction) <= 1e-30) {
        results[segment_index] = 0;
        return;
    }

    for (std::size_t triangle_index = 0;
         triangle_index < triangle_count;
         ++triangle_index) {
        if (segment_intersects_triangle(
                segment,
                triangles[triangle_index]
            )) {
            results[segment_index] = 1;
            return;
        }
    }
    results[segment_index] = 0;
}

}  // namespace medial_axis_3d::cuda_detail
