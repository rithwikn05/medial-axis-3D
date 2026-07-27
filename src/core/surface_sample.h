#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <vector>

#include "mesh.h"
#include "vec3.h"

namespace medial_axis_3d {

struct SurfaceSample {
    Vec3 position{};
    Vec3 normal{};
    std::size_t source_triangle{0};
    std::array<double, 3> barycentric{{0.0, 0.0, 0.0}};
    int boundary_marker{0};
    int source_vertex{-1};
};

struct SurfaceSamplingOptions {
    std::size_t target_sample_count{0};
    bool include_mesh_vertices{true};
    bool interpolate_vertex_normals{false};
};

namespace detail {

inline double radical_inverse_base2(std::size_t value) {
    double result = 0.0;
    double place = 0.5;
    while (value != 0) {
        if ((value & 1U) != 0U) {
            result += place;
        }
        value >>= 1U;
        place *= 0.5;
    }
    return result;
}

inline int face_boundary_marker(const Mesh& mesh, std::size_t face_index) {
    return face_index < mesh.face_boundary_markers.size()
        ? mesh.face_boundary_markers[face_index]
        : 0;
}

}  // namespace detail

inline std::vector<SurfaceSample> sample_surface(
    const Mesh& mesh,
    const SurfaceSamplingOptions& options) {
    std::vector<SurfaceSample> samples;
    if (mesh.faces.empty() ||
        mesh.face_normals.size() != mesh.faces.size() ||
        mesh.vertex_normals.size() != mesh.vertices.size()) {
        return samples;
    }

    std::vector<int> first_incident_face(mesh.vertices.size(), -1);
    std::vector<int> corner_in_first_face(mesh.vertices.size(), -1);
    for (std::size_t face_index = 0; face_index < mesh.faces.size(); ++face_index) {
        for (int corner = 0; corner < 3; ++corner) {
            const int vertex =
                mesh.faces[face_index].vertices[static_cast<std::size_t>(corner)];
            if (first_incident_face[static_cast<std::size_t>(vertex)] < 0) {
                first_incident_face[static_cast<std::size_t>(vertex)] =
                    static_cast<int>(face_index);
                corner_in_first_face[static_cast<std::size_t>(vertex)] = corner;
            }
        }
    }

    if (options.include_mesh_vertices) {
        for (std::size_t vertex = 0; vertex < mesh.vertices.size(); ++vertex) {
            const int face_index = first_incident_face[vertex];
            if (face_index < 0) {
                continue;
            }

            SurfaceSample sample;
            sample.position = mesh.vertices[vertex];
            sample.normal = mesh.vertex_normals[vertex];
            sample.source_triangle = static_cast<std::size_t>(face_index);
            sample.source_vertex = static_cast<int>(vertex);
            sample.barycentric[static_cast<std::size_t>(
                corner_in_first_face[vertex]
            )] = 1.0;
            sample.boundary_marker = detail::face_boundary_marker(
                mesh,
                sample.source_triangle
            );
            samples.push_back(sample);
        }
    }

    const std::size_t target = std::max(
        options.target_sample_count,
        samples.size()
    );
    if (target == samples.size()) {
        return samples;
    }

    std::vector<double> triangle_areas(mesh.faces.size(), 0.0);
    double total_area = 0.0;
    for (std::size_t face_index = 0; face_index < mesh.faces.size(); ++face_index) {
        const auto& face = mesh.faces[face_index].vertices;
        const Vec3& a = mesh.vertices[static_cast<std::size_t>(face[0])];
        const Vec3& b = mesh.vertices[static_cast<std::size_t>(face[1])];
        const Vec3& c = mesh.vertices[static_cast<std::size_t>(face[2])];
        triangle_areas[face_index] = 0.5 * norm(cross(b - a, c - a));
        total_area += triangle_areas[face_index];
    }
    if (total_area <= 0.0) {
        return {};
    }

    const std::size_t generated_count = target - samples.size();
    std::vector<std::size_t> samples_per_face(mesh.faces.size(), 0);
    std::vector<double> fractional_remainders(mesh.faces.size(), 0.0);
    std::size_t allocated = 0;
    for (std::size_t face_index = 0; face_index < mesh.faces.size(); ++face_index) {
        const double exact =
            static_cast<double>(generated_count) *
            triangle_areas[face_index] / total_area;
        const std::size_t count = static_cast<std::size_t>(std::floor(exact));
        samples_per_face[face_index] = count;
        fractional_remainders[face_index] = exact - static_cast<double>(count);
        allocated += count;
    }

    std::vector<std::size_t> remainder_order(mesh.faces.size());
    std::iota(remainder_order.begin(), remainder_order.end(), 0);
    std::stable_sort(
        remainder_order.begin(),
        remainder_order.end(),
        [&](std::size_t a, std::size_t b) {
            return fractional_remainders[a] > fractional_remainders[b];
        }
    );
    for (std::size_t i = allocated; i < generated_count; ++i) {
        ++samples_per_face[remainder_order[i - allocated]];
    }

    samples.reserve(target);
    std::size_t sequence_index = 1;
    for (std::size_t face_index = 0; face_index < mesh.faces.size(); ++face_index) {
        const std::size_t count = samples_per_face[face_index];
        if (count == 0) {
            continue;
        }

        const auto& face = mesh.faces[face_index].vertices;
        const Vec3& a = mesh.vertices[static_cast<std::size_t>(face[0])];
        const Vec3& b = mesh.vertices[static_cast<std::size_t>(face[1])];
        const Vec3& c = mesh.vertices[static_cast<std::size_t>(face[2])];
        const Vec3& normal_a =
            mesh.vertex_normals[static_cast<std::size_t>(face[0])];
        const Vec3& normal_b =
            mesh.vertex_normals[static_cast<std::size_t>(face[1])];
        const Vec3& normal_c =
            mesh.vertex_normals[static_cast<std::size_t>(face[2])];

        for (std::size_t local_index = 0; local_index < count; ++local_index) {
            const double u =
                (static_cast<double>(local_index) + 0.5) /
                static_cast<double>(count);
            const double v = detail::radical_inverse_base2(sequence_index++);
            const double root_u = std::sqrt(u);

            SurfaceSample sample;
            sample.source_triangle = face_index;
            sample.barycentric = {
                1.0 - root_u,
                root_u * (1.0 - v),
                root_u * v
            };
            sample.position =
                a * sample.barycentric[0] +
                b * sample.barycentric[1] +
                c * sample.barycentric[2];
            sample.normal = options.interpolate_vertex_normals
                ? normalized(
                    normal_a * sample.barycentric[0] +
                    normal_b * sample.barycentric[1] +
                    normal_c * sample.barycentric[2]
                )
                : mesh.face_normals[face_index];
            if (squared_norm(sample.normal) == 0.0) {
                sample.normal = mesh.face_normals[face_index];
            }
            sample.boundary_marker =
                detail::face_boundary_marker(mesh, face_index);
            samples.push_back(sample);
        }
    }

    return samples;
}

}  // namespace medial_axis_3d
