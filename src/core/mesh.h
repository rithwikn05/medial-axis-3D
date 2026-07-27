#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <queue>
#include <string>
#include <utility>
#include <vector>

#include "vec3.h"

namespace medial_axis_3d {

struct Triangle {
    std::array<int, 3> vertices{};
};

struct Mesh {
    std::vector<Vec3> vertices;
    std::vector<Triangle> faces;
    std::vector<int> face_boundary_markers;
    std::vector<Vec3> face_normals;
    std::vector<Vec3> vertex_normals;
};

inline bool orient_and_analyze_closed_mesh(Mesh& mesh, std::string& error) {
    error.clear();
    mesh.face_normals.clear();
    mesh.vertex_normals.clear();

    if (mesh.vertices.size() < 4 || mesh.faces.size() < 4) {
        error = "A closed 3D surface requires at least four vertices and four triangles.";
        return false;
    }

    struct EdgeUse {
        std::size_t face{0};
        bool follows_sorted_direction{false};
    };
    using Edge = std::array<int, 2>;
    std::map<Edge, std::vector<EdgeUse>> edge_uses;

    for (std::size_t face_index = 0; face_index < mesh.faces.size(); ++face_index) {
        const auto& vertices = mesh.faces[face_index].vertices;
        for (int vertex : vertices) {
            if (vertex < 0 || vertex >= static_cast<int>(mesh.vertices.size())) {
                error = "Surface triangle " + std::to_string(face_index) +
                        " contains an invalid vertex index.";
                return false;
            }
        }
        if (vertices[0] == vertices[1] ||
            vertices[1] == vertices[2] ||
            vertices[2] == vertices[0]) {
            error = "Surface triangle " + std::to_string(face_index) +
                    " repeats a vertex.";
            return false;
        }

        const Vec3& a = mesh.vertices[static_cast<std::size_t>(vertices[0])];
        const Vec3& b = mesh.vertices[static_cast<std::size_t>(vertices[1])];
        const Vec3& c = mesh.vertices[static_cast<std::size_t>(vertices[2])];
        if (squared_norm(cross(b - a, c - a)) < 1e-24) {
            error = "Surface triangle " + std::to_string(face_index) +
                    " is degenerate.";
            return false;
        }

        for (int edge_index = 0; edge_index < 3; ++edge_index) {
            const int first = vertices[static_cast<std::size_t>(edge_index)];
            const int second = vertices[static_cast<std::size_t>((edge_index + 1) % 3)];
            const Edge edge = first < second
                ? Edge{first, second}
                : Edge{second, first};
            edge_uses[edge].push_back({face_index, first < second});
        }
    }

    std::vector<std::vector<std::pair<std::size_t, bool>>> adjacency(mesh.faces.size());
    for (const auto& [edge, uses] : edge_uses) {
        if (uses.size() != 2) {
            error = "Surface is not watertight/manifold at edge (" +
                    std::to_string(edge[0]) + ", " + std::to_string(edge[1]) +
                    "): expected two incident triangles, found " +
                    std::to_string(uses.size()) + ".";
            return false;
        }

        const bool same_direction =
            uses[0].follows_sorted_direction == uses[1].follows_sorted_direction;
        adjacency[uses[0].face].push_back({uses[1].face, same_direction});
        adjacency[uses[1].face].push_back({uses[0].face, same_direction});
    }

    std::vector<int> flip_state(mesh.faces.size(), -1);
    std::vector<std::vector<std::size_t>> components;
    for (std::size_t seed = 0; seed < mesh.faces.size(); ++seed) {
        if (flip_state[seed] != -1) {
            continue;
        }

        components.emplace_back();
        std::queue<std::size_t> pending;
        pending.push(seed);
        flip_state[seed] = 0;

        while (!pending.empty()) {
            const std::size_t face = pending.front();
            pending.pop();
            components.back().push_back(face);

            for (const auto& [neighbor, same_direction] : adjacency[face]) {
                const int required_state =
                    flip_state[face] ^ (same_direction ? 1 : 0);
                if (flip_state[neighbor] == -1) {
                    flip_state[neighbor] = required_state;
                    pending.push(neighbor);
                } else if (flip_state[neighbor] != required_state) {
                    error = "Surface triangles cannot be oriented consistently.";
                    return false;
                }
            }
        }
    }

    for (std::size_t face_index = 0; face_index < mesh.faces.size(); ++face_index) {
        if (flip_state[face_index] == 1) {
            std::swap(
                mesh.faces[face_index].vertices[1],
                mesh.faces[face_index].vertices[2]
            );
        }
    }

    for (const auto& component : components) {
        double signed_volume_times_six = 0.0;
        for (std::size_t face_index : component) {
            const auto& face = mesh.faces[face_index].vertices;
            const Vec3& a = mesh.vertices[static_cast<std::size_t>(face[0])];
            const Vec3& b = mesh.vertices[static_cast<std::size_t>(face[1])];
            const Vec3& c = mesh.vertices[static_cast<std::size_t>(face[2])];
            signed_volume_times_six += dot(a, cross(b, c));
        }
        if (std::fabs(signed_volume_times_six) < 1e-18) {
            error = "Closed surface component has near-zero signed volume.";
            return false;
        }
        if (signed_volume_times_six < 0.0) {
            for (std::size_t face_index : component) {
                std::swap(
                    mesh.faces[face_index].vertices[1],
                    mesh.faces[face_index].vertices[2]
                );
            }
        }
    }

    mesh.face_normals.reserve(mesh.faces.size());
    mesh.vertex_normals.assign(mesh.vertices.size(), Vec3{});
    for (const Triangle& triangle : mesh.faces) {
        const Vec3& a =
            mesh.vertices[static_cast<std::size_t>(triangle.vertices[0])];
        const Vec3& b =
            mesh.vertices[static_cast<std::size_t>(triangle.vertices[1])];
        const Vec3& c =
            mesh.vertices[static_cast<std::size_t>(triangle.vertices[2])];
        const Vec3 area_normal = cross(b - a, c - a);
        mesh.face_normals.push_back(normalized(area_normal));
        for (int vertex : triangle.vertices) {
            mesh.vertex_normals[static_cast<std::size_t>(vertex)] =
                mesh.vertex_normals[static_cast<std::size_t>(vertex)] + area_normal;
        }
    }
    for (Vec3& normal : mesh.vertex_normals) {
        normal = normalized(normal);
    }

    return true;
}

inline double mesh_winding_number(const Mesh& mesh, const Vec3& point) {
    double solid_angle_sum = 0.0;
    for (const Triangle& triangle : mesh.faces) {
        const Vec3 a =
            mesh.vertices[static_cast<std::size_t>(triangle.vertices[0])] - point;
        const Vec3 b =
            mesh.vertices[static_cast<std::size_t>(triangle.vertices[1])] - point;
        const Vec3 c =
            mesh.vertices[static_cast<std::size_t>(triangle.vertices[2])] - point;

        const double length_a = norm(a);
        const double length_b = norm(b);
        const double length_c = norm(c);
        if (length_a == 0.0 || length_b == 0.0 || length_c == 0.0) {
            return 1.0;
        }

        const double numerator = dot(a, cross(b, c));
        const double denominator =
            length_a * length_b * length_c +
            dot(a, b) * length_c +
            dot(b, c) * length_a +
            dot(c, a) * length_b;
        solid_angle_sum += 2.0 * std::atan2(numerator, denominator);
    }

    constexpr double four_pi = 12.5663706143591729538;
    return solid_angle_sum / four_pi;
}

inline bool point_inside_mesh(const Mesh& mesh, const Vec3& point) {
    return std::fabs(mesh_winding_number(mesh, point)) > 0.5;
}

inline bool segment_intersects_mesh_surface(
    const Mesh& mesh,
    const Vec3& start,
    const Vec3& end) {
    const Vec3 direction = end - start;
    if (squared_norm(direction) <= 1e-30) {
        return false;
    }

    // Moller-Trumbore segment/triangle intersection. Intersections strictly
    // between the endpoints mean that a chord between two interior points
    // leaves or touches the solid boundary.
    constexpr double determinant_epsilon = 1e-12;
    constexpr double barycentric_epsilon = 1e-10;
    constexpr double endpoint_epsilon = 1e-10;
    for (const Triangle& triangle : mesh.faces) {
        const Vec3& a =
            mesh.vertices[static_cast<std::size_t>(triangle.vertices[0])];
        const Vec3& b =
            mesh.vertices[static_cast<std::size_t>(triangle.vertices[1])];
        const Vec3& c =
            mesh.vertices[static_cast<std::size_t>(triangle.vertices[2])];
        const Vec3 edge_ab = b - a;
        const Vec3 edge_ac = c - a;
        const Vec3 p = cross(direction, edge_ac);
        const double determinant = dot(edge_ab, p);
        const double determinant_scale = std::max(
            1e-30,
            norm(direction) * norm(edge_ab) * norm(edge_ac)
        );
        if (std::fabs(determinant) <=
            determinant_epsilon * determinant_scale) {
            continue;
        }

        const double inverse_determinant = 1.0 / determinant;
        const Vec3 from_a = start - a;
        const double u = dot(from_a, p) * inverse_determinant;
        if (u < -barycentric_epsilon ||
            u > 1.0 + barycentric_epsilon) {
            continue;
        }
        const Vec3 q = cross(from_a, edge_ab);
        const double v = dot(direction, q) * inverse_determinant;
        if (v < -barycentric_epsilon ||
            u + v > 1.0 + barycentric_epsilon) {
            continue;
        }
        const double t = dot(edge_ac, q) * inverse_determinant;
        if (t > endpoint_epsilon &&
            t < 1.0 - endpoint_epsilon) {
            return true;
        }
    }
    return false;
}

inline bool segment_inside_mesh(
    const Mesh& mesh,
    const Vec3& start,
    const Vec3& end) {
    return point_inside_mesh(mesh, start) &&
           point_inside_mesh(mesh, end) &&
           point_inside_mesh(mesh, (start + end) / 2.0) &&
           !segment_intersects_mesh_surface(mesh, start, end);
}

}  // namespace medial_axis_3d
