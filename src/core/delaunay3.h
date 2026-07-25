#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <utility>
#include <vector>

#include "vec3.h"

namespace medial_axis_3d {

struct Tetrahedron {
    std::array<int, 4> vertices{};
    std::array<int, 4> neighbors{}; // neighbor tetrahedron indices for each face (-1 if boundary)
    std::vector<std::size_t> parents; // history DAG: indices of tetrahedra that were replaced to create this one

    Tetrahedron() : vertices{}, neighbors{{-1, -1, -1, -1}}, parents() {}
    Tetrahedron(int a, int b, int c, int d) : vertices{a, b, c, d}, neighbors{{-1, -1, -1, -1}}, parents() {}
};

class Delaunay3 {
public:
    bool insert(const Vec3& point) {
        for (const auto& existing : points_) {
            if (existing == point) {
                return false;
            }
        }

        points_.push_back(point);

        if (points_.size() < 4) {
            return true;
        }

        if (points_.size() == 4) {
            tetrahedra_.push_back(Tetrahedron{0, 1, 2, 3});
            return true;
        }

        const int new_point_index = static_cast<int>(points_.size() - 1);
        const std::vector<std::size_t> cavity = cavity_indices(point);

        std::vector<Tetrahedron> kept_tetrahedra;
        kept_tetrahedra.reserve(tetrahedra_.size());
        for (std::size_t index = 0; index < tetrahedra_.size(); ++index) {
            if (std::find(cavity.begin(), cavity.end(), index) == cavity.end()) {
                kept_tetrahedra.push_back(tetrahedra_[index]);
            }
        }

        tetrahedra_ = std::move(kept_tetrahedra);

        if (cavity.empty()) {
            const std::size_t base = points_.size() - 4;
            tetrahedra_.push_back(Tetrahedron{
                static_cast<int>(base),
                static_cast<int>(base + 1),
                static_cast<int>(base + 2),
                static_cast<int>(base + 3)
            });
            rebuild_adjacency();
            return true;
        }

        const std::vector<std::array<int, 3>> boundary_faces = cavity_boundary_faces(cavity);
        for (const auto& face : boundary_faces) {
            tetrahedra_.push_back(Tetrahedron{face[0], face[1], face[2], new_point_index});
            tetrahedra_.back().parents = cavity;
        }

        // rebuild adjacency links after modification
        rebuild_adjacency();

        return true;
    }

    bool circumsphere_contains(std::size_t tetrahedron_index, const Vec3& point) const {
        if (tetrahedron_index >= tetrahedra_.size()) {
            return false;
        }

        const auto& tetrahedron = tetrahedra_[tetrahedron_index];
        if (tetrahedron.vertices[0] < 0 || tetrahedron.vertices[0] >= static_cast<int>(points_.size()) ||
            tetrahedron.vertices[1] < 0 || tetrahedron.vertices[1] >= static_cast<int>(points_.size()) ||
            tetrahedron.vertices[2] < 0 || tetrahedron.vertices[2] >= static_cast<int>(points_.size()) ||
            tetrahedron.vertices[3] < 0 || tetrahedron.vertices[3] >= static_cast<int>(points_.size())) {
            return false;
        }

        const Vec3& a = points_[tetrahedron.vertices[0]];
        const Vec3& b = points_[tetrahedron.vertices[1]];
        const Vec3& c = points_[tetrahedron.vertices[2]];
        const Vec3& d = points_[tetrahedron.vertices[3]];

        Vec3 center{};
        double radius_squared{};
        if (!compute_circumsphere(a, b, c, d, center, radius_squared)) {
            return false;
        }

        return squared_norm(point - center) <= radius_squared + 1e-12;
    }

    // compute the indices of tetrahedra whose circumspheres contain the given point
    std::vector<std::size_t> cavity_indices(const Vec3& point) const {
        std::vector<std::size_t> cavity;
        for (std::size_t index = 0; index < tetrahedra_.size(); ++index) {
            if (circumsphere_contains(index, point)) {
                cavity.push_back(index);
            }
        }
        return cavity;
    }

    const std::vector<Vec3>& points() const {
        return points_;
    }

    const std::vector<Tetrahedron>& tetrahedra() const {
        return tetrahedra_;
    }

    std::size_t point_count() const {
        return points_.size();
    }

    std::size_t tetrahedron_count() const {
        return tetrahedra_.size();
    }

private:
    std::vector<std::array<int, 3>> cavity_boundary_faces(const std::vector<std::size_t>& cavity) const {
        std::vector<std::pair<std::array<int, 3>, int>> face_counts;
        for (std::size_t tetrahedron_index : cavity) {
            const auto& tetrahedron = tetrahedra_[tetrahedron_index];
            std::vector<std::array<int, 3>> faces;
            faces.push_back({tetrahedron.vertices[0], tetrahedron.vertices[1], tetrahedron.vertices[2]});
            faces.push_back({tetrahedron.vertices[0], tetrahedron.vertices[2], tetrahedron.vertices[3]});
            faces.push_back({tetrahedron.vertices[0], tetrahedron.vertices[3], tetrahedron.vertices[1]});
            faces.push_back({tetrahedron.vertices[1], tetrahedron.vertices[2], tetrahedron.vertices[3]});

            for (const auto& face : faces) {
                const std::array<int, 3> normalized = normalize_face(face);
                auto it = std::find_if(face_counts.begin(), face_counts.end(), [&](const auto& entry) {
                    return entry.first == normalized;
                });

                if (it == face_counts.end()) {
                    face_counts.push_back({normalized, 1});
                } else {
                    ++it->second;
                }
            }
        }

        std::vector<std::array<int, 3>> boundary_faces;
        for (const auto& entry : face_counts) {
            if (entry.second == 1) {
                boundary_faces.push_back(entry.first);
            }
        }

        return boundary_faces;
    }

    static std::array<int, 3> normalize_face(const std::array<int, 3>& face) {
        std::array<int, 3> normalized = face;
        std::sort(normalized.begin(), normalized.end());
        return normalized;
    }

    void rebuild_adjacency() {
        // reset neighbors
        for (auto& tet : tetrahedra_) {
            tet.neighbors = {{-1, -1, -1, -1}};
        }

        // map normalized face -> (tet_index, face_id)
        std::map<std::array<int, 3>, std::pair<std::size_t, int>> face_owner;

        for (std::size_t t = 0; t < tetrahedra_.size(); ++t) {
            const auto& tet = tetrahedra_[t];
            const std::array<std::array<int, 3>, 4> faces = {{
                {tet.vertices[0], tet.vertices[1], tet.vertices[2]},
                {tet.vertices[0], tet.vertices[2], tet.vertices[3]},
                {tet.vertices[0], tet.vertices[3], tet.vertices[1]},
                {tet.vertices[1], tet.vertices[2], tet.vertices[3]}
            }};

            for (int f = 0; f < 4; ++f) {
                const auto normalized = normalize_face(faces[f]);
                auto it = face_owner.find(normalized);
                if (it == face_owner.end()) {
                    face_owner[normalized] = {t, f};
                } else {
                    const auto other = it->second;
                    tetrahedra_[t].neighbors[f] = static_cast<int>(other.first);
                    tetrahedra_[other.first].neighbors[other.second] = static_cast<int>(t);
                }
            }
        }
    }

    // compute circumsphere of tetrahedron defined by points a, b, c, d
    static bool compute_circumsphere(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d, Vec3& center, double& radius_squared) {
        const Vec3 ab = b - a;
        const Vec3 ac = c - a;
        const Vec3 ad = d - a;

        const std::array<std::array<double, 3>, 3> system = {{
            {2.0 * ab.x, 2.0 * ab.y, 2.0 * ab.z},
            {2.0 * ac.x, 2.0 * ac.y, 2.0 * ac.z},
            {2.0 * ad.x, 2.0 * ad.y, 2.0 * ad.z}
        }};
        const std::array<double, 3> rhs = {
            dot(b, b) - dot(a, a),
            dot(c, c) - dot(a, a),
            dot(d, d) - dot(a, a)
        };

        std::array<std::array<double, 3>, 3> matrix = system;
        std::array<double, 3> solution{};
        std::array<double, 3> rhs_values = rhs;

        for (int pivot_row = 0; pivot_row < 3; ++pivot_row) {
            int pivot_index = pivot_row;
            double pivot_value = std::fabs(matrix[pivot_row][pivot_row]);
            for (int row = pivot_row + 1; row < 3; ++row) {
                const double candidate = std::fabs(matrix[row][pivot_row]);
                if (candidate > pivot_value) {
                    pivot_value = candidate;
                    pivot_index = row;
                }
            }

            if (pivot_value < 1e-12) {
                return false;
            }

            if (pivot_index != pivot_row) {
                std::swap(matrix[pivot_row], matrix[pivot_index]);
                std::swap(rhs_values[pivot_row], rhs_values[pivot_index]);
            }

            const double pivot = matrix[pivot_row][pivot_row];
            for (int row = pivot_row + 1; row < 3; ++row) {
                const double factor = matrix[row][pivot_row] / pivot;
                for (int column = pivot_row; column < 3; ++column) {
                    matrix[row][column] -= factor * matrix[pivot_row][column];
                }
                rhs_values[row] -= factor * rhs_values[pivot_row];
            }
        }

        for (int row = 2; row >= 0; --row) {
            double sum = rhs_values[row];
            for (int column = row + 1; column < 3; ++column) {
                sum -= matrix[row][column] * solution[column];
            }
            solution[row] = sum / matrix[row][row];
        }

        center = Vec3(solution[0], solution[1], solution[2]);
        radius_squared = squared_norm(center - a);
        return true;
    }

    std::vector<Vec3> points_;
    std::vector<Tetrahedron> tetrahedra_;
};

}  // namespace medial_axis_3d
