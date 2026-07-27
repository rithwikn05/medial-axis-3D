#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <utility>
#include <vector>

#include "predicates.h"
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
    bool build(const std::vector<Vec3>& points) {
        points_.clear();
        tetrahedra_.clear();

        for (const Vec3& point : points) {
            if (std::find(points_.begin(), points_.end(), point) != points_.end()) {
                points_.clear();
                return false;
            }
            points_.push_back(point);
        }

        if (points_.size() < 4) {
            return false;
        }

        Vec3 minimum = points_.front();
        Vec3 maximum = points_.front();
        for (const Vec3& point : points_) {
            minimum.x = std::min(minimum.x, point.x);
            minimum.y = std::min(minimum.y, point.y);
            minimum.z = std::min(minimum.z, point.z);
            maximum.x = std::max(maximum.x, point.x);
            maximum.y = std::max(maximum.y, point.y);
            maximum.z = std::max(maximum.z, point.z);
        }

        const Vec3 center = (minimum + maximum) / 2.0;
        const double extent = std::max({
            maximum.x - minimum.x,
            maximum.y - minimum.y,
            maximum.z - minimum.z
        });
        if (extent <= 0.0) {
            tetrahedra_.clear();
            return false;
        }

        std::vector<Vec3> working_points = points_;
        const double scale = 64.0 * extent;
        const int super_start = static_cast<int>(working_points.size());
        working_points.push_back(center + Vec3(-scale, -scale, -scale));
        working_points.push_back(center + Vec3( scale,  scale, -scale));
        working_points.push_back(center + Vec3( scale, -scale,  scale));
        working_points.push_back(center + Vec3(-scale,  scale,  scale));

        std::vector<Tetrahedron> working_tetrahedra;
        working_tetrahedra.emplace_back(super_start, super_start + 1, super_start + 2, super_start + 3);

        for (int point_index = 0; point_index < static_cast<int>(points_.size()); ++point_index) {
            std::vector<std::size_t> cavity;
            for (std::size_t tetrahedron_index = 0;
                 tetrahedron_index < working_tetrahedra.size();
                 ++tetrahedron_index) {
                const auto& tetrahedron = working_tetrahedra[tetrahedron_index];
                const Vec3& a = working_points[static_cast<std::size_t>(tetrahedron.vertices[0])];
                const Vec3& b = working_points[static_cast<std::size_t>(tetrahedron.vertices[1])];
                const Vec3& c = working_points[static_cast<std::size_t>(tetrahedron.vertices[2])];
                const Vec3& d = working_points[static_cast<std::size_t>(tetrahedron.vertices[3])];
                if (insphere_contains(a, b, c, d, working_points[static_cast<std::size_t>(point_index)])) {
                    cavity.push_back(tetrahedron_index);
                }
            }

            if (cavity.empty()) {
                tetrahedra_.clear();
                return false;
            }

            struct BoundaryFace {
                std::array<int, 3> vertices{};
                int count{0};
            };
            std::map<std::array<int, 3>, BoundaryFace> boundary;
            for (std::size_t tetrahedron_index : cavity) {
                const Tetrahedron& tetrahedron = working_tetrahedra[tetrahedron_index];
                const std::array<std::array<int, 3>, 4> faces{{
                    {{tetrahedron.vertices[0], tetrahedron.vertices[1], tetrahedron.vertices[2]}},
                    {{tetrahedron.vertices[0], tetrahedron.vertices[2], tetrahedron.vertices[3]}},
                    {{tetrahedron.vertices[0], tetrahedron.vertices[3], tetrahedron.vertices[1]}},
                    {{tetrahedron.vertices[1], tetrahedron.vertices[2], tetrahedron.vertices[3]}}
                }};
                for (const auto& face : faces) {
                    const auto key = normalize_face(face);
                    BoundaryFace& entry = boundary[key];
                    if (entry.count == 0) {
                        entry.vertices = face;
                    }
                    ++entry.count;
                }
            }

            std::vector<Tetrahedron> kept;
            kept.reserve(working_tetrahedra.size() + boundary.size());
            for (std::size_t tetrahedron_index = 0;
                 tetrahedron_index < working_tetrahedra.size();
                 ++tetrahedron_index) {
                if (std::find(cavity.begin(), cavity.end(), tetrahedron_index) == cavity.end()) {
                    kept.push_back(working_tetrahedra[tetrahedron_index]);
                }
            }

            for (const auto& [key, face] : boundary) {
                (void)key;
                if (face.count != 1) {
                    continue;
                }
                const Vec3& a = working_points[static_cast<std::size_t>(face.vertices[0])];
                const Vec3& b = working_points[static_cast<std::size_t>(face.vertices[1])];
                const Vec3& c = working_points[static_cast<std::size_t>(face.vertices[2])];
                const Vec3& d = working_points[static_cast<std::size_t>(point_index)];
                if (orient3d(a, b, c, d) != 0.0) {
                    kept.emplace_back(face.vertices[0], face.vertices[1], face.vertices[2], point_index);
                }
            }
            working_tetrahedra = std::move(kept);
        }

        for (const Tetrahedron& tetrahedron : working_tetrahedra) {
            bool uses_super_vertex = false;
            for (int vertex : tetrahedron.vertices) {
                if (vertex >= super_start) {
                    uses_super_vertex = true;
                    break;
                }
            }
            if (!uses_super_vertex) {
                tetrahedra_.push_back(tetrahedron);
            }
        }

        rebuild_adjacency();
        return !tetrahedra_.empty() && is_well_formed();
    }

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

        // The local flip scaffold is only reliable for the first cavity split.
        // Rebuild larger point sets with the robust batch Bowyer-Watson path.
        if (points_.size() > 5) {
            const std::vector<Vec3> all_points = points_;
            return build(all_points);
        }

        const int new_point_index = static_cast<int>(points_.size() - 1);
        const std::vector<std::size_t> cavity = cavity_indices(point);
        const std::vector<std::array<int, 3>> boundary_faces =
            cavity_boundary_faces(cavity);

        std::vector<Tetrahedron> kept_tetrahedra;
        kept_tetrahedra.reserve(tetrahedra_.size());
        for (std::size_t index = 0; index < tetrahedra_.size(); ++index) {
            if (std::find(cavity.begin(), cavity.end(), index) == cavity.end()) {
                kept_tetrahedra.push_back(tetrahedra_[index]);
            }
        }

        tetrahedra_ = std::move(kept_tetrahedra);

        if (cavity.empty()) {
            const std::vector<Vec3> all_points = points_;
            return build(all_points);
        }

        for (const auto& face : boundary_faces) {
            tetrahedra_.push_back(Tetrahedron{face[0], face[1], face[2], new_point_index});
            tetrahedra_.back().parents = cavity;
        }

        // rebuild adjacency links after modification
        rebuild_adjacency();
        repair_local_delaunay();

        return true;
    }

    bool circumcenter(std::size_t tetrahedron_index, Vec3& center) const {
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

        double radius_squared{};
        if (!compute_circumsphere(a, b, c, d, center, radius_squared)) {
            return false;
        }

        return true;
    }

    static std::array<int, 3> face_from_index(const Tetrahedron& tet, int face_index) {
        switch (face_index) {
            case 0: return {tet.vertices[0], tet.vertices[1], tet.vertices[2]};
            case 1: return {tet.vertices[0], tet.vertices[2], tet.vertices[3]};
            case 2: return {tet.vertices[0], tet.vertices[3], tet.vertices[1]};
            case 3: return {tet.vertices[1], tet.vertices[2], tet.vertices[3]};
            default: return {tet.vertices[0], tet.vertices[1], tet.vertices[2]};
        }
    }

    static double signed_volume(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
        return dot(b - a, cross(c - a, d - a));
    }

    static bool points_are_on_opposite_sides(const std::array<int, 3>& face,
                                             const Vec3& first,
                                             const Vec3& second,
                                             const std::vector<Vec3>& points) {
        const double side_first = orient3d(points[face[0]], points[face[1]], points[face[2]], first);
        const double side_second = orient3d(points[face[0]], points[face[1]], points[face[2]], second);
        if (std::fabs(side_first) < 1e-12 || std::fabs(side_second) < 1e-12) {
            return false;
        }
        return side_first * side_second < 0.0;
    }

    static int face_index_for(const Tetrahedron& tetrahedron, const std::array<int, 3>& face) {
        for (int candidate = 0; candidate < 4; ++candidate) {
            if (normalize_face(face_from_index(tetrahedron, candidate)) == normalize_face(face)) {
                return candidate;
            }
        }
        return -1;
    }

    bool flip_if_needed(std::size_t tetrahedron_a, std::size_t tetrahedron_b, int face_index) {
        if (tetrahedron_a >= tetrahedra_.size() || tetrahedron_b >= tetrahedra_.size()) {
            return false;
        }

        if (tetrahedron_a == tetrahedron_b) {
            return false;
        }

        const auto& tet_a = tetrahedra_[tetrahedron_a];
        const auto& tet_b = tetrahedra_[tetrahedron_b];
        const auto shared_face = face_from_index(tet_a, face_index);

        for (int vertex : shared_face) {
            if (std::find(tet_b.vertices.begin(), tet_b.vertices.end(), vertex) == tet_b.vertices.end()) {
                return false;
            }
        }

        const int opposite_a = opposite_vertex(tet_a, shared_face);
        const int opposite_b = opposite_vertex(tet_b, shared_face);
        if (opposite_a < 0 || opposite_b < 0) {
            return false;
        }

        if (!points_are_on_opposite_sides(shared_face, points_[opposite_a], points_[opposite_b], points_)) {
            return false;
        }

        if (!circumsphere_contains(tetrahedron_a, points_[opposite_b]) &&
            !circumsphere_contains(tetrahedron_b, points_[opposite_a])) {
            return false;
        }

        std::vector<Tetrahedron> updated_tetrahedra;
        updated_tetrahedra.reserve(tetrahedra_.size() + 1);
        for (std::size_t index = 0; index < tetrahedra_.size(); ++index) {
            if (index == tetrahedron_a || index == tetrahedron_b) {
                continue;
            }
            updated_tetrahedra.push_back(tetrahedra_[index]);
        }

        Tetrahedron t0{shared_face[0], shared_face[1], opposite_a, opposite_b};
        Tetrahedron t1{shared_face[1], shared_face[2], opposite_a, opposite_b};
        Tetrahedron t2{shared_face[2], shared_face[0], opposite_a, opposite_b};
        t0.parents = {tetrahedron_a, tetrahedron_b};
        t1.parents = {tetrahedron_a, tetrahedron_b};
        t2.parents = {tetrahedron_a, tetrahedron_b};

        updated_tetrahedra.push_back(std::move(t0));
        updated_tetrahedra.push_back(std::move(t1));
        updated_tetrahedra.push_back(std::move(t2));

        const std::size_t new_t0 = updated_tetrahedra.size() - 3;
        const std::size_t new_t1 = updated_tetrahedra.size() - 2;
        const std::size_t new_t2 = updated_tetrahedra.size() - 1;

        for (int old_face_index = 0; old_face_index < 4; ++old_face_index) {
            if (old_face_index == face_index) {
                continue;
            }

            const auto old_face = face_from_index(tet_a, old_face_index);
            const int new_face = face_index_for(updated_tetrahedra[new_t0], old_face);
            if (new_face >= 0) {
                updated_tetrahedra[new_t0].neighbors[new_face] = tet_a.neighbors[old_face_index];
            }
        }

        for (int old_face_index = 0; old_face_index < 4; ++old_face_index) {
            if (old_face_index == face_index) {
                continue;
            }

            const auto old_face = face_from_index(tet_b, old_face_index);
            const int new_face = face_index_for(updated_tetrahedra[new_t0], old_face);
            if (new_face >= 0) {
                updated_tetrahedra[new_t0].neighbors[new_face] = tet_b.neighbors[old_face_index];
            }
        }

        const auto internal_face_t0_t1 = std::array<int, 3>{shared_face[1], opposite_a, opposite_b};
        const auto internal_face_t1_t2 = std::array<int, 3>{shared_face[2], opposite_a, opposite_b};
        const auto internal_face_t2_t0 = std::array<int, 3>{shared_face[0], opposite_a, opposite_b};

        const int face_t0_t1 = face_index_for(updated_tetrahedra[new_t0], internal_face_t0_t1);
        const int face_t1_t0 = face_index_for(updated_tetrahedra[new_t1], internal_face_t0_t1);
        if (face_t0_t1 >= 0 && face_t1_t0 >= 0) {
            updated_tetrahedra[new_t0].neighbors[face_t0_t1] = static_cast<int>(new_t1);
            updated_tetrahedra[new_t1].neighbors[face_t1_t0] = static_cast<int>(new_t0);
        }

        const int face_t1_t2 = face_index_for(updated_tetrahedra[new_t1], internal_face_t1_t2);
        const int face_t2_t1 = face_index_for(updated_tetrahedra[new_t2], internal_face_t1_t2);
        if (face_t1_t2 >= 0 && face_t2_t1 >= 0) {
            updated_tetrahedra[new_t1].neighbors[face_t1_t2] = static_cast<int>(new_t2);
            updated_tetrahedra[new_t2].neighbors[face_t2_t1] = static_cast<int>(new_t1);
        }

        const int face_t2_t0 = face_index_for(updated_tetrahedra[new_t2], internal_face_t2_t0);
        const int face_t0_t2 = face_index_for(updated_tetrahedra[new_t0], internal_face_t2_t0);
        if (face_t2_t0 >= 0 && face_t0_t2 >= 0) {
            updated_tetrahedra[new_t2].neighbors[face_t2_t0] = static_cast<int>(new_t0);
            updated_tetrahedra[new_t0].neighbors[face_t0_t2] = static_cast<int>(new_t2);
        }

        tetrahedra_ = std::move(updated_tetrahedra);
        rebuild_adjacency();
        return true;
    }

    bool is_well_formed() const {
        for (const auto& tetrahedron : tetrahedra_) {
            for (int vertex : tetrahedron.vertices) {
                if (vertex < 0 || vertex >= static_cast<int>(points_.size())) {
                    return false;
                }
            }
            if (tetrahedron.vertices[0] == tetrahedron.vertices[1] ||
                tetrahedron.vertices[0] == tetrahedron.vertices[2] ||
                tetrahedron.vertices[0] == tetrahedron.vertices[3] ||
                tetrahedron.vertices[1] == tetrahedron.vertices[2] ||
                tetrahedron.vertices[1] == tetrahedron.vertices[3] ||
                tetrahedron.vertices[2] == tetrahedron.vertices[3]) {
                return false;
            }
        }
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

        return insphere_contains(a, b, c, d, point);
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

    void repair_local_delaunay() {
        for (int pass = 0; pass < 8; ++pass) {
            rebuild_adjacency();
            bool changed = false;
            for (std::size_t tetrahedron_index = 0; tetrahedron_index < tetrahedra_.size(); ++tetrahedron_index) {
                const auto& tet = tetrahedra_[tetrahedron_index];
                for (int face = 0; face < 4; ++face) {
                    const int neighbor = tet.neighbors[face];
                    if (neighbor < 0 || neighbor <= static_cast<int>(tetrahedron_index)) {
                        continue;
                    }

                    if (flip_if_needed(tetrahedron_index, static_cast<std::size_t>(neighbor), face)) {
                        changed = true;
                        break;
                    }
                }
                if (changed) {
                    break;
                }
            }
            if (!changed) {
                break;
            }
        }
        rebuild_adjacency();
    }

    static int opposite_vertex(const Tetrahedron& tetrahedron, const std::array<int, 3>& face) {
        for (int vertex : tetrahedron.vertices) {
            if (std::find(face.begin(), face.end(), vertex) == face.end()) {
                return vertex;
            }
        }
        return -1;
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
