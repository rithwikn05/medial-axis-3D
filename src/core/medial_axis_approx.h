#pragma once

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include "delaunay3.h"
#include "mesh.h"
#include "vec3.h"

namespace medial_axis_3d {

struct BoundarySurface {
    std::vector<std::array<std::size_t, 3>> faces;
};

struct MedialAxisApproximation {
    std::vector<Vec3> centers;
    std::vector<double> radii;
    std::vector<std::array<std::size_t, 2>> edges;
    std::vector<std::size_t> source_tetrahedra;
    std::size_t rejected_exterior_edge_count{0};
};

struct VoronoiFaceCandidate {
    std::array<int, 2> contact_sites{};
    std::vector<Vec3> vertices;
    // Parallel to vertices. Retaining the generating tetrahedra lets later
    // stages weld the same Voronoi vertex across independently built faces.
    std::vector<std::size_t> source_tetrahedra;
    double contact_angle_degrees{0.0};
};

struct MedialSheetApproximation {
    std::vector<Vec3> vertices;
    std::vector<std::array<std::size_t, 3>> triangles;
    std::vector<double> triangle_contact_angles;
    std::size_t polygon_count{0};
};

struct ValidatedMedialComplexProfile {
    double candidate_generation_seconds{0.0};
    double point_containment_winding_seconds{0.0};
    double surface_distance_queries_seconds{0.0};
    double segment_intersections_seconds{0.0};
    double pole_support_propagation_seconds{0.0};
    double topology_completion_seconds{0.0};
};

namespace detail {

using MedialProfileClock = std::chrono::steady_clock;

inline double medial_profile_elapsed_seconds(
    const MedialProfileClock::time_point& start) {
    return std::chrono::duration<double>(
        MedialProfileClock::now() - start
    ).count();
}

}  // namespace detail

inline bool polygon_fan_inside_mesh(
    const Mesh& mesh,
    const std::vector<Vec3>& polygon,
    ValidatedMedialComplexProfile* profile = nullptr) {
    if (polygon.size() < 3) {
        return false;
    }

    const auto point_is_inside = [&](const Vec3& point) {
        if (profile == nullptr) {
            return point_inside_mesh(mesh, point);
        }
        const auto start = detail::MedialProfileClock::now();
        const bool inside = point_inside_mesh(mesh, point);
        profile->point_containment_winding_seconds +=
            detail::medial_profile_elapsed_seconds(start);
        return inside;
    };
    const auto segment_intersects = [&](const Vec3& start_point,
                                        const Vec3& end_point) {
        if (profile == nullptr) {
            return segment_intersects_mesh_surface(
                mesh,
                start_point,
                end_point
            );
        }
        const auto start = detail::MedialProfileClock::now();
        const bool intersects = segment_intersects_mesh_surface(
            mesh,
            start_point,
            end_point
        );
        profile->segment_intersections_seconds +=
            detail::medial_profile_elapsed_seconds(start);
        return intersects;
    };

    Vec3 center{};
    for (const Vec3& vertex : polygon) {
        if (!point_is_inside(vertex)) {
            return false;
        }
        center = center + vertex;
    }
    center = center / static_cast<double>(polygon.size());
    if (!point_is_inside(center)) {
        return false;
    }

    // Check every true polygon edge plus every diagonal introduced by the
    // fan triangulation. Endpoint tests alone are insufficient in a concave
    // solid: an interior-to-interior chord can cross the exterior twice.
    for (std::size_t vertex = 0; vertex < polygon.size(); ++vertex) {
        if (segment_intersects(
                polygon[vertex],
                polygon[(vertex + 1) % polygon.size()])) {
            return false;
        }
    }
    for (std::size_t vertex = 2; vertex + 1 < polygon.size(); ++vertex) {
        if (segment_intersects(polygon[0], polygon[vertex])) {
            return false;
        }
    }
    for (std::size_t vertex = 1; vertex + 1 < polygon.size(); ++vertex) {
        const Vec3 triangle_center =
            (polygon[0] + polygon[vertex] + polygon[vertex + 1]) / 3.0;
        if (!point_is_inside(triangle_center)) {
            return false;
        }
    }
    return true;
}

inline BoundarySurface extract_boundary_surface(const Delaunay3& delaunay) {
    struct FaceRecord {
        std::array<int, 3> oriented_vertices{};
        int count{0};
    };

    std::map<std::array<int, 3>, FaceRecord> records;
    for (const Tetrahedron& tetrahedron : delaunay.tetrahedra()) {
        const std::array<std::array<int, 3>, 4> faces{{
            Delaunay3::face_from_index(tetrahedron, 0),
            Delaunay3::face_from_index(tetrahedron, 1),
            Delaunay3::face_from_index(tetrahedron, 2),
            Delaunay3::face_from_index(tetrahedron, 3)
        }};

        for (const auto& face : faces) {
            std::array<int, 3> key = face;
            std::sort(key.begin(), key.end());
            FaceRecord& record = records[key];
            if (record.count == 0) {
                record.oriented_vertices = face;
            }
            ++record.count;
        }
    }

    BoundarySurface surface;
    for (const auto& [key, record] : records) {
        (void)key;
        if (record.count == 1) {
            surface.faces.push_back({
                static_cast<std::size_t>(record.oriented_vertices[0]),
                static_cast<std::size_t>(record.oriented_vertices[1]),
                static_cast<std::size_t>(record.oriented_vertices[2])
            });
        }
    }
    return surface;
}

inline bool point_in_tetrahedron(const Vec3& point,
                                 const Tetrahedron& tetrahedron,
                                 const std::vector<Vec3>& points) {
    const Vec3& a = points[static_cast<std::size_t>(tetrahedron.vertices[0])];
    const Vec3& b = points[static_cast<std::size_t>(tetrahedron.vertices[1])];
    const Vec3& c = points[static_cast<std::size_t>(tetrahedron.vertices[2])];
    const Vec3& d = points[static_cast<std::size_t>(tetrahedron.vertices[3])];

    const auto same_side = [](double reference, double candidate) {
        return reference == 0.0 || candidate == 0.0 ||
               (reference > 0.0 && candidate > 0.0) ||
               (reference < 0.0 && candidate < 0.0);
    };

    return same_side(orient3d(a, b, c, d), orient3d(a, b, c, point)) &&
           same_side(orient3d(a, b, d, c), orient3d(a, b, d, point)) &&
           same_side(orient3d(a, c, d, b), orient3d(a, c, d, point)) &&
           same_side(orient3d(b, c, d, a), orient3d(b, c, d, point));
}

inline bool point_in_tetrahedral_volume(const Vec3& point, const Delaunay3& delaunay) {
    for (const Tetrahedron& tetrahedron : delaunay.tetrahedra()) {
        if (point_in_tetrahedron(point, tetrahedron, delaunay.points())) {
            return true;
        }
    }
    return false;
}

inline bool point_in_medial_domain(const Vec3& point,
                                   const Delaunay3& delaunay,
                                   const Mesh* surface_mesh) {
    return surface_mesh != nullptr
        ? point_inside_mesh(*surface_mesh, point)
        : point_in_tetrahedral_volume(point, delaunay);
}

inline MedialAxisApproximation build_medial_axis_approximation(
    const Delaunay3& delaunay,
    const Mesh* surface_mesh = nullptr) {
    MedialAxisApproximation result;
    std::vector<int> candidate_for_tetrahedron(delaunay.tetrahedron_count(), -1);

    for (std::size_t tetrahedron_index = 0;
         tetrahedron_index < delaunay.tetrahedron_count();
         ++tetrahedron_index) {
        Vec3 center{};
        if (!delaunay.circumcenter(tetrahedron_index, center) ||
            !point_in_medial_domain(center, delaunay, surface_mesh)) {
            continue;
        }

        const Tetrahedron& tetrahedron = delaunay.tetrahedra()[tetrahedron_index];
        const Vec3& surface_point =
            delaunay.points()[static_cast<std::size_t>(tetrahedron.vertices[0])];
        candidate_for_tetrahedron[tetrahedron_index] =
            static_cast<int>(result.centers.size());
        result.centers.push_back(center);
        result.radii.push_back(norm(center - surface_point));
        result.source_tetrahedra.push_back(tetrahedron_index);
    }

    for (std::size_t tetrahedron_index = 0;
         tetrahedron_index < delaunay.tetrahedron_count();
         ++tetrahedron_index) {
        const int candidate = candidate_for_tetrahedron[tetrahedron_index];
        if (candidate < 0) {
            continue;
        }

        for (int neighbor : delaunay.tetrahedra()[tetrahedron_index].neighbors) {
            if (neighbor <= static_cast<int>(tetrahedron_index) ||
                neighbor >= static_cast<int>(candidate_for_tetrahedron.size())) {
                continue;
            }
            const int neighbor_candidate =
                candidate_for_tetrahedron[static_cast<std::size_t>(neighbor)];
            if (neighbor_candidate < 0) {
                continue;
            }
            if (surface_mesh != nullptr &&
                !segment_inside_mesh(
                    *surface_mesh,
                    result.centers[static_cast<std::size_t>(candidate)],
                    result.centers[
                        static_cast<std::size_t>(neighbor_candidate)])) {
                ++result.rejected_exterior_edge_count;
                continue;
            }
            result.edges.push_back({
                static_cast<std::size_t>(candidate),
                static_cast<std::size_t>(neighbor_candidate)
            });
        }
    }

    return result;
}

inline std::vector<VoronoiFaceCandidate> build_interior_voronoi_faces(
    const Delaunay3& delaunay,
    const Mesh* surface_mesh = nullptr,
    ValidatedMedialComplexProfile* profile = nullptr) {
    using Edge = std::array<int, 2>;
    const auto candidate_start = detail::MedialProfileClock::now();
    const double containment_before =
        profile != nullptr
            ? profile->point_containment_winding_seconds
            : 0.0;

    const auto normalize_edge = [](int first, int second) {
        return first < second ? Edge{first, second} : Edge{second, first};
    };

    std::set<Edge> boundary_edges;
    std::vector<bool> is_surface_vertex(delaunay.point_count(), false);
    const bool mesh_vertices_match_delaunay =
        surface_mesh != nullptr &&
        surface_mesh->vertices.size() == delaunay.point_count() &&
        std::equal(
            surface_mesh->vertices.begin(),
            surface_mesh->vertices.end(),
            delaunay.points().begin()
        );
    if (mesh_vertices_match_delaunay) {
        for (const Triangle& face : surface_mesh->faces) {
            for (int vertex : face.vertices) {
                is_surface_vertex[static_cast<std::size_t>(vertex)] = true;
            }
            boundary_edges.insert(normalize_edge(face.vertices[0], face.vertices[1]));
            boundary_edges.insert(normalize_edge(face.vertices[1], face.vertices[2]));
            boundary_edges.insert(normalize_edge(face.vertices[2], face.vertices[0]));
        }
    } else {
        const BoundarySurface boundary = extract_boundary_surface(delaunay);
        for (const auto& face : boundary.faces) {
            for (std::size_t vertex : face) {
                is_surface_vertex[vertex] = true;
            }
            boundary_edges.insert(normalize_edge(
                static_cast<int>(face[0]),
                static_cast<int>(face[1])
            ));
            boundary_edges.insert(normalize_edge(
                static_cast<int>(face[1]),
                static_cast<int>(face[2])
            ));
            boundary_edges.insert(normalize_edge(
                static_cast<int>(face[2]),
                static_cast<int>(face[0])
            ));
        }
        if (surface_mesh != nullptr) {
            // A resampled Delaunay point set no longer shares indices with the
            // original mesh, but every generated point is a surface contact.
            std::fill(is_surface_vertex.begin(), is_surface_vertex.end(), true);
        }
    }

    std::map<Edge, std::vector<std::size_t>> incident_tetrahedra;
    for (std::size_t tetrahedron_index = 0;
         tetrahedron_index < delaunay.tetrahedron_count();
         ++tetrahedron_index) {
        const auto& vertices = delaunay.tetrahedra()[tetrahedron_index].vertices;
        for (int first = 0; first < 4; ++first) {
            for (int second = first + 1; second < 4; ++second) {
                incident_tetrahedra[normalize_edge(
                    vertices[static_cast<std::size_t>(first)],
                    vertices[static_cast<std::size_t>(second)]
                )].push_back(tetrahedron_index);
            }
        }
    }

    std::vector<Vec3> circumcenters(delaunay.tetrahedron_count());
    std::vector<bool> center_is_interior(delaunay.tetrahedron_count(), false);
    for (std::size_t tetrahedron_index = 0;
         tetrahedron_index < delaunay.tetrahedron_count();
         ++tetrahedron_index) {
        Vec3 center{};
        if (!delaunay.circumcenter(tetrahedron_index, center)) {
            continue;
        }
        bool center_is_inside = false;
        if (profile != nullptr) {
            const auto containment_start =
                detail::MedialProfileClock::now();
            center_is_inside =
                point_in_medial_domain(center, delaunay, surface_mesh);
            profile->point_containment_winding_seconds +=
                detail::medial_profile_elapsed_seconds(
                    containment_start
                );
        } else {
            center_is_inside =
                point_in_medial_domain(center, delaunay, surface_mesh);
        }
        if (center_is_inside) {
            circumcenters[tetrahedron_index] = center;
            center_is_interior[tetrahedron_index] = true;
        }
    }

    constexpr double radians_to_degrees = 57.2957795130823208768;
    std::vector<VoronoiFaceCandidate> result;
    for (const auto& [edge, incident] : incident_tetrahedra) {
        if (boundary_edges.count(edge) != 0 ||
            !is_surface_vertex[static_cast<std::size_t>(edge[0])] ||
            !is_surface_vertex[static_cast<std::size_t>(edge[1])] ||
            incident.size() < 3) {
            continue;
        }

        struct PolygonVertex {
            Vec3 position{};
            std::size_t source_tetrahedron{0};
        };
        std::vector<PolygonVertex> polygon;
        for (std::size_t tetrahedron_index : incident) {
            if (!center_is_interior[tetrahedron_index]) {
                continue;
            }
            const Vec3 center = circumcenters[tetrahedron_index];
            const bool duplicate = std::any_of(
                polygon.begin(),
                polygon.end(),
                [&](const PolygonVertex& existing) {
                    return squared_norm(existing.position - center) < 1e-20;
                }
            );
            if (!duplicate) {
                polygon.push_back({center, tetrahedron_index});
            }
        }
        if (polygon.size() < 3) {
            continue;
        }

        const Vec3& first_site =
            delaunay.points()[static_cast<std::size_t>(edge[0])];
        const Vec3& second_site =
            delaunay.points()[static_cast<std::size_t>(edge[1])];
        const Vec3 edge_axis = normalized(second_site - first_site);
        const Vec3 helper = std::fabs(edge_axis.x) < 0.9
            ? Vec3(1.0, 0.0, 0.0)
            : Vec3(0.0, 1.0, 0.0);
        const Vec3 basis_u = normalized(cross(edge_axis, helper));
        const Vec3 basis_v = cross(edge_axis, basis_u);
        const Vec3 edge_midpoint = (first_site + second_site) / 2.0;

        std::sort(
            polygon.begin(),
            polygon.end(),
            [&](const PolygonVertex& a, const PolygonVertex& b) {
                const Vec3 offset_a = a.position - edge_midpoint;
                const Vec3 offset_b = b.position - edge_midpoint;
                const double angle_a =
                    std::atan2(dot(offset_a, basis_v), dot(offset_a, basis_u));
                const double angle_b =
                    std::atan2(dot(offset_b, basis_v), dot(offset_b, basis_u));
                return angle_a < angle_b;
            }
        );

        Vec3 centroid{};
        for (const PolygonVertex& vertex : polygon) {
            centroid = centroid + vertex.position;
        }
        centroid = centroid / static_cast<double>(polygon.size());

        const Vec3 first_direction = normalized(first_site - centroid);
        const Vec3 second_direction = normalized(second_site - centroid);
        const double cosine = std::clamp(
            dot(first_direction, second_direction),
            -1.0,
            1.0
        );

        VoronoiFaceCandidate face;
        face.contact_sites = edge;
        face.vertices.reserve(polygon.size());
        face.source_tetrahedra.reserve(polygon.size());
        for (const PolygonVertex& vertex : polygon) {
            face.vertices.push_back(vertex.position);
            face.source_tetrahedra.push_back(vertex.source_tetrahedron);
        }
        face.contact_angle_degrees = std::acos(cosine) * radians_to_degrees;
        result.push_back(std::move(face));
    }

    if (profile != nullptr) {
        const double elapsed =
            detail::medial_profile_elapsed_seconds(candidate_start);
        const double containment =
            profile->point_containment_winding_seconds -
            containment_before;
        profile->candidate_generation_seconds +=
            std::max(0.0, elapsed - containment);
    }
    return result;
}

inline MedialSheetApproximation build_medial_sheet_approximation(
    const Delaunay3& delaunay,
    const Mesh* surface_mesh = nullptr,
    double minimum_contact_angle_degrees = 110.0) {
    MedialSheetApproximation result;
    const auto faces = build_interior_voronoi_faces(delaunay, surface_mesh);

    for (const VoronoiFaceCandidate& face : faces) {
        if (face.contact_angle_degrees < minimum_contact_angle_degrees ||
            face.vertices.size() < 3 ||
            (surface_mesh != nullptr &&
             !polygon_fan_inside_mesh(*surface_mesh, face.vertices))) {
            continue;
        }

        const std::size_t vertex_offset = result.vertices.size();
        result.vertices.insert(
            result.vertices.end(),
            face.vertices.begin(),
            face.vertices.end()
        );

        std::size_t triangle_count_before = result.triangles.size();
        for (std::size_t i = 1; i + 1 < face.vertices.size(); ++i) {
            const Vec3& a = face.vertices[0];
            const Vec3& b = face.vertices[i];
            const Vec3& c = face.vertices[i + 1];
            const Vec3 first_edge = b - a;
            const Vec3 second_edge = c - a;
            if (squared_norm(cross(first_edge, second_edge)) < 1e-20) {
                continue;
            }
            result.triangles.push_back({
                vertex_offset,
                vertex_offset + i,
                vertex_offset + i + 1
            });
            result.triangle_contact_angles.push_back(face.contact_angle_degrees);
        }

        if (result.triangles.size() > triangle_count_before) {
            ++result.polygon_count;
        } else {
            result.vertices.resize(vertex_offset);
        }
    }

    return result;
}

}  // namespace medial_axis_3d
