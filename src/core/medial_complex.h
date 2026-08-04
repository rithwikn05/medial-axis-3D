#pragma once

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <vector>

#include "local_feature_size.h"
#include "medial_axis_approx.h"
#include "medial_candidate.h"

namespace medial_axis_3d {

struct MedialComplexOptions {
    // A weak threshold keeps tapering parts of a medial sheet.  A separate
    // strong threshold supplies high-confidence seeds.  Using the strong
    // threshold as a per-polygon deletion test fragments continuous sheets.
    double minimum_contact_angle_degrees{35.0};
    double strong_contact_angle_degrees{90.0};
    double minimum_triangle_confidence{0.30};
    double propagated_support_decay{0.70};
    double termination_angle_margin_degrees{15.0};
    double maximum_gap_radius_jump{0.85};
    double minimum_gap_normal_alignment{0.35};
    std::size_t maximum_gap_triangles{512};
    double maximum_gap_area_fraction{0.20};
    // Zero selects an automatic, resolution-scaled ring count.
    std::size_t support_ring_count{0};
};

struct MedialBoundaryLoop {
    std::size_t component{0};
    std::vector<std::size_t> vertices;
    double length{0.0};
    bool closed{false};
    bool touches_seam_or_junction{false};
    bool termination_supported{false};
    bool allowed{false};
};

struct MedialComplex {
    std::vector<Vec3> vertices;
    // Exact distance to the triangle surface when a mesh is available;
    // otherwise the generating Delaunay circumsphere radius.
    std::vector<double> vertex_radii;
    std::vector<double> vertex_surface_resolutions;
    std::vector<double> vertex_local_feature_sizes;
    std::vector<double> vertex_sampling_densities;
    // A welded Voronoi vertex may represent multiple co-spherical
    // tetrahedra, so retain every generating tetrahedron.
    std::vector<std::vector<std::size_t>> vertex_source_tetrahedra;
    std::vector<std::array<std::size_t, 3>> triangles;
    // Pole support and confidence are weights. They do not independently
    // remove faces from an otherwise coherent sheet.
    std::vector<double> triangle_support_weights;
    std::vector<double> triangle_confidences;
    std::vector<double> triangle_contact_angles;
    std::vector<double> triangle_radius_jumps;
    std::vector<double> triangle_radius_jump_thresholds;
    std::vector<std::size_t> triangle_components;

    // Incidence-one edges are exposed separately from true non-manifold seams.
    std::vector<std::array<std::size_t, 2>> boundary_edges;
    std::vector<std::array<std::size_t, 2>> seam_edges;
    std::vector<std::size_t> junction_vertices;
    std::vector<std::array<std::size_t, 2>> termination_edges;
    std::vector<std::array<std::size_t, 2>> artificial_boundary_edges;
    std::vector<MedialBoundaryLoop> boundary_loops;
    // Preserve seam/junction provenance when neighboring whole sheets are
    // removed; a former singular curve remains a valid sheet boundary.
    std::vector<std::size_t> permitted_boundary_vertices;

    // Rejected face centroids are useful as a disabled diagnostic layer.
    std::vector<Vec3> rejected_face_centers;
    std::vector<double> rejected_face_confidences;
    std::vector<std::array<std::size_t, 3>>
        topology_restored_candidate_triangles;

    std::size_t source_polygon_count{0};
    std::size_t accepted_polygon_count{0};
    std::size_t topology_restored_candidate_patch_count{0};
    std::size_t component_count{0};
};

namespace detail {

using ComplexEdge = std::array<std::size_t, 2>;

inline ComplexEdge normalized_complex_edge(std::size_t first,
                                           std::size_t second) {
    return first < second
        ? ComplexEdge{first, second}
        : ComplexEdge{second, first};
}

inline double relative_radius_jump(double first,
                                   double second,
                                   double minimum_scale = 1e-12) {
    const double scale = std::max(
        {std::fabs(first), std::fabs(second), minimum_scale}
    );
    return std::fabs(first - second) / scale;
}

inline double triangle_max_radius_jump(
    const MedialComplex& complex,
    const std::array<std::size_t, 3>& triangle) {
    if (complex.vertex_radii.size() != complex.vertices.size()) {
        return 0.0;
    }
    return std::max({
        relative_radius_jump(
            complex.vertex_radii[triangle[0]],
            complex.vertex_radii[triangle[1]]
        ),
        relative_radius_jump(
            complex.vertex_radii[triangle[1]],
            complex.vertex_radii[triangle[2]]
        ),
        relative_radius_jump(
            complex.vertex_radii[triangle[2]],
            complex.vertex_radii[triangle[0]]
        )
    });
}

inline Vec3 polygon_center(const VoronoiFaceCandidate& face) {
    Vec3 center{};
    for (const Vec3& vertex : face.vertices) {
        center = center + vertex;
    }
    return face.vertices.empty()
        ? center
        : center / static_cast<double>(face.vertices.size());
}

inline std::size_t shared_source_tetrahedron_count(
    const VoronoiFaceCandidate& first,
    const VoronoiFaceCandidate& second) {
    std::size_t count = 0;
    for (std::size_t source : first.source_tetrahedra) {
        if (std::find(
                second.source_tetrahedra.begin(),
                second.source_tetrahedra.end(),
                source) != second.source_tetrahedra.end()) {
            ++count;
        }
    }
    return count;
}

inline bool shares_contact_site(const VoronoiFaceCandidate& first,
                                const VoronoiFaceCandidate& second) {
    return first.contact_sites[0] == second.contact_sites[0] ||
           first.contact_sites[0] == second.contact_sites[1] ||
           first.contact_sites[1] == second.contact_sites[0] ||
           first.contact_sites[1] == second.contact_sites[1];
}

inline int triangle_edge_direction(
    const std::array<std::size_t, 3>& triangle,
    const ComplexEdge& edge) {
    for (std::size_t i = 0; i < 3; ++i) {
        const std::size_t first = triangle[i];
        const std::size_t second = triangle[(i + 1) % 3];
        if (first == edge[0] && second == edge[1]) {
            return 1;
        }
        if (first == edge[1] && second == edge[0]) {
            return -1;
        }
    }
    return 0;
}

inline std::map<ComplexEdge, std::vector<std::size_t>>
triangle_edge_incidence(
    const std::vector<std::array<std::size_t, 3>>& triangles) {
    std::map<ComplexEdge, std::vector<std::size_t>> result;
    for (std::size_t triangle_index = 0;
         triangle_index < triangles.size();
         ++triangle_index) {
        const auto& triangle = triangles[triangle_index];
        for (std::size_t edge = 0; edge < 3; ++edge) {
            result[normalized_complex_edge(
                triangle[edge],
                triangle[(edge + 1) % 3]
            )].push_back(triangle_index);
        }
    }
    return result;
}

inline Vec3 medial_triangle_normal(
    const MedialComplex& complex,
    const std::array<std::size_t, 3>& triangle) {
    const Vec3 normal = cross(
        complex.vertices[triangle[1]] - complex.vertices[triangle[0]],
        complex.vertices[triangle[2]] - complex.vertices[triangle[0]]
    );
    const double length = norm(normal);
    return length > 1e-15 ? normal / length : Vec3{};
}

inline void orient_and_label_sheet_components(MedialComplex& complex) {
    const auto incidence = triangle_edge_incidence(complex.triangles);
    std::vector<std::vector<std::pair<std::size_t, bool>>> adjacency(
        complex.triangles.size()
    );
    for (const auto& [edge, incident] : incidence) {
        if (incident.size() != 2) {
            continue;
        }
        const std::size_t first = incident[0];
        const std::size_t second = incident[1];
        const bool same_direction =
            triangle_edge_direction(complex.triangles[first], edge) ==
            triangle_edge_direction(complex.triangles[second], edge);
        adjacency[first].push_back({second, same_direction});
        adjacency[second].push_back({first, same_direction});
    }

    const std::size_t unassigned = std::numeric_limits<std::size_t>::max();
    complex.triangle_components.assign(complex.triangles.size(), unassigned);
    std::vector<bool> flipped(complex.triangles.size(), false);
    std::size_t component = 0;
    for (std::size_t seed = 0; seed < complex.triangles.size(); ++seed) {
        if (complex.triangle_components[seed] != unassigned) {
            continue;
        }
        std::queue<std::size_t> pending;
        complex.triangle_components[seed] = component;
        pending.push(seed);
        while (!pending.empty()) {
            const std::size_t current = pending.front();
            pending.pop();
            for (const auto& [neighbor, same_direction] : adjacency[current]) {
                const bool neighbor_flipped =
                    flipped[current] != same_direction;
                if (complex.triangle_components[neighbor] == unassigned) {
                    complex.triangle_components[neighbor] = component;
                    flipped[neighbor] = neighbor_flipped;
                    pending.push(neighbor);
                }
            }
        }
        ++component;
    }

    for (std::size_t triangle = 0;
         triangle < complex.triangles.size();
         ++triangle) {
        if (flipped[triangle]) {
            std::swap(
                complex.triangles[triangle][1],
                complex.triangles[triangle][2]
            );
        }
    }
    complex.component_count = component;
}

inline void classify_complex_topology(
    MedialComplex& complex,
    double termination_contact_angle_degrees = 105.0) {
    complex.boundary_edges.clear();
    complex.seam_edges.clear();
    complex.junction_vertices.clear();
    complex.termination_edges.clear();
    complex.artificial_boundary_edges.clear();
    complex.boundary_loops.clear();

    const auto incidence = triangle_edge_incidence(complex.triangles);
    std::map<std::size_t, std::size_t> seam_degree;
    for (const auto& [edge, incident] : incidence) {
        if (incident.size() == 1) {
            complex.boundary_edges.push_back(edge);
        } else if (incident.size() > 2) {
            complex.seam_edges.push_back(edge);
            ++seam_degree[edge[0]];
            ++seam_degree[edge[1]];
        }
    }
    for (const auto& [vertex, degree] : seam_degree) {
        if (degree != 2) {
            complex.junction_vertices.push_back(vertex);
        }
    }

    if (complex.boundary_edges.empty() ||
        complex.triangle_components.size() != complex.triangles.size()) {
        return;
    }

    std::set<std::size_t> seam_vertices;
    for (const ComplexEdge& edge : complex.seam_edges) {
        seam_vertices.insert(edge[0]);
        seam_vertices.insert(edge[1]);
    }
    const std::set<std::size_t> junction_set(
        complex.junction_vertices.begin(),
        complex.junction_vertices.end()
    );
    seam_vertices.insert(
        complex.permitted_boundary_vertices.begin(),
        complex.permitted_boundary_vertices.end()
    );
    seam_vertices.insert(junction_set.begin(), junction_set.end());
    complex.permitted_boundary_vertices.assign(
        seam_vertices.begin(),
        seam_vertices.end()
    );

    std::map<std::size_t, std::vector<ComplexEdge>>
        component_boundary_edges;
    std::map<ComplexEdge, bool> termination_supported;
    for (const ComplexEdge& edge : complex.boundary_edges) {
        const auto found = incidence.find(edge);
        if (found == incidence.end() || found->second.size() != 1) {
            continue;
        }
        const std::size_t triangle = found->second.front();
        if (triangle >= complex.triangle_components.size()) {
            continue;
        }
        component_boundary_edges[
            complex.triangle_components[triangle]
        ].push_back(edge);

        bool is_termination = false;
        if (triangle < complex.triangle_contact_angles.size()) {
            is_termination =
                complex.triangle_contact_angles[triangle] <=
                    termination_contact_angle_degrees;
        }
        if (!is_termination &&
            complex.vertex_radii.size() == complex.vertices.size() &&
            complex.vertex_local_feature_sizes.size() ==
                complex.vertices.size()) {
            double ratio = 0.0;
            std::size_t count = 0;
            for (std::size_t vertex : edge) {
                const double lfs =
                    complex.vertex_local_feature_sizes[vertex];
                if (lfs > 1e-12) {
                    ratio += complex.vertex_radii[vertex] / lfs;
                    ++count;
                }
            }
            is_termination =
                count > 0 && ratio / static_cast<double>(count) <= 0.25;
        }
        termination_supported[edge] = is_termination;
        if (is_termination) {
            complex.termination_edges.push_back(edge);
        }
    }

    for (const auto& [component, edges] : component_boundary_edges) {
        std::map<std::size_t, std::vector<std::size_t>> edge_indices_at_vertex;
        for (std::size_t edge_index = 0;
             edge_index < edges.size();
             ++edge_index) {
            edge_indices_at_vertex[edges[edge_index][0]].push_back(edge_index);
            edge_indices_at_vertex[edges[edge_index][1]].push_back(edge_index);
        }

        std::vector<bool> visited(edges.size(), false);
        for (std::size_t seed = 0; seed < edges.size(); ++seed) {
            if (visited[seed]) {
                continue;
            }
            std::vector<std::size_t> group;
            std::queue<std::size_t> pending;
            visited[seed] = true;
            pending.push(seed);
            while (!pending.empty()) {
                const std::size_t edge_index = pending.front();
                pending.pop();
                group.push_back(edge_index);
                for (std::size_t vertex : edges[edge_index]) {
                    for (std::size_t neighbor :
                         edge_indices_at_vertex[vertex]) {
                        if (!visited[neighbor]) {
                            visited[neighbor] = true;
                            pending.push(neighbor);
                        }
                    }
                }
            }

            MedialBoundaryLoop loop;
            loop.component = component;
            std::map<std::size_t, std::size_t> degrees;
            std::set<std::size_t> vertices;
            std::size_t termination_edge_count = 0;
            for (std::size_t edge_index : group) {
                const ComplexEdge& edge = edges[edge_index];
                ++degrees[edge[0]];
                ++degrees[edge[1]];
                vertices.insert(edge[0]);
                vertices.insert(edge[1]);
                loop.length += norm(
                    complex.vertices[edge[1]] -
                    complex.vertices[edge[0]]
                );
                if (termination_supported[edge]) {
                    ++termination_edge_count;
                }
            }
            loop.closed =
                group.size() >= 3 &&
                std::all_of(
                    degrees.begin(),
                    degrees.end(),
                    [](const auto& entry) {
                        return entry.second == 2;
                    }
                );
            loop.vertices.assign(vertices.begin(), vertices.end());
            loop.touches_seam_or_junction = std::any_of(
                vertices.begin(),
                vertices.end(),
                [&](std::size_t vertex) {
                    return seam_vertices.count(vertex) != 0 ||
                           junction_set.count(vertex) != 0;
                }
            );
            loop.termination_supported =
                termination_edge_count * 2 >= group.size();
            loop.allowed =
                loop.touches_seam_or_junction ||
                loop.termination_supported;
            if (!loop.allowed) {
                for (std::size_t edge_index : group) {
                    complex.artificial_boundary_edges.push_back(
                        edges[edge_index]
                    );
                }
            }
            complex.boundary_loops.push_back(std::move(loop));
        }
    }
}

struct CandidateGapRepair {
    std::vector<bool> repaired;
    std::size_t repaired_patch_count{0};
};

struct CandidatePolygonGapRepair {
    std::vector<bool> repaired;
    std::size_t repaired_patch_count{0};
};

struct CandidatePolygonStratumCompletion {
    std::vector<bool> repaired;
    std::size_t completed_stratum_count{0};
};

inline CandidatePolygonStratumCompletion complete_seeded_polygon_strata(
    const std::vector<VoronoiFaceCandidate>& faces,
    const std::vector<bool>& eligible,
    const std::vector<bool>& inside,
    std::vector<bool>& keep) {
    CandidatePolygonStratumCompletion result;
    result.repaired.assign(faces.size(), false);
    if (eligible.size() != faces.size() ||
        inside.size() != faces.size() ||
        keep.size() != faces.size()) {
        return result;
    }

    std::map<ComplexEdge, std::vector<std::size_t>> edge_incidence;
    for (std::size_t face_index = 0;
         face_index < faces.size();
         ++face_index) {
        if (!eligible[face_index] || !inside[face_index]) {
            continue;
        }
        const auto& sources = faces[face_index].source_tetrahedra;
        for (std::size_t vertex = 0;
             vertex < sources.size();
             ++vertex) {
            edge_incidence[normalized_complex_edge(
                sources[vertex],
                sources[(vertex + 1) % sources.size()]
            )].push_back(face_index);
        }
    }

    // Incidence-two edges continue one ordinary 2D sheet. Incidence-three
    // or greater edges are medial seams where separate strata meet; do not
    // flood across them.
    std::vector<std::vector<std::size_t>> adjacency(faces.size());
    for (const auto& [edge, incident] : edge_incidence) {
        (void)edge;
        if (incident.size() != 2 || incident[0] == incident[1]) {
            continue;
        }
        adjacency[incident[0]].push_back(incident[1]);
        adjacency[incident[1]].push_back(incident[0]);
    }

    std::vector<bool> visited(faces.size(), false);
    for (std::size_t seed = 0; seed < faces.size(); ++seed) {
        if (visited[seed] || !eligible[seed] || !inside[seed]) {
            continue;
        }
        std::vector<std::size_t> stratum;
        std::queue<std::size_t> pending;
        bool has_seed = false;
        visited[seed] = true;
        pending.push(seed);
        while (!pending.empty()) {
            const std::size_t face_index = pending.front();
            pending.pop();
            stratum.push_back(face_index);
            has_seed = has_seed || keep[face_index];
            for (std::size_t neighbor : adjacency[face_index]) {
                if (!visited[neighbor]) {
                    visited[neighbor] = true;
                    pending.push(neighbor);
                }
            }
        }
        if (!has_seed) {
            continue;
        }
        bool repaired_any = false;
        for (std::size_t face_index : stratum) {
            if (!keep[face_index]) {
                keep[face_index] = true;
                result.repaired[face_index] = true;
                repaired_any = true;
            }
        }
        if (repaired_any) {
            ++result.completed_stratum_count;
        }
    }
    return result;
}

inline CandidatePolygonGapRepair restore_enclosed_polygon_gaps(
    const std::vector<VoronoiFaceCandidate>& faces,
    const std::vector<bool>& eligible,
    const std::vector<bool>& inside,
    std::vector<bool>& keep) {
    CandidatePolygonGapRepair result;
    result.repaired.assign(faces.size(), false);
    if (eligible.size() != faces.size() ||
        inside.size() != faces.size() ||
        keep.size() != faces.size()) {
        return result;
    }

    // A Voronoi polygon edge joins two consecutive generating
    // tetrahedra. Polygons sharing that pair meet along the same geometric
    // Voronoi edge. Work with these atomic 2-cells so a repair can never
    // restore only part of an arbitrary fan triangulation.
    std::map<ComplexEdge, std::vector<std::size_t>> edge_incidence;
    std::vector<std::vector<ComplexEdge>> polygon_edges(faces.size());
    for (std::size_t face_index = 0;
         face_index < faces.size();
         ++face_index) {
        if (!eligible[face_index] || !inside[face_index]) {
            continue;
        }
        const auto& sources = faces[face_index].source_tetrahedra;
        if (sources.size() < 3) {
            continue;
        }
        for (std::size_t vertex = 0;
             vertex < sources.size();
             ++vertex) {
            const ComplexEdge edge = normalized_complex_edge(
                sources[vertex],
                sources[(vertex + 1) % sources.size()]
            );
            polygon_edges[face_index].push_back(edge);
            edge_incidence[edge].push_back(face_index);
        }
    }

    std::vector<std::vector<std::size_t>> rejected_adjacency(faces.size());
    for (const auto& [edge, incident] : edge_incidence) {
        (void)edge;
        for (std::size_t first = 0; first < incident.size(); ++first) {
            for (std::size_t second = first + 1;
                 second < incident.size();
                 ++second) {
                const std::size_t a = incident[first];
                const std::size_t b = incident[second];
                if (!keep[a] && !keep[b]) {
                    rejected_adjacency[a].push_back(b);
                    rejected_adjacency[b].push_back(a);
                }
            }
        }
    }

    std::vector<bool> visited(faces.size(), false);
    for (std::size_t seed = 0; seed < faces.size(); ++seed) {
        if (visited[seed] || keep[seed] ||
            !eligible[seed] || !inside[seed]) {
            continue;
        }

        std::vector<std::size_t> patch;
        std::queue<std::size_t> pending;
        visited[seed] = true;
        pending.push(seed);
        while (!pending.empty()) {
            const std::size_t face_index = pending.front();
            pending.pop();
            patch.push_back(face_index);
            for (std::size_t neighbor :
                 rejected_adjacency[face_index]) {
                if (!visited[neighbor]) {
                    visited[neighbor] = true;
                    pending.push(neighbor);
                }
            }
        }

        const std::set<std::size_t> patch_set(
            patch.begin(),
            patch.end()
        );
        bool enclosed = !patch.empty();
        std::size_t retained_boundary_edge_count = 0;
        for (std::size_t face_index : patch) {
            for (const ComplexEdge& edge : polygon_edges[face_index]) {
                const auto found = edge_incidence.find(edge);
                if (found == edge_incidence.end()) {
                    enclosed = false;
                    continue;
                }
                bool has_patch_neighbor = false;
                bool has_retained_neighbor = false;
                for (std::size_t incident_face : found->second) {
                    if (incident_face == face_index) {
                        continue;
                    }
                    has_patch_neighbor =
                        has_patch_neighbor ||
                        patch_set.count(incident_face) != 0;
                    has_retained_neighbor =
                        has_retained_neighbor || keep[incident_face];
                }
                if (!has_patch_neighbor) {
                    if (has_retained_neighbor) {
                        ++retained_boundary_edge_count;
                    } else {
                        // The rejected region reaches the boundary of the
                        // candidate complex, so this is a legitimate sheet
                        // termination rather than an enclosed hole.
                        enclosed = false;
                    }
                }
            }
        }
        if (!enclosed || retained_boundary_edge_count < 3) {
            continue;
        }
        for (std::size_t face_index : patch) {
            keep[face_index] = true;
            result.repaired[face_index] = true;
        }
        ++result.repaired_patch_count;
    }
    return result;
}

inline CandidateGapRepair restore_enclosed_candidate_gaps(
    const MedialComplex& reference,
    std::vector<bool>& keep,
    const MedialComplexOptions& options) {
    CandidateGapRepair result;
    result.repaired.assign(reference.triangles.size(), false);
    if (reference.triangles.empty() ||
        keep.size() != reference.triangles.size() ||
        reference.triangle_components.size() !=
            reference.triangles.size()) {
        return result;
    }

    const auto incidence = triangle_edge_incidence(reference.triangles);
    std::vector<std::vector<std::size_t>> rejected_adjacency(
        reference.triangles.size()
    );
    for (const auto& [edge, incident] : incidence) {
        if (incident.size() != 2) {
            continue;
        }
        const std::size_t first = incident[0];
        const std::size_t second = incident[1];
        if (!keep[first] &&
            !keep[second] &&
            reference.triangle_components[first] ==
                reference.triangle_components[second]) {
            rejected_adjacency[first].push_back(second);
            rejected_adjacency[second].push_back(first);
        }
    }

    const auto triangle_area = [&](std::size_t triangle) {
        const auto& face = reference.triangles[triangle];
        return 0.5 * norm(cross(
            reference.vertices[face[1]] - reference.vertices[face[0]],
            reference.vertices[face[2]] - reference.vertices[face[0]]
        ));
    };
    std::vector<double> component_areas(reference.component_count, 0.0);
    for (std::size_t triangle = 0;
         triangle < reference.triangles.size();
         ++triangle) {
        const std::size_t component =
            reference.triangle_components[triangle];
        if (component < component_areas.size()) {
            component_areas[component] += triangle_area(triangle);
        }
    }

    std::vector<bool> visited(reference.triangles.size(), false);
    for (std::size_t seed = 0;
         seed < reference.triangles.size();
         ++seed) {
        if (keep[seed] || visited[seed]) {
            continue;
        }

        std::vector<std::size_t> patch;
        std::queue<std::size_t> pending;
        visited[seed] = true;
        pending.push(seed);
        while (!pending.empty()) {
            const std::size_t triangle = pending.front();
            pending.pop();
            patch.push_back(triangle);
            for (std::size_t neighbor :
                 rejected_adjacency[triangle]) {
                if (!visited[neighbor]) {
                    visited[neighbor] = true;
                    pending.push(neighbor);
                }
            }
        }

        if (patch.empty() ||
            patch.size() > options.maximum_gap_triangles) {
            continue;
        }
        const std::size_t component =
            reference.triangle_components[seed];
        if (component >= component_areas.size() ||
            component_areas[component] <= 0.0) {
            continue;
        }

        std::set<std::size_t> patch_set(patch.begin(), patch.end());
        std::set<std::size_t> boundary_neighbors;
        bool enclosed = true;
        bool touches_seam = false;
        double patch_area = 0.0;
        double maximum_radius_jump = 0.0;
        double minimum_normal_alignment = 1.0;
        for (std::size_t triangle : patch) {
            patch_area += triangle_area(triangle);
            maximum_radius_jump = std::max(
                maximum_radius_jump,
                triangle < reference.triangle_radius_jumps.size()
                    ? reference.triangle_radius_jumps[triangle]
                    : triangle_max_radius_jump(
                        reference,
                        reference.triangles[triangle]
                    )
            );
            const auto& face = reference.triangles[triangle];
            for (std::size_t edge_index = 0;
                 edge_index < 3;
                 ++edge_index) {
                const ComplexEdge edge = normalized_complex_edge(
                    face[edge_index],
                    face[(edge_index + 1) % 3]
                );
                const auto found = incidence.find(edge);
                if (found == incidence.end()) {
                    enclosed = false;
                    continue;
                }
                std::size_t patch_incidence = 0;
                for (std::size_t incident_triangle : found->second) {
                    if (patch_set.count(incident_triangle) != 0) {
                        ++patch_incidence;
                    }
                }
                if (patch_incidence != 1) {
                    continue;
                }
                if (found->second.size() == 1) {
                    enclosed = false;
                    continue;
                }
                if (found->second.size() > 2) {
                    touches_seam = true;
                    enclosed = false;
                }
                bool has_kept_neighbor = false;
                for (std::size_t neighbor : found->second) {
                    if (patch_set.count(neighbor) != 0 ||
                        !keep[neighbor] ||
                        reference.triangle_components[neighbor] !=
                            component) {
                        continue;
                    }
                    has_kept_neighbor = true;
                    boundary_neighbors.insert(neighbor);
                    minimum_normal_alignment = std::min(
                        minimum_normal_alignment,
                        std::fabs(dot(
                            medial_triangle_normal(
                                reference,
                                reference.triangles[triangle]
                            ),
                            medial_triangle_normal(
                                reference,
                                reference.triangles[neighbor]
                            )
                        ))
                    );
                }
                if (!has_kept_neighbor) {
                    enclosed = false;
                }
            }
        }

        const double area_fraction =
            patch_area / component_areas[component];
        if (!enclosed ||
            touches_seam ||
            boundary_neighbors.size() < 3 ||
            area_fraction > options.maximum_gap_area_fraction ||
            maximum_radius_jump > options.maximum_gap_radius_jump ||
            minimum_normal_alignment <
                options.minimum_gap_normal_alignment) {
            continue;
        }

        for (std::size_t triangle : patch) {
            keep[triangle] = true;
            result.repaired[triangle] = true;
        }
        ++result.repaired_patch_count;
    }
    return result;
}

}  // namespace detail

inline std::size_t resolved_medial_support_ring_count(
    std::size_t surface_sample_count,
    const MedialComplexOptions& options) {
    if (options.support_ring_count > 0) {
        return options.support_ring_count;
    }
    // Samples lie on a 2D surface, so their linear spacing scales with the
    // inverse square root of sample count. Grow proportionally more graph
    // rings to cover approximately the same physical neighborhood.
    const double relative_resolution = std::sqrt(
        static_cast<double>(std::max<std::size_t>(
            surface_sample_count,
            1
        )) / 40.0
    );
    return std::clamp<std::size_t>(
        static_cast<std::size_t>(std::ceil(relative_resolution)),
        1,
        4
    );
}

inline MedialComplex build_validated_medial_complex(
    const Delaunay3& delaunay,
    const Mesh* surface_mesh,
    const PoleSelectionResult& pole_selection,
    const MedialComplexOptions& options = {},
    const SurfaceFeatureField* surface_feature_field = nullptr,
    ValidatedMedialComplexProfile* profile = nullptr) {
    const auto total_start = detail::MedialProfileClock::now();
    const ValidatedMedialComplexProfile profile_before =
        profile != nullptr
            ? *profile
            : ValidatedMedialComplexProfile{};
    const auto finish_profile = [&]() {
        if (profile == nullptr) {
            return;
        }
        const double accounted =
            (profile->candidate_generation_seconds -
             profile_before.candidate_generation_seconds) +
            (profile->point_containment_winding_seconds -
             profile_before.point_containment_winding_seconds) +
            (profile->surface_distance_queries_seconds -
             profile_before.surface_distance_queries_seconds) +
            (profile->segment_intersections_seconds -
             profile_before.segment_intersections_seconds) +
            (profile->pole_support_propagation_seconds -
             profile_before.pole_support_propagation_seconds);
        profile->topology_completion_seconds += std::max(
            0.0,
            detail::medial_profile_elapsed_seconds(total_start) -
                accounted
        );
    };

    MedialComplex result;
    const auto faces = build_interior_voronoi_faces(
        delaunay,
        surface_mesh,
        profile
    );
    result.source_polygon_count = faces.size();
    if (faces.empty()) {
        finish_profile();
        return result;
    }

    const auto support_start = detail::MedialProfileClock::now();
    std::vector<double> tetrahedron_support(
        delaunay.tetrahedron_count(),
        0.0
    );
    std::vector<double> sample_support(delaunay.point_count(), 0.0);
    for (const MedialCandidate& pole : pole_selection.poles) {
        if (!pole.validated) {
            continue;
        }
        if (pole.source_tetrahedron < tetrahedron_support.size()) {
            tetrahedron_support[pole.source_tetrahedron] = std::max(
                tetrahedron_support[pole.source_tetrahedron],
                pole.confidence
            );
        }
        for (std::size_t sample : pole.source_samples) {
            if (sample < sample_support.size()) {
                sample_support[sample] =
                    std::max(sample_support[sample], pole.confidence);
            }
        }
    }

    std::vector<double> support(faces.size(), 0.0);
    std::vector<bool> geometry_is_eligible(faces.size(), false);
    for (std::size_t face_index = 0; face_index < faces.size(); ++face_index) {
        const VoronoiFaceCandidate& face = faces[face_index];
        geometry_is_eligible[face_index] =
            face.vertices.size() >= 3 &&
            face.vertices.size() == face.source_tetrahedra.size();
        if (!geometry_is_eligible[face_index]) {
            continue;
        }

        for (std::size_t source : face.source_tetrahedra) {
            if (source < tetrahedron_support.size()) {
                support[face_index] =
                    std::max(support[face_index], tetrahedron_support[source]);
            }
        }

        const std::size_t first =
            static_cast<std::size_t>(face.contact_sites[0]);
        const std::size_t second =
            static_cast<std::size_t>(face.contact_sites[1]);
        if (first < sample_support.size() &&
            second < sample_support.size() &&
            sample_support[first] > 0.0 &&
            sample_support[second] > 0.0) {
            support[face_index] = std::max(
                support[face_index],
                0.5 * (sample_support[first] + sample_support[second])
            );
        }
    }

    // Grow a small number of rings across actual shared Voronoi edges. This
    // stitches a local sheet around each validated pole without allowing the
    // pole evidence to flood the entire unfiltered Voronoi complex.
    const std::size_t support_ring_count =
        resolved_medial_support_ring_count(
            delaunay.point_count(),
            options
        );
    for (std::size_t ring = 0; ring < support_ring_count; ++ring) {
        std::vector<double> next_support = support;
        for (std::size_t candidate = 0;
             candidate < faces.size();
             ++candidate) {
            if (!geometry_is_eligible[candidate] ||
                support[candidate] > 0.0) {
                continue;
            }
            for (std::size_t accepted = 0;
                 accepted < faces.size();
                 ++accepted) {
                if (support[accepted] <= 0.0 ||
                    !detail::shares_contact_site(
                        faces[candidate],
                        faces[accepted]) ||
                    detail::shared_source_tetrahedron_count(
                        faces[candidate],
                        faces[accepted]) < 2) {
                    continue;
                }
                next_support[candidate] = std::max(
                    next_support[candidate],
                    support[accepted] * options.propagated_support_decay
                );
            }
        }
        support.swap(next_support);
    }
    if (profile != nullptr) {
        profile->pole_support_propagation_seconds +=
            detail::medial_profile_elapsed_seconds(support_start);
    }

    std::vector<bool> polygon_is_inside(faces.size(), false);
    std::vector<bool> polygon_is_kept(faces.size(), false);
    for (std::size_t face_index = 0;
         face_index < faces.size();
         ++face_index) {
        if (!geometry_is_eligible[face_index]) {
            continue;
        }
        polygon_is_inside[face_index] =
            surface_mesh == nullptr ||
            polygon_fan_inside_mesh(
                *surface_mesh,
                faces[face_index].vertices,
                profile
            );
        if (!polygon_is_inside[face_index]) {
            continue;
        }
        const double angle =
            faces[face_index].contact_angle_degrees;
        const bool has_strong_medial_evidence =
            support[face_index] > 0.0 ||
            angle >= options.strong_contact_angle_degrees;
        polygon_is_kept[face_index] =
            has_strong_medial_evidence &&
            angle >= options.minimum_contact_angle_degrees;
    }
    const detail::CandidatePolygonGapRepair polygon_gap_repair =
        detail::restore_enclosed_polygon_gaps(
            faces,
            geometry_is_eligible,
            polygon_is_inside,
            polygon_is_kept
        );
    const detail::CandidatePolygonStratumCompletion stratum_completion =
        detail::complete_seeded_polygon_strata(
            faces,
            geometry_is_eligible,
            polygon_is_inside,
            polygon_is_kept
        );

    Vec3 lower{};
    Vec3 upper{};
    if (!delaunay.points().empty()) {
        lower = delaunay.points().front();
        upper = delaunay.points().front();
        for (const Vec3& point : delaunay.points()) {
            lower.x = std::min(lower.x, point.x);
            lower.y = std::min(lower.y, point.y);
            lower.z = std::min(lower.z, point.z);
            upper.x = std::max(upper.x, point.x);
            upper.y = std::max(upper.y, point.y);
            upper.z = std::max(upper.z, point.z);
        }
    }
    const double weld_tolerance =
        std::max(1e-12, norm(upper - lower) * 1e-10);
    const double weld_tolerance_squared = weld_tolerance * weld_tolerance;
    std::vector<std::size_t> vertex_for_tetrahedron(
        delaunay.tetrahedron_count(),
        std::numeric_limits<std::size_t>::max()
    );
    const auto measure_radius = [&](const Vec3& position,
                                    std::size_t source_tetrahedron) {
        if (surface_mesh != nullptr && !surface_mesh->faces.empty()) {
            std::vector<SurfaceContact> contacts;
            if (profile != nullptr) {
                const auto distance_start =
                    detail::MedialProfileClock::now();
                contacts = nearest_surface_contacts(
                    *surface_mesh,
                    position,
                    0.0
                );
                profile->surface_distance_queries_seconds +=
                    detail::medial_profile_elapsed_seconds(
                        distance_start
                    );
            } else {
                contacts = nearest_surface_contacts(
                    *surface_mesh,
                    position,
                    0.0
                );
            }
            if (!contacts.empty()) {
                return contacts.front().distance;
            }
        }
        if (source_tetrahedron < delaunay.tetrahedron_count()) {
            const Tetrahedron& tetrahedron =
                delaunay.tetrahedra()[source_tetrahedron];
            const int site = tetrahedron.vertices[0];
            if (site >= 0 &&
                static_cast<std::size_t>(site) < delaunay.point_count()) {
                return norm(
                    position -
                    delaunay.points()[static_cast<std::size_t>(site)]
                );
            }
        }
        return 0.0;
    };
    struct VertexFeature {
        double resolution{0.0};
        double local_feature_size{0.0};
        double sampling_density{0.0};
        bool available{false};
    };
    const auto feature_for_tetrahedron =
        [&](std::size_t source_tetrahedron) {
            VertexFeature feature;
            if (surface_feature_field == nullptr ||
                surface_feature_field->resolutions.size() !=
                    delaunay.point_count() ||
                surface_feature_field->local_feature_sizes.size() !=
                    delaunay.point_count() ||
                surface_feature_field->sampling_densities.size() !=
                    delaunay.point_count() ||
                source_tetrahedron >= delaunay.tetrahedron_count()) {
                return feature;
            }
            feature.local_feature_size =
                std::numeric_limits<double>::infinity();
            const Tetrahedron& tetrahedron =
                delaunay.tetrahedra()[source_tetrahedron];
            for (int site : tetrahedron.vertices) {
                if (site < 0 ||
                    static_cast<std::size_t>(site) >=
                        delaunay.point_count()) {
                    continue;
                }
                const std::size_t sample =
                    static_cast<std::size_t>(site);
                feature.resolution = std::max(
                    feature.resolution,
                    surface_feature_field->resolutions[sample]
                );
                feature.local_feature_size = std::min(
                    feature.local_feature_size,
                    surface_feature_field->local_feature_sizes[sample]
                );
                feature.sampling_density = std::max(
                    feature.sampling_density,
                    surface_feature_field->sampling_densities[sample]
                );
                feature.available = true;
            }
            if (!feature.available ||
                !std::isfinite(feature.local_feature_size)) {
                feature = VertexFeature{};
            }
            return feature;
        };

    const auto weld_vertex = [&](const Vec3& position,
                                 std::size_t source_tetrahedron,
                                 MedialComplex& complex) {
        const auto record_source = [&](std::size_t vertex) {
            auto& sources = complex.vertex_source_tetrahedra[vertex];
            if (std::find(
                    sources.begin(),
                    sources.end(),
                    source_tetrahedron) != sources.end()) {
                return;
            }
            sources.push_back(source_tetrahedron);
            const VertexFeature feature =
                feature_for_tetrahedron(source_tetrahedron);
            if (feature.available) {
                complex.vertex_surface_resolutions[vertex] = std::max(
                    complex.vertex_surface_resolutions[vertex],
                    feature.resolution
                );
                const double existing_lfs =
                    complex.vertex_local_feature_sizes[vertex];
                complex.vertex_local_feature_sizes[vertex] =
                    existing_lfs > 0.0
                        ? std::min(
                            existing_lfs,
                            feature.local_feature_size
                        )
                        : feature.local_feature_size;
                complex.vertex_sampling_densities[vertex] = std::max(
                    complex.vertex_sampling_densities[vertex],
                    feature.sampling_density
                );
            }
        };
        if (source_tetrahedron < vertex_for_tetrahedron.size()) {
            const std::size_t existing =
                vertex_for_tetrahedron[source_tetrahedron];
            if (existing != std::numeric_limits<std::size_t>::max()) {
                record_source(existing);
                return existing;
            }
        }
        for (std::size_t index = 0; index < complex.vertices.size(); ++index) {
            if (squared_norm(complex.vertices[index] - position) <=
                weld_tolerance_squared) {
                if (source_tetrahedron < vertex_for_tetrahedron.size()) {
                    vertex_for_tetrahedron[source_tetrahedron] = index;
                }
                record_source(index);
                return index;
            }
        }
        const std::size_t index = complex.vertices.size();
        const double radius =
            measure_radius(position, source_tetrahedron);
        const VertexFeature feature =
            feature_for_tetrahedron(source_tetrahedron);
        complex.vertices.push_back(position);
        complex.vertex_radii.push_back(radius);
        complex.vertex_surface_resolutions.push_back(
            feature.available ? feature.resolution : 0.0
        );
        complex.vertex_local_feature_sizes.push_back(
            feature.available ? feature.local_feature_size : radius
        );
        complex.vertex_sampling_densities.push_back(
            feature.available ? feature.sampling_density : 0.0
        );
        complex.vertex_source_tetrahedra.push_back({source_tetrahedron});
        if (source_tetrahedron < vertex_for_tetrahedron.size()) {
            vertex_for_tetrahedron[source_tetrahedron] = index;
        }
        return index;
    };

    std::vector<bool> initially_kept;
    std::vector<std::size_t> triangle_source_polygons;
    std::vector<double> polygon_confidences(faces.size(), 0.0);
    for (std::size_t face_index = 0; face_index < faces.size(); ++face_index) {
        const VoronoiFaceCandidate& face = faces[face_index];
        const double angle_score = std::clamp(
            (face.contact_angle_degrees -
             options.minimum_contact_angle_degrees) /
                std::max(
                    1.0,
                    180.0 - options.minimum_contact_angle_degrees),
            0.0,
            1.0
        );
        const double confidence =
            0.75 * support[face_index] + 0.25 * angle_score;
        polygon_confidences[face_index] = confidence;
        if (!geometry_is_eligible[face_index]) {
            continue;
        }

        std::vector<std::size_t> polygon;
        polygon.reserve(face.vertices.size());
        for (std::size_t vertex = 0; vertex < face.vertices.size(); ++vertex) {
            polygon.push_back(weld_vertex(
                face.vertices[vertex],
                face.source_tetrahedra[vertex],
                result
            ));
        }

        const std::size_t triangle_count_before = result.triangles.size();

        for (std::size_t i = 1; i + 1 < polygon.size(); ++i) {
            const std::array<std::size_t, 3> triangle{{
                polygon[0],
                polygon[i],
                polygon[i + 1]
            }};
            if (triangle[0] == triangle[1] ||
                triangle[1] == triangle[2] ||
                triangle[2] == triangle[0]) {
                continue;
            }
            const Vec3& a = result.vertices[triangle[0]];
            const Vec3& b = result.vertices[triangle[1]];
            const Vec3& c = result.vertices[triangle[2]];
            if (squared_norm(cross(b - a, c - a)) < 1e-20) {
                continue;
            }
            result.triangles.push_back(triangle);
            result.triangle_support_weights.push_back(
                support[face_index]
            );
            result.triangle_confidences.push_back(confidence);
            result.triangle_contact_angles.push_back(
                face.contact_angle_degrees
            );
            // The polygon-level weak/strong decision above is copied to
            // every fan triangle. Pole support, confidence, and radius
            // continuity remain weights and cannot independently punch a
            // triangle-shaped hole through this atomic 2-cell.
            initially_kept.push_back(
                polygon_is_kept[face_index]
            );
            triangle_source_polygons.push_back(face_index);
        }
        (void)triangle_count_before;
    }

    detail::orient_and_label_sheet_components(result);
    result.triangle_radius_jumps.reserve(result.triangles.size());
    for (const auto& triangle : result.triangles) {
        result.triangle_radius_jumps.push_back(
            detail::triangle_max_radius_jump(result, triangle)
        );
    }
    const detail::CandidateGapRepair gap_repair =
        detail::restore_enclosed_candidate_gaps(
            result,
            initially_kept,
            options
        );

    MedialComplex retained = result;
    retained.triangles.clear();
    retained.triangle_support_weights.clear();
    retained.triangle_confidences.clear();
    retained.triangle_contact_angles.clear();
    retained.triangle_radius_jumps.clear();
    retained.triangle_radius_jump_thresholds.clear();
    retained.triangle_components.clear();
    retained.boundary_edges.clear();
    retained.seam_edges.clear();
    retained.junction_vertices.clear();
    retained.termination_edges.clear();
    retained.artificial_boundary_edges.clear();
    retained.boundary_loops.clear();
    retained.rejected_face_centers.clear();
    retained.rejected_face_confidences.clear();
    retained.topology_restored_candidate_triangles.clear();
    retained.topology_restored_candidate_patch_count =
        polygon_gap_repair.repaired_patch_count +
        stratum_completion.completed_stratum_count +
        gap_repair.repaired_patch_count;
    retained.component_count = 0;

    std::vector<bool> polygon_accepted(faces.size(), false);
    for (std::size_t triangle = 0;
         triangle < result.triangles.size();
         ++triangle) {
        if (!initially_kept[triangle]) {
            continue;
        }
        retained.triangles.push_back(result.triangles[triangle]);
        retained.triangle_support_weights.push_back(
            result.triangle_support_weights[triangle]
        );
        retained.triangle_confidences.push_back(
            result.triangle_confidences[triangle]
        );
        retained.triangle_contact_angles.push_back(
            result.triangle_contact_angles[triangle]
        );
        retained.triangle_radius_jumps.push_back(
            result.triangle_radius_jumps[triangle]
        );
        const bool polygon_was_repaired =
            triangle < triangle_source_polygons.size() &&
            triangle_source_polygons[triangle] <
                polygon_gap_repair.repaired.size() &&
            (polygon_gap_repair.repaired[
                 triangle_source_polygons[triangle]
             ] ||
             (triangle_source_polygons[triangle] <
                  stratum_completion.repaired.size() &&
              stratum_completion.repaired[
                  triangle_source_polygons[triangle]
              ]));
        if (polygon_was_repaired ||
            (triangle < gap_repair.repaired.size() &&
             gap_repair.repaired[triangle])) {
            retained.topology_restored_candidate_triangles.push_back(
                result.triangles[triangle]
            );
        }
        if (triangle < triangle_source_polygons.size()) {
            polygon_accepted[
                triangle_source_polygons[triangle]
            ] = true;
        }
    }
    for (std::size_t face_index = 0;
         face_index < faces.size();
         ++face_index) {
        if (polygon_accepted[face_index]) {
            ++retained.accepted_polygon_count;
        } else {
            retained.rejected_face_centers.push_back(
                detail::polygon_center(faces[face_index])
            );
            retained.rejected_face_confidences.push_back(
                polygon_confidences[face_index]
            );
        }
    }

    if (!retained.triangles.empty()) {
        detail::orient_and_label_sheet_components(retained);
        detail::classify_complex_topology(
            retained,
            options.minimum_contact_angle_degrees +
                options.termination_angle_margin_degrees
        );
    }
    finish_profile();
    return retained;
}

}  // namespace medial_axis_3d
