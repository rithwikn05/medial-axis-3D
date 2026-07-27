#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <vector>

#include "filters.h"
#include "medial_complex.h"

namespace medial_axis_3d {

struct ResolutionMedialComplex {
    std::size_t sample_count{0};
    MedialComplex complex;
    std::vector<MedialComponentMetrics> component_metrics;
};

struct CrossResolutionStabilityOptions {
    double maximum_lfs_normalized_vertex_distance{1.0};
    double maximum_relative_radius_difference{0.50};
    double minimum_component_overlap{0.15};
    double minimum_component_normal_alignment{0.40};
    double minimum_component_match_score{0.35};
    double maximum_lfs_normalized_triangle_distance{1.0};
    double minimum_triangle_normal_alignment{0.40};
    double minimum_stability{0.50};
    // Cross-resolution matching is a confidence diagnostic by default.
    // Coarse runs can be orders of magnitude smaller than a vertex-preserving
    // displayed run, so using them as a deletion mask destroys valid detailed
    // sheets. Explicit callers may still opt into destructive pruning.
    bool remove_unstable_components{false};
    bool repair_enclosed_holes{true};
    std::size_t maximum_repaired_hole_triangles{24};
    double maximum_repaired_hole_area_fraction{0.05};
    double minimum_repaired_hole_boundary_stability{0.50};
    double minimum_boundary_complexity_ratio_for_repair{1.0};
    bool preserve_supported_component_topology{true};
};

struct CrossResolutionStabilityResult {
    MedialComplex retained;
    std::vector<double> vertex_stabilities;
    std::vector<double> retained_triangle_stabilities;
    std::vector<double> retained_component_stabilities;
    std::vector<std::size_t> retained_triangle_source_components;
    std::vector<double> retained_triangle_repair_flags;

    std::vector<std::array<std::size_t, 3>> removed_triangles;
    std::vector<double> removed_triangle_stabilities;
    std::vector<std::size_t> removed_triangle_source_components;

    std::vector<std::array<std::size_t, 3>> topology_repaired_triangles;
    std::vector<double> topology_repaired_triangle_stabilities;
    std::vector<std::size_t>
        topology_repaired_triangle_source_components;
    std::size_t topology_repaired_patch_count{0};

    std::vector<double> source_component_stabilities;
    std::vector<std::size_t> compared_sample_counts;
};

namespace detail {

struct StabilityComponentDescriptor {
    Vec3 centroid{};
    Vec3 normal{};
    double area{0.0};
    std::size_t triangle_count{0};
    std::size_t supporting_pole_count{0};
    std::size_t boundary_edge_count{0};
    std::size_t seam_edge_count{0};
    std::size_t junction_count{0};
    std::set<std::size_t> vertices;
};

inline Vec3 stability_triangle_centroid(
    const MedialComplex& complex,
    const std::array<std::size_t, 3>& triangle) {
    return (
        complex.vertices[triangle[0]] +
        complex.vertices[triangle[1]] +
        complex.vertices[triangle[2]]
    ) / 3.0;
}

inline Vec3 stability_triangle_normal(
    const MedialComplex& complex,
    const std::array<std::size_t, 3>& triangle) {
    return normalized(cross(
        complex.vertices[triangle[1]] -
            complex.vertices[triangle[0]],
        complex.vertices[triangle[2]] -
            complex.vertices[triangle[0]]
    ));
}

inline double stability_triangle_average(
    const std::vector<double>& values,
    const std::array<std::size_t, 3>& triangle,
    double fallback) {
    if (values.empty() ||
        triangle[0] >= values.size() ||
        triangle[1] >= values.size() ||
        triangle[2] >= values.size()) {
        return fallback;
    }
    return (
        values[triangle[0]] +
        values[triangle[1]] +
        values[triangle[2]]
    ) / 3.0;
}

inline double stability_relative_difference(double first,
                                            double second) {
    return std::fabs(first - second) /
        std::max({std::fabs(first), std::fabs(second), 1e-12});
}

inline double stability_vertex_scale(const MedialComplex& first,
                                     std::size_t first_vertex,
                                     const MedialComplex& second,
                                     std::size_t second_vertex) {
    double first_lfs =
        first_vertex < first.vertex_local_feature_sizes.size()
            ? first.vertex_local_feature_sizes[first_vertex]
            : 0.0;
    double second_lfs =
        second_vertex < second.vertex_local_feature_sizes.size()
            ? second.vertex_local_feature_sizes[second_vertex]
            : 0.0;
    if (first_lfs <= 0.0 && first_vertex < first.vertex_radii.size()) {
        first_lfs = first.vertex_radii[first_vertex];
    }
    if (second_lfs <= 0.0 && second_vertex < second.vertex_radii.size()) {
        second_lfs = second.vertex_radii[second_vertex];
    }
    return std::max(1e-12, 0.5 * (first_lfs + second_lfs));
}

inline std::vector<std::set<std::size_t>>
stability_vertex_components(const MedialComplex& complex) {
    std::vector<std::set<std::size_t>> result(complex.vertices.size());
    for (std::size_t triangle = 0;
         triangle < complex.triangles.size();
         ++triangle) {
        if (triangle >= complex.triangle_components.size()) {
            continue;
        }
        const std::size_t component =
            complex.triangle_components[triangle];
        for (std::size_t vertex : complex.triangles[triangle]) {
            if (vertex < result.size()) {
                result[vertex].insert(component);
            }
        }
    }
    return result;
}

inline std::vector<StabilityComponentDescriptor>
stability_component_descriptors(
    const ResolutionMedialComplex& run) {
    const MedialComplex& complex = run.complex;
    std::vector<StabilityComponentDescriptor> result(
        complex.component_count
    );
    for (std::size_t triangle = 0;
         triangle < complex.triangles.size();
         ++triangle) {
        if (triangle >= complex.triangle_components.size()) {
            continue;
        }
        const std::size_t component =
            complex.triangle_components[triangle];
        if (component >= result.size()) {
            continue;
        }
        const auto& face = complex.triangles[triangle];
        const double area = medial_triangle_area(complex, face);
        StabilityComponentDescriptor& descriptor = result[component];
        descriptor.centroid =
            descriptor.centroid +
            stability_triangle_centroid(complex, face) * area;
        descriptor.normal =
            descriptor.normal +
            stability_triangle_normal(complex, face) * area;
        descriptor.area += area;
        ++descriptor.triangle_count;
        descriptor.vertices.insert(face.begin(), face.end());
    }

    const auto incidence =
        triangle_edge_incidence(complex.triangles);
    for (const auto& [edge, incident] : incidence) {
        std::set<std::size_t> components;
        for (std::size_t triangle : incident) {
            if (triangle < complex.triangle_components.size()) {
                components.insert(complex.triangle_components[triangle]);
            }
        }
        for (std::size_t component : components) {
            if (component >= result.size()) {
                continue;
            }
            if (incident.size() == 1) {
                ++result[component].boundary_edge_count;
            } else if (incident.size() > 2) {
                ++result[component].seam_edge_count;
            }
        }
    }

    for (std::size_t junction : complex.junction_vertices) {
        if (junction >= complex.vertices.size()) {
            continue;
        }
        for (std::size_t triangle = 0;
             triangle < complex.triangles.size();
             ++triangle) {
            if (std::find(
                    complex.triangles[triangle].begin(),
                    complex.triangles[triangle].end(),
                    junction) == complex.triangles[triangle].end()) {
                continue;
            }
            const std::size_t component =
                complex.triangle_components[triangle];
            if (component < result.size()) {
                ++result[component].junction_count;
            }
        }
    }

    for (std::size_t component = 0;
         component < result.size();
         ++component) {
        StabilityComponentDescriptor& descriptor = result[component];
        if (descriptor.area > 0.0) {
            descriptor.centroid =
                descriptor.centroid / descriptor.area;
        }
        descriptor.normal = normalized(descriptor.normal);
        if (component < run.component_metrics.size()) {
            descriptor.supporting_pole_count =
                run.component_metrics[component].supporting_pole_count;
        }
    }
    return result;
}

inline double stability_count_similarity(std::size_t first,
                                         std::size_t second) {
    return static_cast<double>(std::min(first, second) + 1) /
        static_cast<double>(std::max(first, second) + 1);
}

inline std::vector<int> match_stability_vertices(
    const MedialComplex& base,
    const MedialComplex& comparison,
    const CrossResolutionStabilityOptions& options) {
    std::vector<int> result(base.vertices.size(), -1);
    for (std::size_t base_vertex = 0;
         base_vertex < base.vertices.size();
         ++base_vertex) {
        double best_score = -std::numeric_limits<double>::infinity();
        for (std::size_t comparison_vertex = 0;
             comparison_vertex < comparison.vertices.size();
             ++comparison_vertex) {
            const double scale = stability_vertex_scale(
                base,
                base_vertex,
                comparison,
                comparison_vertex
            );
            const double normalized_distance = norm(
                base.vertices[base_vertex] -
                comparison.vertices[comparison_vertex]
            ) / scale;
            const double base_radius =
                base_vertex < base.vertex_radii.size()
                    ? base.vertex_radii[base_vertex]
                    : 0.0;
            const double comparison_radius =
                comparison_vertex < comparison.vertex_radii.size()
                    ? comparison.vertex_radii[comparison_vertex]
                    : 0.0;
            const double radius_difference = stability_relative_difference(
                base_radius,
                comparison_radius
            );
            if (normalized_distance >
                    options.maximum_lfs_normalized_vertex_distance ||
                radius_difference >
                    options.maximum_relative_radius_difference) {
                continue;
            }
            const double score =
                1.0 -
                0.65 * normalized_distance /
                    options.maximum_lfs_normalized_vertex_distance -
                0.35 * radius_difference /
                    options.maximum_relative_radius_difference;
            if (score > best_score) {
                best_score = score;
                result[base_vertex] =
                    static_cast<int>(comparison_vertex);
            }
        }
    }
    return result;
}

inline std::vector<int> match_stability_components(
    const ResolutionMedialComplex& base,
    const ResolutionMedialComplex& comparison,
    const std::vector<int>& vertex_matches,
    const CrossResolutionStabilityOptions& options) {
    const auto base_descriptors =
        stability_component_descriptors(base);
    const auto comparison_descriptors =
        stability_component_descriptors(comparison);
    const auto comparison_vertex_components =
        stability_vertex_components(comparison.complex);
    std::vector<int> result(base_descriptors.size(), -1);

    for (std::size_t base_component = 0;
         base_component < base_descriptors.size();
         ++base_component) {
        const auto& first = base_descriptors[base_component];
        double best_score = -1.0;
        for (std::size_t comparison_component = 0;
             comparison_component < comparison_descriptors.size();
             ++comparison_component) {
            const auto& second =
                comparison_descriptors[comparison_component];
            std::size_t overlap_count = 0;
            for (std::size_t vertex : first.vertices) {
                if (vertex >= vertex_matches.size() ||
                    vertex_matches[vertex] < 0) {
                    continue;
                }
                const std::size_t matched =
                    static_cast<std::size_t>(vertex_matches[vertex]);
                if (matched < comparison_vertex_components.size() &&
                    comparison_vertex_components[matched].count(
                        comparison_component) != 0) {
                    ++overlap_count;
                }
            }
            const double overlap = first.vertices.empty()
                ? 0.0
                : static_cast<double>(overlap_count) /
                    static_cast<double>(first.vertices.size());
            const double normal_alignment = std::fabs(
                dot(first.normal, second.normal)
            );
            if (overlap < options.minimum_component_overlap ||
                normal_alignment <
                    options.minimum_component_normal_alignment) {
                continue;
            }
            const double pole_similarity = stability_count_similarity(
                first.supporting_pole_count,
                second.supporting_pole_count
            );
            const double topology_similarity =
                (
                    stability_count_similarity(
                        first.boundary_edge_count,
                        second.boundary_edge_count
                    ) +
                    stability_count_similarity(
                        first.seam_edge_count,
                        second.seam_edge_count
                    ) +
                    stability_count_similarity(
                        first.junction_count,
                        second.junction_count
                    )
                ) / 3.0;
            const double score =
                0.50 * overlap +
                0.20 * normal_alignment +
                0.15 * pole_similarity +
                0.15 * topology_similarity;
            if (score >= options.minimum_component_match_score &&
                score > best_score) {
                best_score = score;
                result[base_component] =
                    static_cast<int>(comparison_component);
            }
        }
    }
    return result;
}

inline bool stability_triangle_supported(
    const MedialComplex& base,
    std::size_t base_triangle,
    const MedialComplex& comparison,
    std::size_t comparison_component,
    const std::vector<int>& vertex_matches,
    const std::vector<std::set<std::size_t>>&
        comparison_vertex_components,
    const CrossResolutionStabilityOptions& options) {
    const auto& triangle = base.triangles[base_triangle];
    std::size_t directly_matched = 0;
    for (std::size_t vertex : triangle) {
        if (vertex >= vertex_matches.size() ||
            vertex_matches[vertex] < 0) {
            continue;
        }
        const std::size_t matched =
            static_cast<std::size_t>(vertex_matches[vertex]);
        if (matched < comparison_vertex_components.size() &&
            comparison_vertex_components[matched].count(
                comparison_component) != 0) {
            ++directly_matched;
        }
    }
    if (directly_matched >= 2) {
        return true;
    }

    const Vec3 base_centroid =
        stability_triangle_centroid(base, triangle);
    const Vec3 base_normal =
        stability_triangle_normal(base, triangle);
    const double base_lfs = stability_triangle_average(
        base.vertex_local_feature_sizes,
        triangle,
        stability_triangle_average(base.vertex_radii, triangle, 1.0)
    );
    const double base_radius =
        stability_triangle_average(base.vertex_radii, triangle, 0.0);

    for (std::size_t candidate = 0;
         candidate < comparison.triangles.size();
         ++candidate) {
        if (candidate >= comparison.triangle_components.size() ||
            comparison.triangle_components[candidate] !=
                comparison_component) {
            continue;
        }
        const auto& comparison_triangle =
            comparison.triangles[candidate];
        const double comparison_lfs = stability_triangle_average(
            comparison.vertex_local_feature_sizes,
            comparison_triangle,
            stability_triangle_average(
                comparison.vertex_radii,
                comparison_triangle,
                1.0
            )
        );
        const double scale = std::max(
            1e-12,
            0.5 * (base_lfs + comparison_lfs)
        );
        const Vec3 closest = closest_point_on_triangle(
            base_centroid,
            comparison.vertices[comparison_triangle[0]],
            comparison.vertices[comparison_triangle[1]],
            comparison.vertices[comparison_triangle[2]]
        );
        const double normalized_distance =
            norm(base_centroid - closest) / scale;
        if (normalized_distance >
                options.maximum_lfs_normalized_triangle_distance) {
            continue;
        }
        const double normal_alignment = std::fabs(dot(
            base_normal,
            stability_triangle_normal(
                comparison,
                comparison_triangle
            )
        ));
        if (normal_alignment <
                options.minimum_triangle_normal_alignment) {
            continue;
        }
        const double comparison_radius = stability_triangle_average(
            comparison.vertex_radii,
            comparison_triangle,
            0.0
        );
        if (stability_relative_difference(
                base_radius,
                comparison_radius) <=
            options.maximum_relative_radius_difference) {
            return true;
        }
    }
    return false;
}

struct StabilityTopologyRepair {
    std::vector<bool> repaired;
    std::size_t repaired_patch_count{0};
};

inline StabilityTopologyRepair repair_stability_holes(
    const MedialComplex& complex,
    const std::vector<double>& triangle_stabilities,
    const std::vector<double>& component_stabilities,
    std::vector<bool>& keep,
    const CrossResolutionStabilityOptions& options) {
    StabilityTopologyRepair result;
    result.repaired.assign(complex.triangles.size(), false);
    if (!options.repair_enclosed_holes ||
        complex.triangles.empty() ||
        keep.size() != complex.triangles.size() ||
        triangle_stabilities.size() != complex.triangles.size()) {
        return result;
    }

    const auto incidence = triangle_edge_incidence(complex.triangles);
    std::vector<std::vector<std::size_t>> rejected_adjacency(
        complex.triangles.size()
    );
    for (const auto& [edge, incident] : incidence) {
        // A medial seam is not a hole boundary. Only use ordinary manifold
        // sheet adjacency when grouping rejected patches.
        if (incident.size() != 2) {
            continue;
        }
        const std::size_t first = incident[0];
        const std::size_t second = incident[1];
        const bool same_component =
            first < complex.triangle_components.size() &&
            second < complex.triangle_components.size() &&
            complex.triangle_components[first] ==
                complex.triangle_components[second];
        if (!keep[first] && !keep[second] && same_component) {
            rejected_adjacency[first].push_back(second);
            rejected_adjacency[second].push_back(first);
        }
    }

    std::vector<double> component_areas(complex.component_count, 0.0);
    for (std::size_t triangle = 0;
         triangle < complex.triangles.size();
         ++triangle) {
        if (triangle < complex.triangle_components.size() &&
            complex.triangle_components[triangle] <
                component_areas.size()) {
            component_areas[complex.triangle_components[triangle]] +=
                medial_triangle_area(
                    complex,
                    complex.triangles[triangle]
                );
        }
    }

    std::vector<bool> visited(complex.triangles.size(), false);
    for (std::size_t seed = 0;
         seed < complex.triangles.size();
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

        if (patch.empty()) {
            continue;
        }
        const std::size_t source_component =
            seed < complex.triangle_components.size()
                ? complex.triangle_components[seed]
                : complex.component_count;
        if (source_component >= component_areas.size() ||
            component_areas[source_component] <= 0.0) {
            continue;
        }
        const bool component_supported =
            source_component < component_stabilities.size() &&
            component_stabilities[source_component] >=
                options.minimum_stability;
        const bool preserves_matched_component =
            options.preserve_supported_component_topology &&
            component_supported;
        if (!preserves_matched_component &&
            patch.size() >
                options.maximum_repaired_hole_triangles) {
            continue;
        }

        std::set<std::size_t> patch_set(patch.begin(), patch.end());
        std::set<std::size_t> boundary_neighbors;
        bool enclosed_by_retained_sheet = true;
        bool touches_medial_seam = false;
        std::size_t boundary_edge_count = 0;
        double patch_area = 0.0;
        double retained_interface_length = 0.0;
        double original_boundary_length = 0.0;
        for (std::size_t triangle : patch) {
            patch_area += medial_triangle_area(
                complex,
                complex.triangles[triangle]
            );
            const auto& face = complex.triangles[triangle];
            for (std::size_t edge_index = 0;
                 edge_index < 3;
                 ++edge_index) {
                const ComplexEdge edge = normalized_complex_edge(
                    face[edge_index],
                    face[(edge_index + 1) % 3]
                );
                const auto found = incidence.find(edge);
                if (found == incidence.end()) {
                    enclosed_by_retained_sheet = false;
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

                ++boundary_edge_count;
                const double edge_length =
                    medial_edge_length(complex, edge);
                if (found->second.size() == 1) {
                    original_boundary_length += edge_length;
                    enclosed_by_retained_sheet = false;
                    continue;
                }
                if (found->second.size() > 2) {
                    touches_medial_seam = true;
                    enclosed_by_retained_sheet = false;
                }
                bool has_retained_neighbor = false;
                for (std::size_t neighbor : found->second) {
                    if (patch_set.count(neighbor) != 0 ||
                        !keep[neighbor] ||
                        neighbor >=
                            complex.triangle_components.size() ||
                        complex.triangle_components[neighbor] !=
                            source_component) {
                        continue;
                    }
                    has_retained_neighbor = true;
                    boundary_neighbors.insert(neighbor);
                }
                if (has_retained_neighbor) {
                    retained_interface_length += edge_length;
                } else {
                    enclosed_by_retained_sheet = false;
                }
            }
        }

        const double area_fraction =
            patch_area / component_areas[source_component];
        const bool closes_enclosed_hole =
            enclosed_by_retained_sheet &&
            boundary_edge_count >= 3 &&
            boundary_neighbors.size() >= 3;
        const bool avoids_boundary_notch =
            component_supported &&
            boundary_neighbors.size() >= 2 &&
            retained_interface_length >
                options.minimum_boundary_complexity_ratio_for_repair *
                    original_boundary_length;
        const bool preserves_supported_seam =
            component_supported &&
            touches_medial_seam &&
            boundary_neighbors.size() >= 2;
        if ((!closes_enclosed_hole &&
             !avoids_boundary_notch &&
             !preserves_supported_seam &&
             !preserves_matched_component) ||
            (!preserves_matched_component &&
             area_fraction >
                options.maximum_repaired_hole_area_fraction)) {
            continue;
        }

        double boundary_stability = 0.0;
        for (std::size_t neighbor : boundary_neighbors) {
            boundary_stability += triangle_stabilities[neighbor];
        }
        boundary_stability /=
            static_cast<double>(boundary_neighbors.size());
        if (boundary_stability <
                options.minimum_repaired_hole_boundary_stability) {
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

inline CrossResolutionStabilityResult
analyze_cross_resolution_stability(
    const std::vector<ResolutionMedialComplex>& runs,
    std::size_t base_run = 0,
    const CrossResolutionStabilityOptions& options = {}) {
    CrossResolutionStabilityResult result;
    if (runs.empty() || base_run >= runs.size()) {
        return result;
    }
    const ResolutionMedialComplex& base = runs[base_run];
    result.retained = base.complex;
    result.compared_sample_counts.reserve(runs.size());
    for (const ResolutionMedialComplex& run : runs) {
        result.compared_sample_counts.push_back(run.sample_count);
    }

    std::vector<std::size_t> vertex_support(
        base.complex.vertices.size(),
        1
    );
    std::vector<std::size_t> triangle_support(
        base.complex.triangles.size(),
        1
    );
    std::vector<std::size_t> component_support(
        base.complex.component_count,
        1
    );

    for (std::size_t run_index = 0;
         run_index < runs.size();
         ++run_index) {
        if (run_index == base_run || runs[run_index].complex.triangles.empty()) {
            continue;
        }
        const ResolutionMedialComplex& comparison = runs[run_index];
        const auto vertex_matches = detail::match_stability_vertices(
            base.complex,
            comparison.complex,
            options
        );
        for (std::size_t vertex = 0;
             vertex < vertex_matches.size();
             ++vertex) {
            if (vertex_matches[vertex] >= 0) {
                ++vertex_support[vertex];
            }
        }

        const auto component_matches =
            detail::match_stability_components(
                base,
                comparison,
                vertex_matches,
                options
            );
        for (std::size_t component = 0;
             component < component_matches.size();
             ++component) {
            if (component_matches[component] >= 0) {
                ++component_support[component];
            }
        }

        const auto comparison_vertex_components =
            detail::stability_vertex_components(comparison.complex);
        for (std::size_t triangle = 0;
             triangle < base.complex.triangles.size();
             ++triangle) {
            if (triangle >= base.complex.triangle_components.size()) {
                continue;
            }
            const std::size_t component =
                base.complex.triangle_components[triangle];
            if (component >= component_matches.size() ||
                component_matches[component] < 0) {
                continue;
            }
            if (detail::stability_triangle_supported(
                    base.complex,
                    triangle,
                    comparison.complex,
                    static_cast<std::size_t>(
                        component_matches[component]
                    ),
                    vertex_matches,
                    comparison_vertex_components,
                    options)) {
                ++triangle_support[triangle];
            }
        }
    }

    const double run_count = static_cast<double>(runs.size());
    result.vertex_stabilities.reserve(vertex_support.size());
    for (std::size_t support : vertex_support) {
        result.vertex_stabilities.push_back(
            static_cast<double>(support) / run_count
        );
    }
    result.source_component_stabilities.reserve(component_support.size());
    for (std::size_t support : component_support) {
        result.source_component_stabilities.push_back(
            static_cast<double>(support) / run_count
        );
    }

    std::vector<double> triangle_stabilities(
        base.complex.triangles.size(),
        0.0
    );
    std::vector<bool> keep(base.complex.triangles.size(), false);
    for (std::size_t triangle = 0;
         triangle < base.complex.triangles.size();
         ++triangle) {
        triangle_stabilities[triangle] =
            static_cast<double>(triangle_support[triangle]) / run_count;
        const std::size_t component =
            triangle < base.complex.triangle_components.size()
                ? base.complex.triangle_components[triangle]
                : base.complex.component_count;
        // Stability removes complete sheet components. Per-triangle
        // stability is retained as a diagnostic, but cannot puncture a
        // component which is supported across resolutions.
        keep[triangle] =
            !options.remove_unstable_components ||
            (component < result.source_component_stabilities.size() &&
             result.source_component_stabilities[component] >=
                 options.minimum_stability);
    }
    detail::StabilityTopologyRepair topology_repair;
    topology_repair.repaired.assign(
        base.complex.triangles.size(),
        false
    );

    result.retained.triangles.clear();
    result.retained.triangle_support_weights.clear();
    result.retained.triangle_confidences.clear();
    result.retained.triangle_contact_angles.clear();
    result.retained.triangle_radius_jumps.clear();
    result.retained.triangle_radius_jump_thresholds.clear();
    result.retained.triangle_components.clear();
    result.retained.boundary_edges.clear();
    result.retained.seam_edges.clear();
    result.retained.junction_vertices.clear();
    result.retained.termination_edges.clear();
    result.retained.artificial_boundary_edges.clear();
    result.retained.boundary_loops.clear();
    result.retained.component_count = 0;

    for (std::size_t triangle = 0;
         triangle < base.complex.triangles.size();
         ++triangle) {
        const double stability = triangle_stabilities[triangle];
        const std::size_t source_component =
            triangle < base.complex.triangle_components.size()
                ? base.complex.triangle_components[triangle]
                : 0;
        if (keep[triangle]) {
            result.retained.triangles.push_back(
                base.complex.triangles[triangle]
            );
            if (triangle <
                    base.complex.triangle_support_weights.size()) {
                result.retained.triangle_support_weights.push_back(
                    base.complex.triangle_support_weights[triangle]
                );
            }
            result.retained.triangle_confidences.push_back(
                base.complex.triangle_confidences[triangle]
            );
            result.retained.triangle_contact_angles.push_back(
                base.complex.triangle_contact_angles[triangle]
            );
            result.retained.triangle_radius_jumps.push_back(
                base.complex.triangle_radius_jumps[triangle]
            );
            result.retained.triangle_radius_jump_thresholds.push_back(
                triangle <
                        base.complex.triangle_radius_jump_thresholds.size()
                    ? base.complex.
                        triangle_radius_jump_thresholds[triangle]
                    : 0.0
            );
            result.retained_triangle_stabilities.push_back(stability);
            result.retained_component_stabilities.push_back(
                source_component <
                        result.source_component_stabilities.size()
                    ? result.source_component_stabilities[source_component]
                    : 0.0
            );
            result.retained_triangle_source_components.push_back(
                source_component
            );
            const bool repaired =
                triangle < topology_repair.repaired.size() &&
                topology_repair.repaired[triangle];
            result.retained_triangle_repair_flags.push_back(
                repaired ? 1.0 : 0.0
            );
            if (repaired) {
                result.topology_repaired_triangles.push_back(
                    base.complex.triangles[triangle]
                );
                result.topology_repaired_triangle_stabilities.push_back(
                    stability
                );
                result.
                    topology_repaired_triangle_source_components.push_back(
                        source_component
                    );
            }
        } else {
            result.removed_triangles.push_back(
                base.complex.triangles[triangle]
            );
            result.removed_triangle_stabilities.push_back(stability);
            result.removed_triangle_source_components.push_back(
                source_component
            );
        }
    }

    if (!result.retained.triangles.empty()) {
        detail::orient_and_label_sheet_components(result.retained);
        detail::classify_complex_topology(result.retained);
    }
    return result;
}

}  // namespace medial_axis_3d
