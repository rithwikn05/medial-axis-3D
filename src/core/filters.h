#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <iterator>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <vector>

#include "medial_complex.h"

namespace medial_axis_3d {

struct FilterSettings {
    double min_radius{0.0};
    double max_radius{std::numeric_limits<double>::infinity()};
};

struct RadiusContinuityFilterOptions {
    bool use_adaptive_lfs_thresholds{true};
    double minimum_relative_radius_jump{0.35};
    double sampling_density_scale{0.60};
    double maximum_relative_radius_jump{0.85};
    bool preserve_lowest_jump_triangle_if_all_removed{true};
    bool preserve_sheet_topology{true};
    std::size_t maximum_topology_repair_triangles{24};
    double maximum_topology_repair_area_fraction{0.05};
};

struct RadiusContinuityFilterResult {
    MedialComplex retained;
    std::vector<std::array<std::size_t, 3>> flagged_triangles;
    std::vector<double> flagged_triangle_radius_jumps;
    std::vector<double> flagged_triangle_radius_jump_thresholds;
    std::vector<std::array<std::size_t, 3>> removed_triangles;
    std::vector<double> removed_triangle_confidences;
    std::vector<double> removed_triangle_contact_angles;
    std::vector<double> removed_triangle_radius_jumps;
    std::vector<double> removed_triangle_radius_jump_thresholds;
    std::vector<std::array<std::size_t, 3>>
        topology_repaired_triangles;
    std::vector<double> topology_repaired_triangle_radius_jumps;
    std::vector<double>
        topology_repaired_triangle_radius_jump_thresholds;
    std::size_t topology_repaired_patch_count{0};
    double maximum_observed_radius_jump{0.0};
    double minimum_applied_threshold{0.0};
    double maximum_applied_threshold{0.0};
};

struct MedialComponentMetrics {
    std::size_t component{0};
    std::size_t triangle_count{0};
    std::size_t vertex_count{0};
    std::size_t supporting_pole_count{0};
    double area{0.0};
    double area_fraction{0.0};
    double mean_confidence{0.0};
    double minimum_confidence{0.0};
    double mean_radius_jump{0.0};
    double maximum_radius_jump{0.0};
    double boundary_length{0.0};
    double artificial_boundary_length{0.0};
    double artificial_boundary_fraction{0.0};
    double seam_length{0.0};
    std::size_t boundary_loop_count{0};
    std::size_t artificial_boundary_loop_count{0};
    bool retained{false};
};

struct MedialComponentFilterOptions {
    // Keep even small pole-supported pieces for cross-resolution matching.
    // Stability analysis, rather than mesh-resolution-dependent size, decides
    // whether a sheet is persistent.
    std::size_t minimum_triangle_count{1};
    double minimum_area_fraction{0.0};
    // Confidence is diagnostic by default. Destructive confidence pruning
    // can remove a complete stratum after topology completion and reopen
    // visible gaps.
    double minimum_mean_confidence{0.0};
    // Propagated support can validate a neighboring stratum even when that
    // stratum does not contain a pole's source tetrahedron. Requiring direct
    // provenance by default cuts supported pieces out of the complex.
    // Callers can still request direct poles explicitly.
    std::size_t minimum_supporting_poles{0};
    // Boundary classification is likewise diagnostic by default. Callers
    // may explicitly enable destructive boundary validation.
    bool require_valid_sheet_boundaries{false};
    // Tiny unresolved loops are commonly caused by floating-point surface
    // classification. Do not discard an otherwise coherent sheet unless the
    // unresolved boundary is a material part of its total boundary.
    double maximum_artificial_boundary_fraction{0.02};
    bool preserve_largest_component_if_all_removed{true};
};

struct MedialComplexFilterResult {
    MedialComplex retained;
    std::vector<std::array<std::size_t, 3>> removed_triangles;
    std::vector<double> removed_triangle_confidences;
    std::vector<double> removed_triangle_contact_angles;
    std::vector<double> removed_triangle_radius_jumps;
    std::vector<double> removed_triangle_radius_jump_thresholds;
    std::vector<std::size_t> removed_triangle_source_components;
    std::vector<MedialComponentMetrics> source_component_metrics;
    std::vector<MedialComponentMetrics> retained_component_metrics;
};

inline RadiusContinuityFilterResult
filter_medial_complex_radius_continuity(
    const MedialComplex& input,
    const RadiusContinuityFilterOptions& options = {}) {
    RadiusContinuityFilterResult result;
    if (input.triangles.empty() ||
        input.vertex_radii.size() != input.vertices.size()) {
        result.retained = input;
        return result;
    }

    std::vector<double> jumps;
    std::vector<double> thresholds;
    jumps.reserve(input.triangles.size());
    thresholds.reserve(input.triangles.size());
    const bool has_adaptive_field =
        options.use_adaptive_lfs_thresholds &&
        input.vertex_sampling_densities.size() == input.vertices.size() &&
        input.vertex_surface_resolutions.size() == input.vertices.size() &&
        std::any_of(
            input.vertex_surface_resolutions.begin(),
            input.vertex_surface_resolutions.end(),
            [](double resolution) {
                return resolution > 0.0;
            }
        );
    result.minimum_applied_threshold =
        std::numeric_limits<double>::infinity();
    for (std::size_t triangle = 0;
         triangle < input.triangles.size();
         ++triangle) {
        const double jump =
            triangle < input.triangle_radius_jumps.size()
                ? input.triangle_radius_jumps[triangle]
                : detail::triangle_max_radius_jump(
                    input,
                    input.triangles[triangle]
                );
        jumps.push_back(jump);
        result.maximum_observed_radius_jump =
            std::max(result.maximum_observed_radius_jump, jump);
        double threshold = options.maximum_relative_radius_jump;
        if (has_adaptive_field) {
            const auto& vertices = input.triangles[triangle];
            const double density = std::max({
                input.vertex_sampling_densities[vertices[0]],
                input.vertex_sampling_densities[vertices[1]],
                input.vertex_sampling_densities[vertices[2]]
            });
            const double lower = std::min(
                options.minimum_relative_radius_jump,
                options.maximum_relative_radius_jump
            );
            const double bounded_density =
                density / (1.0 + std::max(0.0, density));
            threshold = std::clamp(
                options.minimum_relative_radius_jump +
                    options.sampling_density_scale * bounded_density,
                lower,
                options.maximum_relative_radius_jump
            );
        }
        thresholds.push_back(threshold);
        result.minimum_applied_threshold =
            std::min(result.minimum_applied_threshold, threshold);
        result.maximum_applied_threshold =
            std::max(result.maximum_applied_threshold, threshold);
    }

    std::vector<bool> keep(input.triangles.size(), true);
    for (std::size_t triangle = 0; triangle < jumps.size(); ++triangle) {
        if (jumps[triangle] > thresholds[triangle]) {
            result.flagged_triangles.push_back(
                input.triangles[triangle]
            );
            result.flagged_triangle_radius_jumps.push_back(
                jumps[triangle]
            );
            result.flagged_triangle_radius_jump_thresholds.push_back(
                thresholds[triangle]
            );
        }
    }
    if (options.preserve_lowest_jump_triangle_if_all_removed &&
        !keep.empty() &&
        std::none_of(keep.begin(), keep.end(), [](bool value) {
            return value;
        })) {
        const auto lowest =
            std::min_element(jumps.begin(), jumps.end());
        keep[static_cast<std::size_t>(
            std::distance(jumps.begin(), lowest)
        )] = true;
    }

    std::vector<bool> topology_repaired(
        input.triangles.size(),
        false
    );
    if (options.preserve_sheet_topology &&
        input.triangle_components.size() == input.triangles.size()) {
        const auto incidence =
            detail::triangle_edge_incidence(input.triangles);
        std::vector<std::vector<std::size_t>> rejected_adjacency(
            input.triangles.size()
        );
        for (const auto& [edge, incident] : incidence) {
            if (incident.size() != 2) {
                continue;
            }
            const std::size_t first = incident[0];
            const std::size_t second = incident[1];
            if (!keep[first] &&
                !keep[second] &&
                input.triangle_components[first] ==
                    input.triangle_components[second]) {
                rejected_adjacency[first].push_back(second);
                rejected_adjacency[second].push_back(first);
            }
        }

        const auto triangle_area = [&](std::size_t triangle) {
            const auto& face = input.triangles[triangle];
            return 0.5 * norm(cross(
                input.vertices[face[1]] - input.vertices[face[0]],
                input.vertices[face[2]] - input.vertices[face[0]]
            ));
        };
        std::vector<double> component_areas(
            input.component_count,
            0.0
        );
        for (std::size_t triangle = 0;
             triangle < input.triangles.size();
             ++triangle) {
            const std::size_t component =
                input.triangle_components[triangle];
            if (component < component_areas.size()) {
                component_areas[component] += triangle_area(triangle);
            }
        }

        std::vector<bool> visited(input.triangles.size(), false);
        for (std::size_t seed = 0;
             seed < input.triangles.size();
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
                patch.size() >
                    options.maximum_topology_repair_triangles) {
                continue;
            }

            const std::size_t source_component =
                input.triangle_components[seed];
            if (source_component >= component_areas.size() ||
                component_areas[source_component] <= 0.0) {
                continue;
            }
            double patch_area = 0.0;
            for (std::size_t triangle : patch) {
                patch_area += triangle_area(triangle);
            }
            if (patch_area / component_areas[source_component] >
                options.maximum_topology_repair_area_fraction) {
                continue;
            }

            std::set<std::size_t> patch_set(
                patch.begin(),
                patch.end()
            );
            std::set<std::size_t> boundary_neighbors;
            bool enclosed = true;
            bool touches_seam = false;
            std::size_t boundary_edge_count = 0;
            double retained_interface_length = 0.0;
            double original_boundary_length = 0.0;
            for (std::size_t triangle : patch) {
                const auto& face = input.triangles[triangle];
                for (std::size_t edge_index = 0;
                     edge_index < 3;
                     ++edge_index) {
                    const detail::ComplexEdge edge =
                        detail::normalized_complex_edge(
                            face[edge_index],
                            face[(edge_index + 1) % 3]
                        );
                    const auto found = incidence.find(edge);
                    if (found == incidence.end()) {
                        enclosed = false;
                        continue;
                    }
                    std::size_t patch_incidence = 0;
                    for (std::size_t incident_triangle :
                         found->second) {
                        if (patch_set.count(incident_triangle) != 0) {
                            ++patch_incidence;
                        }
                    }
                    if (patch_incidence != 1) {
                        continue;
                    }
                    ++boundary_edge_count;
                    const double edge_length = norm(
                        input.vertices[edge[1]] -
                        input.vertices[edge[0]]
                    );
                    if (found->second.size() == 1) {
                        original_boundary_length += edge_length;
                        enclosed = false;
                        continue;
                    }
                    if (found->second.size() > 2) {
                        touches_seam = true;
                        enclosed = false;
                    }
                    bool has_retained_neighbor = false;
                    for (std::size_t neighbor : found->second) {
                        if (patch_set.count(neighbor) != 0 ||
                            !keep[neighbor] ||
                            input.triangle_components[neighbor] !=
                                source_component) {
                            continue;
                        }
                        has_retained_neighbor = true;
                        boundary_neighbors.insert(neighbor);
                    }
                    if (has_retained_neighbor) {
                        retained_interface_length += edge_length;
                    } else {
                        enclosed = false;
                    }
                }
            }

            const bool closes_hole =
                enclosed &&
                boundary_edge_count >= 3 &&
                boundary_neighbors.size() >= 3;
            const bool avoids_notch =
                boundary_neighbors.size() >= 2 &&
                retained_interface_length >
                    original_boundary_length;
            const bool preserves_seam =
                touches_seam &&
                boundary_neighbors.size() >= 2;
            if (!closes_hole && !avoids_notch && !preserves_seam) {
                continue;
            }
            for (std::size_t triangle : patch) {
                keep[triangle] = true;
                topology_repaired[triangle] = true;
            }
            ++result.topology_repaired_patch_count;
        }
    }

    result.retained = input;
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
         triangle < input.triangles.size();
         ++triangle) {
        const double confidence =
            triangle < input.triangle_confidences.size()
                ? input.triangle_confidences[triangle]
                : 0.0;
        const double contact_angle =
            triangle < input.triangle_contact_angles.size()
                ? input.triangle_contact_angles[triangle]
                : 0.0;
        const double support =
            triangle < input.triangle_support_weights.size()
                ? input.triangle_support_weights[triangle]
                : 0.0;
        const double continuity_weight =
            jumps[triangle] <= thresholds[triangle] ||
                    jumps[triangle] <= 1e-12
                ? 1.0
                : std::clamp(
                    thresholds[triangle] / jumps[triangle],
                    0.0,
                    1.0
                );
        if (keep[triangle]) {
            result.retained.triangles.push_back(input.triangles[triangle]);
            result.retained.triangle_support_weights.push_back(support);
            result.retained.triangle_confidences.push_back(
                confidence * continuity_weight
            );
            result.retained.triangle_contact_angles.push_back(contact_angle);
            result.retained.triangle_radius_jumps.push_back(jumps[triangle]);
            result.retained.triangle_radius_jump_thresholds.push_back(
                thresholds[triangle]
            );
            if (topology_repaired[triangle]) {
                result.topology_repaired_triangles.push_back(
                    input.triangles[triangle]
                );
                result.
                    topology_repaired_triangle_radius_jumps.push_back(
                        jumps[triangle]
                    );
                result.
                    topology_repaired_triangle_radius_jump_thresholds.
                    push_back(thresholds[triangle]);
            }
        } else {
            result.removed_triangles.push_back(input.triangles[triangle]);
            result.removed_triangle_confidences.push_back(confidence);
            result.removed_triangle_contact_angles.push_back(contact_angle);
            result.removed_triangle_radius_jumps.push_back(jumps[triangle]);
            result.removed_triangle_radius_jump_thresholds.push_back(
                thresholds[triangle]
            );
        }
    }

    if (!result.retained.triangles.empty()) {
        detail::orient_and_label_sheet_components(result.retained);
        detail::classify_complex_topology(result.retained);
    }
    return result;
}

inline double medial_triangle_area(
    const MedialComplex& complex,
    const std::array<std::size_t, 3>& triangle) {
    const Vec3& a = complex.vertices[triangle[0]];
    const Vec3& b = complex.vertices[triangle[1]];
    const Vec3& c = complex.vertices[triangle[2]];
    return 0.5 * norm(cross(b - a, c - a));
}

inline double medial_edge_length(
    const MedialComplex& complex,
    const std::array<std::size_t, 2>& edge) {
    return norm(complex.vertices[edge[1]] - complex.vertices[edge[0]]);
}

struct TriangularHoleSealResult {
    MedialComplex sealed;
    std::vector<std::array<std::size_t, 3>> restored_triangles;
    std::vector<std::size_t> reference_triangle_indices;
};

inline std::array<std::size_t, 3> canonical_triangle(
    std::array<std::size_t, 3> triangle) {
    std::sort(triangle.begin(), triangle.end());
    return triangle;
}

inline TriangularHoleSealResult seal_isolated_triangular_holes(
    const MedialComplex& reference,
    const MedialComplex& input) {
    TriangularHoleSealResult result;
    result.sealed = input;
    if (reference.vertices.size() != input.vertices.size() ||
        reference.triangles.empty()) {
        return result;
    }

    std::set<std::array<std::size_t, 3>> retained_faces;
    for (const auto& triangle : input.triangles) {
        retained_faces.insert(canonical_triangle(triangle));
    }
    const auto retained_incidence =
        detail::triangle_edge_incidence(input.triangles);

    for (std::size_t triangle_index = 0;
         triangle_index < reference.triangles.size();
         ++triangle_index) {
        const auto& triangle = reference.triangles[triangle_index];
        if (retained_faces.count(canonical_triangle(triangle)) != 0) {
            continue;
        }

        bool is_isolated_triangular_hole = true;
        for (std::size_t edge_index = 0;
             edge_index < 3;
             ++edge_index) {
            const detail::ComplexEdge edge =
                detail::normalized_complex_edge(
                    triangle[edge_index],
                    triangle[(edge_index + 1) % 3]
                );
            const auto found = retained_incidence.find(edge);
            // Each side of an isolated missing triangle is a boundary edge
            // of the retained sheet. Requiring incidence one prevents this
            // pass from bridging genuine outer boundaries or seams.
            if (found == retained_incidence.end() ||
                found->second.size() != 1) {
                is_isolated_triangular_hole = false;
                break;
            }
        }
        if (!is_isolated_triangular_hole) {
            continue;
        }

        result.sealed.triangles.push_back(triangle);
        if (reference.triangle_support_weights.size() ==
                reference.triangles.size() &&
            result.sealed.triangle_support_weights.size() + 1 ==
                result.sealed.triangles.size()) {
            result.sealed.triangle_support_weights.push_back(
                reference.triangle_support_weights[triangle_index]
            );
        }
        result.sealed.triangle_confidences.push_back(
            triangle_index < reference.triangle_confidences.size()
                ? reference.triangle_confidences[triangle_index]
                : 0.0
        );
        result.sealed.triangle_contact_angles.push_back(
            triangle_index < reference.triangle_contact_angles.size()
                ? reference.triangle_contact_angles[triangle_index]
                : 0.0
        );
        result.sealed.triangle_radius_jumps.push_back(
            triangle_index < reference.triangle_radius_jumps.size()
                ? reference.triangle_radius_jumps[triangle_index]
                : detail::triangle_max_radius_jump(
                    reference,
                    triangle
                )
        );
        result.sealed.triangle_radius_jump_thresholds.push_back(
            triangle_index <
                    reference.triangle_radius_jump_thresholds.size()
                ? reference.
                    triangle_radius_jump_thresholds[triangle_index]
                : 0.0
        );
        result.restored_triangles.push_back(triangle);
        result.reference_triangle_indices.push_back(triangle_index);
    }

    if (!result.restored_triangles.empty()) {
        result.sealed.triangle_components.clear();
        result.sealed.boundary_edges.clear();
        result.sealed.seam_edges.clear();
        result.sealed.junction_vertices.clear();
        result.sealed.termination_edges.clear();
        result.sealed.artificial_boundary_edges.clear();
        result.sealed.boundary_loops.clear();
        result.sealed.component_count = 0;
        detail::orient_and_label_sheet_components(result.sealed);
        detail::classify_complex_topology(result.sealed);
    }
    return result;
}

inline std::vector<MedialComponentMetrics> analyze_medial_components(
    const MedialComplex& complex,
    const PoleSelectionResult& pole_selection) {
    if (complex.triangle_components.size() != complex.triangles.size() ||
        complex.component_count == 0) {
        return {};
    }

    std::vector<MedialComponentMetrics> metrics(complex.component_count);
    std::vector<std::set<std::size_t>> component_vertices(
        complex.component_count
    );
    std::vector<std::set<std::size_t>> component_supporting_poles(
        complex.component_count
    );
    double total_area = 0.0;

    for (std::size_t component = 0;
         component < metrics.size();
         ++component) {
        metrics[component].component = component;
        metrics[component].minimum_confidence =
            std::numeric_limits<double>::infinity();
    }

    for (std::size_t triangle_index = 0;
         triangle_index < complex.triangles.size();
         ++triangle_index) {
        const std::size_t component =
            complex.triangle_components[triangle_index];
        if (component >= metrics.size()) {
            continue;
        }
        MedialComponentMetrics& current = metrics[component];
        ++current.triangle_count;
        const double area =
            medial_triangle_area(complex, complex.triangles[triangle_index]);
        current.area += area;
        total_area += area;

        const double confidence =
            triangle_index < complex.triangle_confidences.size()
                ? complex.triangle_confidences[triangle_index]
                : 0.0;
        current.mean_confidence += confidence;
        current.minimum_confidence =
            std::min(current.minimum_confidence, confidence);
        const double radius_jump =
            triangle_index < complex.triangle_radius_jumps.size()
                ? complex.triangle_radius_jumps[triangle_index]
                : detail::triangle_max_radius_jump(
                    complex,
                    complex.triangles[triangle_index]
                );
        current.mean_radius_jump += radius_jump;
        current.maximum_radius_jump =
            std::max(current.maximum_radius_jump, radius_jump);

        for (std::size_t vertex : complex.triangles[triangle_index]) {
            component_vertices[component].insert(vertex);
        }
    }

    for (std::size_t component = 0;
         component < metrics.size();
         ++component) {
        for (std::size_t vertex : component_vertices[component]) {
            if (vertex >= complex.vertex_source_tetrahedra.size()) {
                continue;
            }
            const auto& sources =
                complex.vertex_source_tetrahedra[vertex];
            for (std::size_t pole_index = 0;
                 pole_index < pole_selection.poles.size();
                 ++pole_index) {
                const MedialCandidate& pole =
                    pole_selection.poles[pole_index];
                if (pole.validated &&
                    std::find(
                        sources.begin(),
                        sources.end(),
                        pole.source_tetrahedron) != sources.end()) {
                    component_supporting_poles[component].insert(pole_index);
                }
            }
        }
    }

    const auto incidence =
        detail::triangle_edge_incidence(complex.triangles);
    for (const auto& [edge, incident] : incidence) {
        const double length = medial_edge_length(complex, edge);
        if (incident.size() == 1) {
            const std::size_t component =
                complex.triangle_components[incident.front()];
            if (component < metrics.size()) {
                metrics[component].boundary_length += length;
            }
        } else if (incident.size() > 2) {
            std::set<std::size_t> incident_components;
            for (std::size_t triangle : incident) {
                incident_components.insert(
                    complex.triangle_components[triangle]
                );
            }
            for (std::size_t component : incident_components) {
                if (component < metrics.size()) {
                    metrics[component].seam_length += length;
                }
            }
        }
    }

    for (const MedialBoundaryLoop& loop : complex.boundary_loops) {
        if (loop.component >= metrics.size()) {
            continue;
        }
        ++metrics[loop.component].boundary_loop_count;
        if (!loop.allowed) {
            ++metrics[loop.component].artificial_boundary_loop_count;
            metrics[loop.component].artificial_boundary_length +=
                loop.length;
        }
    }

    for (std::size_t component = 0;
         component < metrics.size();
         ++component) {
        MedialComponentMetrics& current = metrics[component];
        current.vertex_count = component_vertices[component].size();
        current.supporting_pole_count =
            component_supporting_poles[component].size();
        current.area_fraction =
            total_area > 0.0 ? current.area / total_area : 0.0;
        current.artificial_boundary_fraction =
            current.boundary_length > 1e-12
                ? current.artificial_boundary_length /
                      current.boundary_length
                : 0.0;
        if (current.triangle_count > 0) {
            current.mean_confidence /=
                static_cast<double>(current.triangle_count);
            current.mean_radius_jump /=
                static_cast<double>(current.triangle_count);
        } else {
            current.minimum_confidence = 0.0;
        }
    }
    return metrics;
}

inline MedialComplexFilterResult filter_medial_complex_components(
    const MedialComplex& input,
    const PoleSelectionResult& pole_selection,
    const MedialComponentFilterOptions& options = {}) {
    MedialComplexFilterResult result;
    MedialComplex working = input;
    if (!working.triangles.empty() &&
        working.triangle_components.size() != working.triangles.size()) {
        detail::orient_and_label_sheet_components(working);
    }

    result.source_component_metrics =
        analyze_medial_components(working, pole_selection);
    std::vector<bool> keep(result.source_component_metrics.size(), false);
    const auto has_valid_boundaries =
        [&](const MedialComponentMetrics& component) {
            return !options.require_valid_sheet_boundaries ||
                   component.artificial_boundary_loop_count == 0 ||
                   component.artificial_boundary_fraction <=
                       options.maximum_artificial_boundary_fraction;
        };
    for (MedialComponentMetrics& component :
         result.source_component_metrics) {
        component.retained =
            component.triangle_count >= options.minimum_triangle_count &&
            component.area_fraction >= options.minimum_area_fraction &&
            component.mean_confidence >=
                options.minimum_mean_confidence &&
            component.supporting_pole_count >=
                options.minimum_supporting_poles &&
            has_valid_boundaries(component);
        keep[component.component] = component.retained;
    }

    if (!working.triangles.empty() &&
        options.preserve_largest_component_if_all_removed &&
        std::none_of(keep.begin(), keep.end(), [](bool value) {
            return value;
        }) &&
        !result.source_component_metrics.empty()) {
        const auto largest = std::max_element(
            result.source_component_metrics.begin(),
            result.source_component_metrics.end(),
            [&](const MedialComponentMetrics& first,
                const MedialComponentMetrics& second) {
                const bool first_valid =
                    has_valid_boundaries(first);
                const bool second_valid =
                    has_valid_boundaries(second);
                if (first_valid != second_valid) {
                    return !first_valid;
                }
                return first.area < second.area;
            }
        );
        if (has_valid_boundaries(*largest)) {
            largest->retained = true;
            keep[largest->component] = true;
        }
    }

    result.retained = working;
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
         triangle < working.triangles.size();
         ++triangle) {
        const std::size_t component =
            working.triangle_components[triangle];
        const bool retained =
            component < keep.size() && keep[component];
        const double confidence =
            triangle < working.triangle_confidences.size()
                ? working.triangle_confidences[triangle]
                : 0.0;
        const double support =
            triangle < working.triangle_support_weights.size()
                ? working.triangle_support_weights[triangle]
                : 0.0;
        const double contact_angle =
            triangle < working.triangle_contact_angles.size()
                ? working.triangle_contact_angles[triangle]
                : 0.0;
        const double radius_jump =
            triangle < working.triangle_radius_jumps.size()
                ? working.triangle_radius_jumps[triangle]
                : detail::triangle_max_radius_jump(
                    working,
                    working.triangles[triangle]
                );
        const double radius_jump_threshold =
            triangle < working.triangle_radius_jump_thresholds.size()
                ? working.triangle_radius_jump_thresholds[triangle]
                : 0.0;
        if (retained) {
            result.retained.triangles.push_back(
                working.triangles[triangle]
            );
            result.retained.triangle_support_weights.push_back(support);
            result.retained.triangle_confidences.push_back(confidence);
            result.retained.triangle_contact_angles.push_back(contact_angle);
            result.retained.triangle_radius_jumps.push_back(radius_jump);
            result.retained.triangle_radius_jump_thresholds.push_back(
                radius_jump_threshold
            );
        } else {
            result.removed_triangles.push_back(working.triangles[triangle]);
            result.removed_triangle_confidences.push_back(confidence);
            result.removed_triangle_contact_angles.push_back(contact_angle);
            result.removed_triangle_radius_jumps.push_back(radius_jump);
            result.removed_triangle_radius_jump_thresholds.push_back(
                radius_jump_threshold
            );
            result.removed_triangle_source_components.push_back(component);
        }
    }

    if (!result.retained.triangles.empty()) {
        detail::orient_and_label_sheet_components(result.retained);
        detail::classify_complex_topology(result.retained);
        result.retained_component_metrics =
            analyze_medial_components(result.retained, pole_selection);
        for (MedialComponentMetrics& component :
             result.retained_component_metrics) {
            component.retained = true;
        }
    }
    return result;
}

}  // namespace medial_axis_3d
