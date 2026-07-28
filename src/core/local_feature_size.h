#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

#include "medial_candidate.h"
#include "vec3.h"

namespace medial_axis_3d {

struct LocalFeatureSizeOptions {
    std::size_t resolution_neighbor_count{6};
    std::size_t interpolation_anchor_count{4};
    // This is only a numerical floor. Values h/LFS > 1 are intentionally
    // preserved because they identify locally undersampled surface regions.
    double minimum_lfs_resolution_ratio{0.25};
    double maximum_lfs_diagonal_ratio{2.0};
    double minimum_normal_angle_radians{1e-3};
};

struct SurfaceFeatureField {
    std::vector<double> resolutions;
    std::vector<double> local_feature_sizes;
    std::vector<double> sampling_densities;
    std::vector<bool> pole_anchored;
};

struct AdaptiveSamplingWeightOptions {
    // Squaring h/LFS gives meaningfully more budget to locally undersampled
    // regions while keeping the redistribution bounded.
    double density_exponent{2.0};
    double minimum_triangle_importance{0.25};
    double maximum_triangle_importance{4.0};
};

namespace detail {

inline double median_value(std::vector<double> values) {
    if (values.empty()) {
        return 0.0;
    }
    std::sort(values.begin(), values.end());
    const std::size_t middle = values.size() / 2;
    return values.size() % 2 == 0
        ? 0.5 * (values[middle - 1] + values[middle])
        : values[middle];
}

inline double point_set_diagonal(const std::vector<Vec3>& points) {
    if (points.empty()) {
        return 0.0;
    }
    Vec3 lower = points.front();
    Vec3 upper = points.front();
    for (const Vec3& point : points) {
        lower.x = std::min(lower.x, point.x);
        lower.y = std::min(lower.y, point.y);
        lower.z = std::min(lower.z, point.z);
        upper.x = std::max(upper.x, point.x);
        upper.y = std::max(upper.y, point.y);
        upper.z = std::max(upper.z, point.z);
    }
    return norm(upper - lower);
}

inline void assign_lfs_anchor(std::vector<double>& anchors,
                              std::vector<bool>& anchored,
                              std::size_t sample,
                              double radius) {
    if (sample >= anchors.size() || !std::isfinite(radius) ||
        radius <= 0.0) {
        return;
    }
    anchors[sample] = anchored[sample]
        ? std::min(anchors[sample], radius)
        : radius;
    anchored[sample] = true;
}

}  // namespace detail

inline SurfaceFeatureField estimate_surface_feature_field(
    const std::vector<Vec3>& points,
    const std::vector<Vec3>& normals,
    const PoleSelectionResult& pole_selection,
    const LocalFeatureSizeOptions& options = {}) {
    SurfaceFeatureField result;
    const std::size_t sample_count = points.size();
    result.resolutions.assign(sample_count, 0.0);
    result.local_feature_sizes.assign(sample_count, 0.0);
    result.sampling_densities.assign(sample_count, 0.0);
    result.pole_anchored.assign(sample_count, false);
    if (sample_count == 0) {
        return result;
    }

    const double diagonal =
        std::max(1e-12, detail::point_set_diagonal(points));
    const double maximum_lfs =
        diagonal * std::max(1.0, options.maximum_lfs_diagonal_ratio);
    const std::size_t neighbor_count = std::min(
        options.resolution_neighbor_count,
        sample_count > 0 ? sample_count - 1 : 0
    );

    using Neighbor = std::pair<double, std::size_t>;
    std::vector<std::vector<Neighbor>> neighborhoods(sample_count);
    std::vector<double> curvature_fallback(sample_count, maximum_lfs);
    for (std::size_t sample = 0; sample < sample_count; ++sample) {
        std::vector<Neighbor> neighbors;
        neighbors.reserve(sample_count > 0 ? sample_count - 1 : 0);
        for (std::size_t candidate = 0;
             candidate < sample_count;
             ++candidate) {
            if (candidate == sample) {
                continue;
            }
            const double distance = norm(points[candidate] - points[sample]);
            if (distance > 1e-14) {
                neighbors.push_back({distance, candidate});
            }
        }
        std::sort(
            neighbors.begin(),
            neighbors.end(),
            [](const Neighbor& first, const Neighbor& second) {
                return first.first < second.first;
            }
        );
        if (neighbors.size() > neighbor_count) {
            neighbors.resize(neighbor_count);
        }
        neighborhoods[sample] = neighbors;

        std::vector<double> distances;
        distances.reserve(neighbors.size());
        for (const Neighbor& neighbor : neighbors) {
            distances.push_back(neighbor.first);
        }
        double resolution = detail::median_value(std::move(distances));
        if (resolution <= 0.0) {
            resolution = diagonal /
                std::max(1.0, std::sqrt(static_cast<double>(sample_count)));
        }
        result.resolutions[sample] = resolution;

        if (normals.size() == sample_count &&
            norm(normals[sample]) > 0.0) {
            std::vector<double> curvature_radii;
            for (const Neighbor& neighbor : neighbors) {
                if (norm(normals[neighbor.second]) <= 0.0) {
                    continue;
                }
                const double cosine = std::clamp(
                    dot(
                        normalized(normals[sample]),
                        normalized(normals[neighbor.second])
                    ),
                    -1.0,
                    1.0
                );
                const double angle = std::acos(cosine);
                if (angle >= options.minimum_normal_angle_radians) {
                    curvature_radii.push_back(
                        neighbor.first / angle
                    );
                }
            }
            if (!curvature_radii.empty()) {
                curvature_fallback[sample] = std::clamp(
                    detail::median_value(std::move(curvature_radii)),
                    options.minimum_lfs_resolution_ratio * resolution,
                    maximum_lfs
                );
            }
        }
    }

    std::vector<double> anchors(sample_count, 0.0);
    for (const MedialCandidate& pole : pole_selection.poles) {
        if (!pole.validated || pole.radius <= 0.0) {
            continue;
        }
        for (std::size_t sample : pole.source_samples) {
            detail::assign_lfs_anchor(
                anchors,
                result.pole_anchored,
                sample,
                pole.radius
            );
        }
        for (const SurfaceContact& contact : pole.contacts) {
            std::size_t closest_sample = 0;
            double closest_distance =
                std::numeric_limits<double>::infinity();
            for (std::size_t sample = 0;
                 sample < sample_count;
                 ++sample) {
                const double distance =
                    squared_norm(points[sample] - contact.position);
                if (distance < closest_distance) {
                    closest_distance = distance;
                    closest_sample = sample;
                }
            }
            detail::assign_lfs_anchor(
                anchors,
                result.pole_anchored,
                closest_sample,
                pole.radius
            );
        }
    }

    std::vector<std::size_t> anchor_indices;
    for (std::size_t sample = 0; sample < sample_count; ++sample) {
        if (result.pole_anchored[sample]) {
            anchor_indices.push_back(sample);
        }
    }

    for (std::size_t sample = 0; sample < sample_count; ++sample) {
        double pole_estimate = maximum_lfs;
        if (result.pole_anchored[sample]) {
            pole_estimate = anchors[sample];
        } else if (!anchor_indices.empty()) {
            std::vector<Neighbor> nearest_anchors;
            nearest_anchors.reserve(anchor_indices.size());
            for (std::size_t anchor : anchor_indices) {
                nearest_anchors.push_back({
                    norm(points[sample] - points[anchor]),
                    anchor
                });
            }
            std::sort(
                nearest_anchors.begin(),
                nearest_anchors.end(),
                [](const Neighbor& first, const Neighbor& second) {
                    return first.first < second.first;
                }
            );
            if (nearest_anchors.size() >
                options.interpolation_anchor_count) {
                nearest_anchors.resize(
                    options.interpolation_anchor_count
                );
            }
            double weighted_radius = 0.0;
            double total_weight = 0.0;
            for (const Neighbor& anchor : nearest_anchors) {
                const double weight =
                    1.0 / std::max(anchor.first, 1e-12);
                weighted_radius += weight * anchors[anchor.second];
                total_weight += weight;
            }
            if (total_weight > 0.0) {
                pole_estimate = weighted_radius / total_weight;
            }
        }

        const double lower_bound =
            options.minimum_lfs_resolution_ratio *
            result.resolutions[sample];
        result.local_feature_sizes[sample] = std::clamp(
            std::min(pole_estimate, curvature_fallback[sample]),
            lower_bound,
            maximum_lfs
        );
        result.sampling_densities[sample] =
            result.resolutions[sample] /
            std::max(result.local_feature_sizes[sample], 1e-12);
    }
    return result;
}

inline std::vector<double> lfs_adaptive_triangle_importance(
    const Mesh& mesh,
    const SurfaceFeatureField& feature_field,
    const AdaptiveSamplingWeightOptions& options = {}) {
    std::vector<double> importance(mesh.faces.size(), 1.0);
    if (feature_field.sampling_densities.size() != mesh.vertices.size() ||
        !std::isfinite(options.density_exponent) ||
        options.density_exponent < 0.0 ||
        !std::isfinite(options.minimum_triangle_importance) ||
        !std::isfinite(options.maximum_triangle_importance) ||
        options.minimum_triangle_importance < 0.0 ||
        options.maximum_triangle_importance <
            options.minimum_triangle_importance) {
        return importance;
    }

    for (std::size_t face_index = 0;
         face_index < mesh.faces.size();
         ++face_index) {
        const auto& face = mesh.faces[face_index].vertices;
        double maximum_density = 0.0;
        bool valid_face = true;
        for (int vertex : face) {
            if (vertex < 0 ||
                static_cast<std::size_t>(vertex) >=
                    feature_field.sampling_densities.size()) {
                valid_face = false;
                break;
            }
            const double density =
                feature_field.sampling_densities[
                    static_cast<std::size_t>(vertex)
                ];
            if (!std::isfinite(density) || density < 0.0) {
                valid_face = false;
                break;
            }
            maximum_density = std::max(maximum_density, density);
        }
        if (!valid_face) {
            continue;
        }
        const double weighted_density =
            std::pow(maximum_density, options.density_exponent);
        importance[face_index] = std::clamp(
            weighted_density,
            options.minimum_triangle_importance,
            options.maximum_triangle_importance
        );
    }
    return importance;
}

}  // namespace medial_axis_3d
