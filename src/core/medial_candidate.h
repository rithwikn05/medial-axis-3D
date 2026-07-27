#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "delaunay3.h"
#include "mesh.h"
#include "surface_sample.h"
#include "vec3.h"

namespace medial_axis_3d {

struct SurfaceContact {
    Vec3 position{};
    Vec3 normal{};
    std::size_t triangle{0};
    double distance{0.0};
};

struct MedialCandidate {
    Vec3 center{};
    double radius{0.0};
    double voronoi_radius{0.0};
    double radius_agreement{0.0};
    double maximum_contact_angle_degrees{0.0};
    double confidence{0.0};
    std::size_t source_tetrahedron{0};
    std::vector<std::size_t> source_samples;
    std::vector<SurfaceContact> contacts;
    bool validated{false};
};

struct PoleSelectionOptions {
    double minimum_inward_cosine{0.05};
    double contact_relative_tolerance{0.08};
    double minimum_contact_angle_degrees{90.0};
    double minimum_radius_agreement{0.65};
    double minimum_radius{1e-8};
};

struct PoleSelectionResult {
    std::vector<MedialCandidate> poles;

    std::size_t validated_count() const {
        return static_cast<std::size_t>(std::count_if(
            poles.begin(),
            poles.end(),
            [](const MedialCandidate& pole) {
                return pole.validated;
            }
        ));
    }
};

inline Vec3 closest_point_on_triangle(const Vec3& point,
                                      const Vec3& a,
                                      const Vec3& b,
                                      const Vec3& c) {
    const Vec3 ab = b - a;
    const Vec3 ac = c - a;
    const Vec3 ap = point - a;
    const double d1 = dot(ab, ap);
    const double d2 = dot(ac, ap);
    if (d1 <= 0.0 && d2 <= 0.0) {
        return a;
    }

    const Vec3 bp = point - b;
    const double d3 = dot(ab, bp);
    const double d4 = dot(ac, bp);
    if (d3 >= 0.0 && d4 <= d3) {
        return b;
    }

    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        const double v = d1 / (d1 - d3);
        return a + ab * v;
    }

    const Vec3 cp = point - c;
    const double d5 = dot(ab, cp);
    const double d6 = dot(ac, cp);
    if (d6 >= 0.0 && d5 <= d6) {
        return c;
    }

    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        const double w = d2 / (d2 - d6);
        return a + ac * w;
    }

    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        const Vec3 bc = c - b;
        const double w = (d4 - d3) /
                         ((d4 - d3) + (d5 - d6));
        return b + bc * w;
    }

    const double denominator = 1.0 / (va + vb + vc);
    const double v = vb * denominator;
    const double w = vc * denominator;
    return a + ab * v + ac * w;
}

inline std::vector<SurfaceContact> nearest_surface_contacts(
    const Mesh& mesh,
    const Vec3& center,
    double relative_tolerance) {
    struct ContactCandidate {
        SurfaceContact contact;
        double squared_distance{0.0};
    };

    std::vector<ContactCandidate> candidates;
    candidates.reserve(mesh.faces.size());
    double minimum_squared_distance = std::numeric_limits<double>::infinity();

    for (std::size_t face_index = 0; face_index < mesh.faces.size(); ++face_index) {
        const auto& face = mesh.faces[face_index].vertices;
        const Vec3& a = mesh.vertices[static_cast<std::size_t>(face[0])];
        const Vec3& b = mesh.vertices[static_cast<std::size_t>(face[1])];
        const Vec3& c = mesh.vertices[static_cast<std::size_t>(face[2])];
        const Vec3 closest = closest_point_on_triangle(center, a, b, c);
        const double distance_squared = squared_norm(closest - center);

        ContactCandidate candidate;
        candidate.contact.position = closest;
        candidate.contact.normal = mesh.face_normals[face_index];
        candidate.contact.triangle = face_index;
        candidate.contact.distance = std::sqrt(distance_squared);
        candidate.squared_distance = distance_squared;
        candidates.push_back(candidate);
        minimum_squared_distance =
            std::min(minimum_squared_distance, distance_squared);
    }

    if (!std::isfinite(minimum_squared_distance)) {
        return {};
    }

    const double minimum_distance = std::sqrt(minimum_squared_distance);
    const double maximum_distance =
        minimum_distance * (1.0 + std::max(0.0, relative_tolerance)) + 1e-12;
    const double duplicate_tolerance_squared =
        std::max(1e-24, minimum_squared_distance * 1e-12);

    std::vector<SurfaceContact> result;
    for (const ContactCandidate& candidate : candidates) {
        if (candidate.contact.distance > maximum_distance) {
            continue;
        }
        const bool duplicate = std::any_of(
            result.begin(),
            result.end(),
            [&](const SurfaceContact& existing) {
                return squared_norm(
                    existing.position - candidate.contact.position
                ) <= duplicate_tolerance_squared;
            }
        );
        if (!duplicate) {
            result.push_back(candidate.contact);
        }
    }
    return result;
}

inline double maximum_contact_angle_degrees(
    const Vec3& center,
    const std::vector<SurfaceContact>& contacts) {
    double maximum_angle = 0.0;
    constexpr double radians_to_degrees = 57.2957795130823208768;
    for (std::size_t first = 0; first < contacts.size(); ++first) {
        const Vec3 first_direction =
            normalized(contacts[first].position - center);
        for (std::size_t second = first + 1; second < contacts.size(); ++second) {
            const Vec3 second_direction =
                normalized(contacts[second].position - center);
            const double cosine = std::clamp(
                dot(first_direction, second_direction),
                -1.0,
                1.0
            );
            maximum_angle = std::max(
                maximum_angle,
                std::acos(cosine) * radians_to_degrees
            );
        }
    }
    return maximum_angle;
}

inline PoleSelectionResult select_inward_poles(
    const Delaunay3& delaunay,
    const std::vector<SurfaceSample>& samples,
    const Mesh& mesh,
    const PoleSelectionOptions& options = {}) {
    PoleSelectionResult result;
    if (samples.size() != delaunay.point_count()) {
        return result;
    }

    std::vector<Vec3> circumcenters(delaunay.tetrahedron_count());
    std::vector<bool> is_interior(delaunay.tetrahedron_count(), false);
    for (std::size_t tetrahedron_index = 0;
         tetrahedron_index < delaunay.tetrahedron_count();
         ++tetrahedron_index) {
        Vec3 center{};
        if (delaunay.circumcenter(tetrahedron_index, center) &&
            point_inside_mesh(mesh, center)) {
            circumcenters[tetrahedron_index] = center;
            is_interior[tetrahedron_index] = true;
        }
    }

    for (std::size_t sample_index = 0;
         sample_index < samples.size();
         ++sample_index) {
        const SurfaceSample& sample = samples[sample_index];
        const double normal_length = norm(sample.normal);
        if (normal_length == 0.0) {
            continue;
        }

        std::size_t best_tetrahedron = delaunay.tetrahedron_count();
        double best_distance_squared = -1.0;
        for (std::size_t tetrahedron_index = 0;
             tetrahedron_index < delaunay.tetrahedron_count();
             ++tetrahedron_index) {
            if (!is_interior[tetrahedron_index]) {
                continue;
            }
            const Tetrahedron& tetrahedron =
                delaunay.tetrahedra()[tetrahedron_index];
            if (std::find(
                    tetrahedron.vertices.begin(),
                    tetrahedron.vertices.end(),
                    static_cast<int>(sample_index)
                ) == tetrahedron.vertices.end()) {
                continue;
            }

            const Vec3 offset =
                circumcenters[tetrahedron_index] - sample.position;
            const double distance_squared = squared_norm(offset);
            if (distance_squared <= 0.0) {
                continue;
            }
            const double inward_cosine =
                -dot(normalized(offset), sample.normal / normal_length);
            if (inward_cosine < options.minimum_inward_cosine) {
                continue;
            }
            if (distance_squared > best_distance_squared) {
                best_distance_squared = distance_squared;
                best_tetrahedron = tetrahedron_index;
            }
        }

        if (best_tetrahedron == delaunay.tetrahedron_count()) {
            continue;
        }

        const Vec3 center = circumcenters[best_tetrahedron];
        auto existing = std::find_if(
            result.poles.begin(),
            result.poles.end(),
            [&](const MedialCandidate& pole) {
                return squared_norm(pole.center - center) < 1e-20;
            }
        );
        if (existing != result.poles.end()) {
            existing->source_samples.push_back(sample_index);
            continue;
        }

        MedialCandidate pole;
        pole.center = center;
        pole.voronoi_radius = std::sqrt(best_distance_squared);
        pole.source_tetrahedron = best_tetrahedron;
        pole.source_samples.push_back(sample_index);
        pole.contacts = nearest_surface_contacts(
            mesh,
            center,
            options.contact_relative_tolerance
        );
        if (!pole.contacts.empty()) {
            pole.radius = std::min_element(
                pole.contacts.begin(),
                pole.contacts.end(),
                [](const SurfaceContact& a, const SurfaceContact& b) {
                    return a.distance < b.distance;
                }
            )->distance;
        }
        pole.radius_agreement = pole.voronoi_radius > 0.0
            ? pole.radius / pole.voronoi_radius
            : 0.0;
        pole.maximum_contact_angle_degrees =
            maximum_contact_angle_degrees(center, pole.contacts);

        const double angle_confidence = std::clamp(
            pole.maximum_contact_angle_degrees / 180.0,
            0.0,
            1.0
        );
        const double radius_confidence =
            std::clamp(pole.radius_agreement, 0.0, 1.0);
        const double contact_confidence =
            std::min(1.0, static_cast<double>(pole.contacts.size()) / 3.0);
        pole.confidence =
            angle_confidence * radius_confidence * contact_confidence;
        pole.validated =
            pole.radius >= options.minimum_radius &&
            pole.contacts.size() >= 2 &&
            pole.maximum_contact_angle_degrees >=
                options.minimum_contact_angle_degrees &&
            pole.radius_agreement >= options.minimum_radius_agreement;
        result.poles.push_back(std::move(pole));
    }

    return result;
}

}  // namespace medial_axis_3d
