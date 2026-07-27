#include "../core/delaunay3.h"
#include "../core/cross_resolution_stability.h"
#include "../core/filters.h"
#include "../core/medial_axis_approx.h"
#include "../core/medial_candidate.h"
#include "../core/medial_complex.h"
#include "../core/surface_sample.h"
#include "../io/tetgen_io.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <polyscope/curve_network.h>
#include <polyscope/point_cloud.h>
#include <polyscope/point_cloud_scalar_quantity.h>
#include <polyscope/point_cloud_vector_quantity.h>
#include <polyscope/polyscope.h>
#include <polyscope/screenshot.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/surface_vector_quantity.h>

namespace {

using namespace medial_axis_3d;

std::vector<glm::vec3> to_glm_points(const std::vector<Vec3>& points) {
    std::vector<glm::vec3> result;
    result.reserve(points.size());
    for (const Vec3& point : points) {
        result.emplace_back(
            static_cast<float>(point.x),
            static_cast<float>(point.y),
            static_cast<float>(point.z)
        );
    }
    return result;
}

std::vector<std::array<std::size_t, 3>> to_surface_faces(const Mesh& mesh) {
    std::vector<std::array<std::size_t, 3>> result;
    result.reserve(mesh.faces.size());
    for (const Triangle& triangle : mesh.faces) {
        result.push_back({
            static_cast<std::size_t>(triangle.vertices[0]),
            static_cast<std::size_t>(triangle.vertices[1]),
            static_cast<std::size_t>(triangle.vertices[2])
        });
    }
    return result;
}

bool samples_match_points(const std::vector<SurfaceSample>& samples,
                          const std::vector<Vec3>& points) {
    if (samples.size() != points.size()) {
        return false;
    }
    for (std::size_t i = 0; i < samples.size(); ++i) {
        if (!(samples[i].position == points[i])) {
            return false;
        }
    }
    return true;
}

bool build_resolution_stability_run(
    const Mesh& surface,
    std::size_t target_sample_count,
    const MedialComplexOptions& medial_options,
    const RadiusContinuityFilterOptions& radius_options,
    const MedialComponentFilterOptions& component_options,
    ResolutionMedialComplex& result,
    std::string& error) {
    SurfaceSamplingOptions sampling_options;
    sampling_options.target_sample_count = target_sample_count;
    // Stability resolutions must actually differ even when the input mesh
    // has more vertices than the requested count.
    sampling_options.include_mesh_vertices = false;
    const std::vector<SurfaceSample> samples =
        sample_surface(surface, sampling_options);
    if (samples.size() < 4) {
        error = "Stability resampling produced fewer than four points.";
        return false;
    }

    std::vector<Vec3> points;
    std::vector<Vec3> normals;
    points.reserve(samples.size());
    normals.reserve(samples.size());
    for (const SurfaceSample& sample : samples) {
        points.push_back(sample.position);
        normals.push_back(sample.normal);
    }

    Delaunay3 delaunay;
    if (!delaunay.build(points)) {
        error = "Could not tetrahedralize a stability-resolution sample set.";
        return false;
    }
    const PoleSelectionResult poles =
        select_inward_poles(delaunay, samples, surface);
    const SurfaceFeatureField features =
        estimate_surface_feature_field(points, normals, poles);
    const MedialComplex supported = build_validated_medial_complex(
        delaunay,
        &surface,
        poles,
        medial_options,
        &features
    );
    const RadiusContinuityFilterResult radius_filtered =
        filter_medial_complex_radius_continuity(
            supported,
            radius_options
        );
    const MedialComplexFilterResult component_filtered =
        filter_medial_complex_components(
            radius_filtered.retained,
            poles,
            component_options
        );
    result.sample_count = samples.size();
    result.complex = component_filtered.retained;
    result.component_metrics =
        component_filtered.retained_component_metrics;
    return true;
}

void print_usage(const char* program) {
    std::cerr << "Usage: " << program
              << " <input.node> [output.ele] [--samples N] "
                 "[--no-gui | --screenshot <image.png>]\n"
              << "       [--min-component-triangles N] "
                 "[--min-component-area-fraction X]\n"
              << "       [--min-component-confidence X] "
                 "[--min-component-poles N]\n"
              << "       [--max-relative-radius-jump X]\n"
              << "       [--fixed-radius-jump]\n"
              << "       [--support-rings N]\n"
              << "       [--no-cross-resolution]\n"
              << "If output.ele is supplied, the generated tetrahedra are written in TetGen format.\n"
              << "--samples N generates N deterministic area-weighted surface samples.\n"
              << "--no-gui computes and prints geometry statistics without opening Polyscope.\n"
              << "--screenshot renders one image and exits.\n";
}

}  // namespace

int main(int argc, char** argv) {
    bool no_gui = false;
    bool cross_resolution_enabled = true;
    std::size_t target_sample_count = 0;
    RadiusContinuityFilterOptions radius_filter_options;
    MedialComponentFilterOptions component_filter_options;
    MedialComplexOptions medial_complex_options;
    std::filesystem::path screenshot_path;
    std::vector<std::string> positional_arguments;
    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];
        if (argument == "--no-gui") {
            no_gui = true;
        } else if (argument == "--no-cross-resolution") {
            cross_resolution_enabled = false;
        } else if (argument == "--fixed-radius-jump") {
            radius_filter_options.use_adaptive_lfs_thresholds = false;
        } else if (argument == "--support-rings") {
            if (i + 1 >= argc) {
                print_usage(argv[0]);
                return 2;
            }
            const std::string value = argv[++i];
            try {
                std::size_t parsed = 0;
                medial_complex_options.support_ring_count =
                    static_cast<std::size_t>(
                        std::stoull(value, &parsed)
                    );
                if (parsed != value.size()) {
                    throw std::invalid_argument("trailing characters");
                }
            } catch (const std::exception&) {
                std::cerr << argument
                          << " requires a nonnegative integer; zero uses "
                             "automatic scaling.\n";
                return 2;
            }
        } else if (argument == "--samples") {
            if (i + 1 >= argc) {
                print_usage(argv[0]);
                return 2;
            }
            try {
                target_sample_count =
                    static_cast<std::size_t>(std::stoull(argv[++i]));
            } catch (const std::exception&) {
                std::cerr << "--samples requires a positive integer.\n";
                return 2;
            }
            if (target_sample_count < 4) {
                std::cerr << "--samples must be at least 4.\n";
                return 2;
            }
        } else if (argument == "--screenshot") {
            if (i + 1 >= argc) {
                print_usage(argv[0]);
                return 2;
            }
            screenshot_path = argv[++i];
        } else if (argument == "--min-component-triangles") {
            if (i + 1 >= argc) {
                print_usage(argv[0]);
                return 2;
            }
            const std::string value = argv[++i];
            if (value.empty() || value.front() == '-') {
                std::cerr << argument << " requires a nonnegative integer.\n";
                return 2;
            }
            try {
                std::size_t parsed = 0;
                component_filter_options.minimum_triangle_count =
                    static_cast<std::size_t>(std::stoull(value, &parsed));
                if (parsed != value.size()) {
                    throw std::invalid_argument("trailing characters");
                }
            } catch (const std::exception&) {
                std::cerr << argument << " requires a nonnegative integer.\n";
                return 2;
            }
        } else if (argument == "--min-component-area-fraction") {
            if (i + 1 >= argc) {
                print_usage(argv[0]);
                return 2;
            }
            const std::string value = argv[++i];
            try {
                std::size_t parsed = 0;
                component_filter_options.minimum_area_fraction =
                    std::stod(value, &parsed);
                if (parsed != value.size()) {
                    throw std::invalid_argument("trailing characters");
                }
            } catch (const std::exception&) {
                std::cerr << argument << " requires a number from 0 to 1.\n";
                return 2;
            }
            if (!std::isfinite(
                    component_filter_options.minimum_area_fraction) ||
                component_filter_options.minimum_area_fraction < 0.0 ||
                component_filter_options.minimum_area_fraction > 1.0) {
                std::cerr << argument << " must be between 0 and 1.\n";
                return 2;
            }
        } else if (argument == "--min-component-confidence") {
            if (i + 1 >= argc) {
                print_usage(argv[0]);
                return 2;
            }
            const std::string value = argv[++i];
            try {
                std::size_t parsed = 0;
                component_filter_options.minimum_mean_confidence =
                    std::stod(value, &parsed);
                if (parsed != value.size()) {
                    throw std::invalid_argument("trailing characters");
                }
            } catch (const std::exception&) {
                std::cerr << argument << " requires a number from 0 to 1.\n";
                return 2;
            }
            if (!std::isfinite(
                    component_filter_options.minimum_mean_confidence) ||
                component_filter_options.minimum_mean_confidence < 0.0 ||
                component_filter_options.minimum_mean_confidence > 1.0) {
                std::cerr << argument << " must be between 0 and 1.\n";
                return 2;
            }
        } else if (argument == "--min-component-poles") {
            if (i + 1 >= argc) {
                print_usage(argv[0]);
                return 2;
            }
            const std::string value = argv[++i];
            if (value.empty() || value.front() == '-') {
                std::cerr << argument << " requires a nonnegative integer.\n";
                return 2;
            }
            try {
                std::size_t parsed = 0;
                component_filter_options.minimum_supporting_poles =
                    static_cast<std::size_t>(std::stoull(value, &parsed));
                if (parsed != value.size()) {
                    throw std::invalid_argument("trailing characters");
                }
            } catch (const std::exception&) {
                std::cerr << argument << " requires a nonnegative integer.\n";
                return 2;
            }
        } else if (argument == "--max-relative-radius-jump") {
            if (i + 1 >= argc) {
                print_usage(argv[0]);
                return 2;
            }
            const std::string value = argv[++i];
            try {
                std::size_t parsed = 0;
                radius_filter_options.maximum_relative_radius_jump =
                    std::stod(value, &parsed);
                if (parsed != value.size()) {
                    throw std::invalid_argument("trailing characters");
                }
            } catch (const std::exception&) {
                std::cerr << argument << " requires a number from 0 to 1.\n";
                return 2;
            }
            if (!std::isfinite(
                    radius_filter_options.maximum_relative_radius_jump) ||
                radius_filter_options.maximum_relative_radius_jump < 0.0 ||
                radius_filter_options.maximum_relative_radius_jump > 1.0) {
                std::cerr << argument << " must be between 0 and 1.\n";
                return 2;
            }
        } else {
            positional_arguments.push_back(argument);
        }
    }

    if (positional_arguments.empty() || positional_arguments.size() > 2) {
        print_usage(argv[0]);
        return 2;
    }

    const std::filesystem::path input_path = positional_arguments[0];
    if (input_path.extension() != ".node") {
        std::cerr << "Input must be a TetGen-style .node file.\n";
        return 2;
    }

    TetGenNodeData node_data;
    std::string error;
    if (!read_tetgen_node(input_path, node_data, error)) {
        std::cerr << error << '\n';
        return 1;
    }

    Mesh explicit_surface;
    const Mesh* surface_mesh = nullptr;
    std::filesystem::path face_path = input_path;
    face_path.replace_extension(".face");
    if (std::filesystem::exists(face_path)) {
        TetGenFaceData face_data;
        if (!read_tetgen_face(face_path, face_data, error) ||
            !build_surface_mesh(node_data, face_data, explicit_surface, error)) {
            std::cerr << error << '\n';
            return 1;
        }
        surface_mesh = &explicit_surface;
        std::cout << "Loaded explicit closed surface from "
                  << face_path.string() << ".\n";
    } else {
        std::cout << "No companion " << face_path.string()
                  << " found; using the convex Delaunay boundary.\n";
    }

    std::vector<SurfaceSample> surface_samples;
    std::vector<Vec3> delaunay_points = node_data.points;
    std::vector<Vec3> sample_normals;
    std::vector<int> output_node_ids = node_data.ids;
    std::vector<int> output_boundary_markers = node_data.boundary_markers;
    int output_index_base = node_data.index_base;
    bool used_resampling = false;

    if (target_sample_count > 0) {
        if (surface_mesh == nullptr) {
            std::cerr << "--samples requires a valid companion .face surface.\n";
            return 2;
        }

        SurfaceSamplingOptions sampling_options;
        sampling_options.target_sample_count = target_sample_count;
        surface_samples = sample_surface(*surface_mesh, sampling_options);
        if (surface_samples.size() < 4) {
            std::cerr << "Surface resampling did not produce enough points.\n";
            return 1;
        }

        delaunay_points.clear();
        sample_normals.clear();
        output_boundary_markers.clear();
        output_node_ids.clear();
        delaunay_points.reserve(surface_samples.size());
        sample_normals.reserve(surface_samples.size());
        output_boundary_markers.reserve(surface_samples.size());
        output_node_ids.reserve(surface_samples.size());
        output_index_base = 1;
        for (std::size_t i = 0; i < surface_samples.size(); ++i) {
            delaunay_points.push_back(surface_samples[i].position);
            sample_normals.push_back(surface_samples[i].normal);
            output_boundary_markers.push_back(
                surface_samples[i].boundary_marker
            );
            output_node_ids.push_back(static_cast<int>(i + 1));
        }
        used_resampling = true;
        std::cout << "Generated " << surface_samples.size()
                  << " deterministic surface samples from "
                  << surface_mesh->faces.size() << " triangles.\n";
    } else if (surface_mesh != nullptr) {
        SurfaceSamplingOptions vertex_sampling_options;
        surface_samples =
            sample_surface(*surface_mesh, vertex_sampling_options);
        if (samples_match_points(surface_samples, node_data.points)) {
            sample_normals.reserve(surface_samples.size());
            for (const SurfaceSample& sample : surface_samples) {
                sample_normals.push_back(sample.normal);
            }
        } else {
            surface_samples.clear();
            std::cout << "Pole selection is disabled because the .node file "
                         "contains points not represented by the .face surface.\n";
        }
    }

    Delaunay3 delaunay;
    if (!delaunay.build(delaunay_points)) {
        std::cerr << "Could not tetrahedralize the input points. "
                  << "Check for duplicates, coplanarity, or severe degeneracy.\n";
        return 1;
    }

    if (positional_arguments.size() == 2) {
        const std::filesystem::path output_path = positional_arguments[1];
        if (output_path.extension() != ".ele") {
            std::cerr << "Output must use the .ele suffix.\n";
            return 2;
        }
        if (used_resampling) {
            std::filesystem::path sampled_node_path = output_path;
            sampled_node_path.replace_extension(".node");
            if (std::filesystem::absolute(sampled_node_path).lexically_normal() ==
                std::filesystem::absolute(input_path).lexically_normal()) {
                std::cerr << "Refusing to overwrite the input .node file with "
                             "resampled points; choose a different output base name.\n";
                return 2;
            }
            if (!write_tetgen_node(
                    sampled_node_path,
                    delaunay.points(),
                    output_boundary_markers,
                    output_index_base,
                    error)) {
                std::cerr << error << '\n';
                return 1;
            }
            std::cout << "Wrote resampled nodes to "
                      << sampled_node_path.string() << ".\n";
        }
        if (!write_tetgen_ele(
                output_path,
                delaunay,
                output_node_ids,
                output_index_base,
                error)) {
            std::cerr << error << '\n';
            return 1;
        }
        std::cout << "Wrote " << delaunay.tetrahedron_count()
                  << " tetrahedra to " << output_path.string() << ".\n";
    }

    const BoundarySurface reconstructed_boundary =
        extract_boundary_surface(delaunay);
    const auto object_faces = surface_mesh != nullptr
        ? to_surface_faces(*surface_mesh)
        : reconstructed_boundary.faces;
    const MedialAxisApproximation medial =
        build_medial_axis_approximation(delaunay, surface_mesh);
    const MedialSheetApproximation sheets =
        build_medial_sheet_approximation(delaunay, surface_mesh);
    PoleSelectionResult pole_selection;
    if (surface_mesh != nullptr &&
        samples_match_points(surface_samples, delaunay.points())) {
        pole_selection = select_inward_poles(
            delaunay,
            surface_samples,
            *surface_mesh
        );
    }
    const SurfaceFeatureField surface_features =
        estimate_surface_feature_field(
            delaunay.points(),
            sample_normals,
            pole_selection
        );
    const MedialComplex supported_complex =
        build_validated_medial_complex(
            delaunay,
            surface_mesh,
            pole_selection,
            medial_complex_options,
            &surface_features
        );
    const RadiusContinuityFilterResult radius_filter =
        filter_medial_complex_radius_continuity(
            supported_complex,
            radius_filter_options
        );
    const MedialComplexFilterResult component_filter =
        filter_medial_complex_components(
            radius_filter.retained,
            pole_selection,
            component_filter_options
        );
    const MedialComplex& component_filtered_complex =
        component_filter.retained;

    std::vector<ResolutionMedialComplex> stability_runs;
    ResolutionMedialComplex base_stability_run;
    base_stability_run.sample_count = delaunay.point_count();
    base_stability_run.complex = component_filtered_complex;
    base_stability_run.component_metrics =
        component_filter.retained_component_metrics;
    stability_runs.push_back(std::move(base_stability_run));

    if (cross_resolution_enabled &&
        surface_mesh != nullptr &&
        target_sample_count > 0) {
        const std::array<std::size_t, 3> stability_resolutions{
            40, 80, 160
        };
        for (std::size_t resolution : stability_resolutions) {
            if (resolution == delaunay.point_count()) {
                continue;
            }
            ResolutionMedialComplex run;
            if (!build_resolution_stability_run(
                    *surface_mesh,
                    resolution,
                    medial_complex_options,
                    radius_filter_options,
                    component_filter_options,
                    run,
                    error)) {
                std::cerr << error << '\n';
                return 1;
            }
            stability_runs.push_back(std::move(run));
        }
    }
    CrossResolutionStabilityResult stability =
        analyze_cross_resolution_stability(stability_runs);
    const TriangularHoleSealResult triangular_hole_seal =
        seal_isolated_triangular_holes(
            supported_complex,
            stability.retained
        );
    if (!triangular_hole_seal.restored_triangles.empty()) {
        stability.retained = triangular_hole_seal.sealed;
        for (std::size_t reference_triangle :
             triangular_hole_seal.reference_triangle_indices) {
            const std::size_t source_component =
                reference_triangle <
                        supported_complex.triangle_components.size()
                    ? supported_complex.
                        triangle_components[reference_triangle]
                    : 0;
            stability.retained_triangle_stabilities.push_back(0.0);
            stability.retained_component_stabilities.push_back(0.0);
            stability.retained_triangle_source_components.push_back(
                source_component
            );
            stability.retained_triangle_repair_flags.push_back(1.0);
        }
    }
    const MedialComplex& validated_complex = stability.retained;
    const std::vector<MedialComponentMetrics> final_component_metrics =
        analyze_medial_components(validated_complex, pole_selection);

    if (object_faces.empty()) {
        std::cerr << "The tetrahedralization has no boundary surface to display.\n";
        return 1;
    }

    const auto object_points = to_glm_points(
        surface_mesh != nullptr ? surface_mesh->vertices : delaunay.points()
    );
    const auto sample_points = to_glm_points(delaunay.points());
    const auto medial_points = to_glm_points(medial.centers);
    const auto sheet_points = to_glm_points(sheets.vertices);
    const auto validated_sheet_points =
        to_glm_points(validated_complex.vertices);

    std::cout << "Loaded " << delaunay.point_count() << " samples, generated "
              << delaunay.tetrahedron_count() << " tetrahedra, "
              << object_faces.size() << " boundary triangles, "
              << medial.centers.size() << " interior Voronoi vertices, "
              << medial.edges.size() << " contained Voronoi edges (rejected "
              << medial.rejected_exterior_edge_count
              << " exterior-crossing edges), and "
              << sheets.polygon_count << " medial-sheet polygons ("
              << sheets.triangles.size() << " triangles).\n";
    std::cout << "Selected " << pole_selection.poles.size()
              << " unique inward poles; "
              << pole_selection.validated_count()
              << " passed medial-ball validation.\n";
    if (!surface_features.resolutions.empty()) {
        const auto resolution_range = std::minmax_element(
            surface_features.resolutions.begin(),
            surface_features.resolutions.end()
        );
        const auto lfs_range = std::minmax_element(
            surface_features.local_feature_sizes.begin(),
            surface_features.local_feature_sizes.end()
        );
        const auto density_range = std::minmax_element(
            surface_features.sampling_densities.begin(),
            surface_features.sampling_densities.end()
        );
        const std::size_t anchor_count =
            static_cast<std::size_t>(std::count(
                surface_features.pole_anchored.begin(),
                surface_features.pole_anchored.end(),
                true
            ));
        std::cout << "Estimated local surface resolution ["
                  << *resolution_range.first << ", "
                  << *resolution_range.second << "], LFS ["
                  << *lfs_range.first << ", " << *lfs_range.second
                  << "], and h/LFS density [" << *density_range.first
                  << ", " << *density_range.second << "] using "
                  << anchor_count << " pole-anchored samples.\n";
    }
    std::cout << "Built " << supported_complex.accepted_polygon_count
              << " geometrically eligible sheet polygons ("
              << supported_complex.triangles.size() << " triangles in "
              << supported_complex.component_count << " component"
              << (supported_complex.component_count == 1 ? "" : "s")
              << ") using "
              << resolved_medial_support_ring_count(
                     delaunay.point_count(),
                     medial_complex_options
                 )
              << " support ring"
              << (resolved_medial_support_ring_count(
                      delaunay.point_count(),
                      medial_complex_options
                  ) == 1
                      ? ""
                      : "s")
              << "; pole evidence is stored as a weight. Restored "
              << supported_complex.
                     topology_restored_candidate_triangles.size()
              << " triangles from "
              << supported_complex.
                     topology_restored_candidate_patch_count
              << " enclosed rejected-face patch"
              << (supported_complex.
                          topology_restored_candidate_patch_count == 1
                      ? ""
                      : "es")
              << ".\n";
    std::cout << "Radius-continuity filtering measured a maximum relative "
                 "jump of "
              << radius_filter.maximum_observed_radius_jump
              << " and down-weighted "
              << radius_filter.flagged_triangles.size()
              << " triangles without deleting them, using "
              << (radius_filter_options.use_adaptive_lfs_thresholds
                      ? "adaptive"
                      : "fixed")
              << " thresholds ["
              << radius_filter.minimum_applied_threshold << ", "
              << radius_filter.maximum_applied_threshold
              << "].\n";
    for (const MedialComponentMetrics& component :
         component_filter.source_component_metrics) {
        std::cout << "  Component " << component.component << ": "
                  << component.triangle_count << " triangles, area "
                  << component.area << " ("
                  << 100.0 * component.area_fraction
                  << "%), mean confidence "
                  << component.mean_confidence
                  << ", mean/max radius jump "
                  << component.mean_radius_jump << "/"
                  << component.maximum_radius_jump << ", "
                  << component.supporting_pole_count
                  << " supporting poles, boundary length "
                  << component.boundary_length << ", seam length "
                  << component.seam_length << ", boundary loops "
                  << component.boundary_loop_count << " ("
                  << component.artificial_boundary_loop_count
                  << " artificial, "
                  << 100.0 * component.artificial_boundary_fraction
                  << "% of boundary length) -> "
                  << (component.retained ? "retained" : "removed")
                  << ".\n";
    }
    std::cout << "Component filtering retained "
              << component_filtered_complex.component_count << " component"
              << (component_filtered_complex.component_count == 1 ? "" : "s")
              << " and "
              << component_filtered_complex.triangles.size()
              << " triangles; removed "
              << component_filter.removed_triangles.size()
              << " triangles. The retained complex has "
              << component_filtered_complex.boundary_edges.size()
              << " boundary edges, "
              << component_filtered_complex.seam_edges.size()
              << " non-manifold seam edges, and "
              << component_filtered_complex.junction_vertices.size()
              << " junction vertices, "
              << component_filtered_complex.boundary_loops.size()
              << " detected boundary loops, and "
              << component_filtered_complex.
                     artificial_boundary_edges.size()
              << " unresolved artificial boundary edges.\n";
    std::cout << "Cross-resolution stability compared sample counts";
    for (std::size_t count : stability.compared_sample_counts) {
        std::cout << ' ' << count;
    }
    std::cout << "; retained " << validated_complex.triangles.size()
              << " triangles and removed "
              << stability.removed_triangles.size()
              << " triangles belonging to whole sheet components that "
                 "lacked support in at least half of the resolutions. "
                 "Topology repair restored "
              << stability.topology_repaired_triangles.size()
              << " triangles in "
              << stability.topology_repaired_patch_count
              << " topology-damaging gap"
              << (stability.topology_repaired_patch_count == 1 ? "" : "s")
              << ".\n";
    if (!stability.retained_triangle_stabilities.empty()) {
        const auto triangle_stability_range = std::minmax_element(
            stability.retained_triangle_stabilities.begin(),
            stability.retained_triangle_stabilities.end()
        );
        const auto component_stability_range = std::minmax_element(
            stability.retained_component_stabilities.begin(),
            stability.retained_component_stabilities.end()
        );
        std::cout << "Retained triangle stability range ["
                  << *triangle_stability_range.first << ", "
                  << *triangle_stability_range.second
                  << "] and source-component stability range ["
                  << *component_stability_range.first << ", "
                  << *component_stability_range.second << "].\n";
    }
    std::cout << "Final triangular-hole invariant restored "
              << triangular_hole_seal.restored_triangles.size()
              << " isolated source triangle"
              << (triangular_hole_seal.restored_triangles.size() == 1
                      ? ""
                      : "s")
              << ".\n";
    if (sheets.triangles.empty()) {
        std::cout << "No sheet passed the current 110-degree contact-angle filter.\n";
    }
    std::cout << std::flush;

    if (no_gui) {
        return 0;
    }

    polyscope::options::programName = "Medial Axis 3D Viewer";
    polyscope::options::autocenterStructures = true;
    polyscope::options::transparencyMode = polyscope::TransparencyMode::Pretty;
    polyscope::init();

    auto* surface = polyscope::registerSurfaceMesh(
        "reconstructed object",
        object_points,
        object_faces
    );
    surface->setSurfaceColor(glm::vec3(0.62f, 0.76f, 0.93f));
    surface->setEdgeWidth(0.7);
    std::vector<double> transparency(object_points.size(), 0.32);
    auto* transparency_quantity =
        surface->addVertexScalarQuantity("surface transparency", transparency);
    transparency_quantity->setMapRange({0.0, 1.0});
    surface->setTransparencyQuantity(transparency_quantity);
    if (surface_mesh != nullptr) {
        surface->addVertexVectorQuantity(
            "surface normals",
            to_glm_points(surface_mesh->vertex_normals)
        )->setEnabled(false);
    }

    auto* samples = polyscope::registerPointCloud("surface samples", sample_points);
    samples->setPointColor(glm::vec3(0.12f, 0.28f, 0.52f));
    samples->setPointRadius(0.004, false);
    samples->setEnabled(false);
    if (sample_normals.size() == delaunay.point_count()) {
        samples->addVectorQuantity(
            "sample normals",
            to_glm_points(sample_normals)
        )->setEnabled(false);
    }
    if (surface_features.resolutions.size() == delaunay.point_count()) {
        samples->addScalarQuantity(
            "surface resolution",
            surface_features.resolutions
        );
        samples->addScalarQuantity(
            "estimated LFS",
            surface_features.local_feature_sizes
        );
        samples->addScalarQuantity(
            "sampling density h/LFS",
            surface_features.sampling_densities
        );
    }

    if (!sheet_points.empty() && !sheets.triangles.empty()) {
        auto* medial_sheets = polyscope::registerSurfaceMesh(
            "medial sheets",
            sheet_points,
            sheets.triangles
        );
        medial_sheets->setSurfaceColor(glm::vec3(0.96f, 0.20f, 0.08f));
        medial_sheets->setEdgeColor(glm::vec3(0.45f, 0.03f, 0.01f));
        medial_sheets->setEdgeWidth(0.8);
        medial_sheets->setBackFacePolicy(polyscope::BackFacePolicy::Identical);
        medial_sheets->setEnabled(validated_complex.triangles.empty());
    }

    if (!validated_sheet_points.empty() &&
        !validated_complex.triangles.empty()) {
        auto* validated_sheets = polyscope::registerSurfaceMesh(
            "pole-supported medial sheets",
            validated_sheet_points,
            validated_complex.triangles
        );
        validated_sheets->setSurfaceColor(glm::vec3(0.98f, 0.48f, 0.06f));
        validated_sheets->setEdgeColor(glm::vec3(0.48f, 0.12f, 0.01f));
        validated_sheets->setEdgeWidth(0.9);
        validated_sheets->setBackFacePolicy(
            polyscope::BackFacePolicy::Identical
        );
        validated_sheets->addFaceScalarQuantity(
            "sheet confidence",
            validated_complex.triangle_confidences
        );
        if (validated_complex.triangle_support_weights.size() ==
            validated_complex.triangles.size()) {
            validated_sheets->addFaceScalarQuantity(
                "pole support weight",
                validated_complex.triangle_support_weights
            );
        }
        validated_sheets->addFaceScalarQuantity(
            "contact angle",
            validated_complex.triangle_contact_angles
        );
        validated_sheets->addVertexScalarQuantity(
            "medial radius",
            validated_complex.vertex_radii
        );
        validated_sheets->addVertexScalarQuantity(
            "surface resolution",
            validated_complex.vertex_surface_resolutions
        );
        validated_sheets->addVertexScalarQuantity(
            "estimated LFS",
            validated_complex.vertex_local_feature_sizes
        );
        validated_sheets->addVertexScalarQuantity(
            "sampling density h/LFS",
            validated_complex.vertex_sampling_densities
        );
        validated_sheets->addFaceScalarQuantity(
            "relative radius jump",
            validated_complex.triangle_radius_jumps
        );
        validated_sheets->addFaceScalarQuantity(
            "adaptive radius-jump threshold",
            validated_complex.triangle_radius_jump_thresholds
        );
        validated_sheets->addVertexScalarQuantity(
            "cross-resolution vertex stability",
            stability.vertex_stabilities
        );
        validated_sheets->addFaceScalarQuantity(
            "cross-resolution stability",
            stability.retained_triangle_stabilities
        );
        validated_sheets->addFaceScalarQuantity(
            "source component stability",
            stability.retained_component_stabilities
        );
        validated_sheets->addFaceScalarQuantity(
            "topology-repaired gap",
            stability.retained_triangle_repair_flags
        );
        std::vector<double> sheet_components;
        sheet_components.reserve(
            validated_complex.triangle_components.size()
        );
        for (std::size_t component :
             validated_complex.triangle_components) {
            sheet_components.push_back(static_cast<double>(component));
        }
        validated_sheets->addFaceScalarQuantity(
            "sheet component",
            sheet_components
        );
        std::vector<double> component_area_fractions;
        std::vector<double> component_mean_confidences;
        std::vector<double> component_pole_counts;
        std::vector<double> component_mean_radius_jumps;
        std::vector<double> component_maximum_radius_jumps;
        component_area_fractions.reserve(
            validated_complex.triangles.size()
        );
        component_mean_confidences.reserve(
            validated_complex.triangles.size()
        );
        component_pole_counts.reserve(
            validated_complex.triangles.size()
        );
        component_mean_radius_jumps.reserve(
            validated_complex.triangles.size()
        );
        component_maximum_radius_jumps.reserve(
            validated_complex.triangles.size()
        );
        for (std::size_t component :
             validated_complex.triangle_components) {
            const MedialComponentMetrics& metrics =
                final_component_metrics[component];
            component_area_fractions.push_back(metrics.area_fraction);
            component_mean_confidences.push_back(metrics.mean_confidence);
            component_pole_counts.push_back(
                static_cast<double>(metrics.supporting_pole_count)
            );
            component_mean_radius_jumps.push_back(
                metrics.mean_radius_jump
            );
            component_maximum_radius_jumps.push_back(
                metrics.maximum_radius_jump
            );
        }
        validated_sheets->addFaceScalarQuantity(
            "component area fraction",
            component_area_fractions
        );
        validated_sheets->addFaceScalarQuantity(
            "component mean confidence",
            component_mean_confidences
        );
        validated_sheets->addFaceScalarQuantity(
            "component supporting poles",
            component_pole_counts
        );
        validated_sheets->addFaceScalarQuantity(
            "component mean radius jump",
            component_mean_radius_jumps
        );
        validated_sheets->addFaceScalarQuantity(
            "component maximum radius jump",
            component_maximum_radius_jumps
        );

        if (!validated_complex.boundary_edges.empty()) {
            auto* boundaries = polyscope::registerCurveNetwork(
                "medial sheet boundaries",
                validated_sheet_points,
                validated_complex.boundary_edges
            );
            boundaries->setColor(glm::vec3(0.08f, 0.55f, 0.72f));
            boundaries->setRadius(0.003, false);
            boundaries->setEnabled(false);
        }
        if (!validated_complex.termination_edges.empty()) {
            auto* terminations = polyscope::registerCurveNetwork(
                "validated sheet terminations",
                validated_sheet_points,
                validated_complex.termination_edges
            );
            terminations->setColor(glm::vec3(0.10f, 0.72f, 0.32f));
            terminations->setRadius(0.004, false);
            terminations->setEnabled(false);
        }
        if (!validated_complex.artificial_boundary_edges.empty()) {
            auto* artificial_boundaries =
                polyscope::registerCurveNetwork(
                    "unresolved artificial sheet boundaries",
                    validated_sheet_points,
                    validated_complex.artificial_boundary_edges
                );
            artificial_boundaries->setColor(
                glm::vec3(0.92f, 0.08f, 0.08f)
            );
            artificial_boundaries->setRadius(0.006, false);
            artificial_boundaries->setEnabled(false);
        }
        if (!validated_complex.seam_edges.empty()) {
            auto* seams = polyscope::registerCurveNetwork(
                "medial sheet seams",
                validated_sheet_points,
                validated_complex.seam_edges
            );
            seams->setColor(glm::vec3(0.82f, 0.05f, 0.62f));
            seams->setRadius(0.005, false);
        }
    }

    if (!validated_sheet_points.empty() &&
        !validated_complex.triangles.empty() &&
        stability.compared_sample_counts.size() > 1) {
        auto* stability_sheets = polyscope::registerSurfaceMesh(
            "cross-resolution stability",
            validated_sheet_points,
            validated_complex.triangles
        );
        stability_sheets->setSurfaceColor(
            glm::vec3(0.18f, 0.72f, 0.42f)
        );
        stability_sheets->setEdgeColor(glm::vec3(0.03f, 0.28f, 0.13f));
        stability_sheets->setEdgeWidth(0.8);
        stability_sheets->setBackFacePolicy(
            polyscope::BackFacePolicy::Identical
        );
        stability_sheets->addFaceScalarQuantity(
            "triangle stability",
            stability.retained_triangle_stabilities
        );
        stability_sheets->addFaceScalarQuantity(
            "component stability",
            stability.retained_component_stabilities
        );
        stability_sheets->addFaceScalarQuantity(
            "topology-repaired gap",
            stability.retained_triangle_repair_flags
        );
        stability_sheets->addVertexScalarQuantity(
            "vertex stability",
            stability.vertex_stabilities
        );
        stability_sheets->setEnabled(false);
    }

    if (!stability.topology_repaired_triangles.empty()) {
        auto* repaired_holes = polyscope::registerSurfaceMesh(
            "topology-repaired stability gaps",
            to_glm_points(component_filtered_complex.vertices),
            stability.topology_repaired_triangles
        );
        repaired_holes->setSurfaceColor(
            glm::vec3(0.10f, 0.68f, 0.82f)
        );
        repaired_holes->setEdgeColor(glm::vec3(0.02f, 0.24f, 0.32f));
        repaired_holes->setEdgeWidth(0.9);
        repaired_holes->setBackFacePolicy(
            polyscope::BackFacePolicy::Identical
        );
        repaired_holes->addFaceScalarQuantity(
            "original triangle stability",
            stability.topology_repaired_triangle_stabilities
        );
        std::vector<double> repaired_source_components;
        repaired_source_components.reserve(
            stability.
                topology_repaired_triangle_source_components.size()
        );
        for (std::size_t component :
             stability.topology_repaired_triangle_source_components) {
            repaired_source_components.push_back(
                static_cast<double>(component)
            );
        }
        repaired_holes->addFaceScalarQuantity(
            "source component",
            repaired_source_components
        );
        repaired_holes->setEnabled(false);
    }

    if (!triangular_hole_seal.restored_triangles.empty()) {
        auto* sealed_triangular_holes =
            polyscope::registerSurfaceMesh(
                "sealed triangular holes",
                to_glm_points(supported_complex.vertices),
                triangular_hole_seal.restored_triangles
            );
        sealed_triangular_holes->setSurfaceColor(
            glm::vec3(0.02f, 0.78f, 0.88f)
        );
        sealed_triangular_holes->setEdgeColor(
            glm::vec3(0.01f, 0.22f, 0.28f)
        );
        sealed_triangular_holes->setEdgeWidth(1.0);
        sealed_triangular_holes->setBackFacePolicy(
            polyscope::BackFacePolicy::Identical
        );
        sealed_triangular_holes->setEnabled(false);
    }

    if (!supported_complex.
            topology_restored_candidate_triangles.empty()) {
        auto* restored_candidates = polyscope::registerSurfaceMesh(
            "restored rejected Voronoi faces",
            to_glm_points(supported_complex.vertices),
            supported_complex.topology_restored_candidate_triangles
        );
        restored_candidates->setSurfaceColor(
            glm::vec3(0.05f, 0.78f, 0.62f)
        );
        restored_candidates->setEdgeColor(
            glm::vec3(0.01f, 0.26f, 0.20f)
        );
        restored_candidates->setEdgeWidth(1.0);
        restored_candidates->setBackFacePolicy(
            polyscope::BackFacePolicy::Identical
        );
        restored_candidates->setEnabled(false);
    }

    if (!stability.removed_triangles.empty()) {
        auto* removed_resolution_sheets =
            polyscope::registerSurfaceMesh(
                "removed resolution-unstable sheets",
                to_glm_points(component_filtered_complex.vertices),
                stability.removed_triangles
            );
        removed_resolution_sheets->setSurfaceColor(
            glm::vec3(0.52f, 0.16f, 0.72f)
        );
        removed_resolution_sheets->setEdgeColor(
            glm::vec3(0.20f, 0.03f, 0.30f)
        );
        removed_resolution_sheets->setEdgeWidth(0.9);
        removed_resolution_sheets->setBackFacePolicy(
            polyscope::BackFacePolicy::Identical
        );
        removed_resolution_sheets->addFaceScalarQuantity(
            "triangle stability",
            stability.removed_triangle_stabilities
        );
        std::vector<double> removed_source_components;
        removed_source_components.reserve(
            stability.removed_triangle_source_components.size()
        );
        for (std::size_t component :
             stability.removed_triangle_source_components) {
            removed_source_components.push_back(
                static_cast<double>(component)
            );
        }
        removed_resolution_sheets->addFaceScalarQuantity(
            "source component",
            removed_source_components
        );
        removed_resolution_sheets->setEnabled(false);
    }

    if (!radius_filter.topology_repaired_triangles.empty()) {
        auto* repaired_radius_gaps = polyscope::registerSurfaceMesh(
            "topology-repaired radius gaps",
            to_glm_points(supported_complex.vertices),
            radius_filter.topology_repaired_triangles
        );
        repaired_radius_gaps->setSurfaceColor(
            glm::vec3(0.08f, 0.72f, 0.78f)
        );
        repaired_radius_gaps->setEdgeColor(
            glm::vec3(0.02f, 0.25f, 0.30f)
        );
        repaired_radius_gaps->setEdgeWidth(0.9);
        repaired_radius_gaps->setBackFacePolicy(
            polyscope::BackFacePolicy::Identical
        );
        repaired_radius_gaps->addFaceScalarQuantity(
            "relative radius jump",
            radius_filter.topology_repaired_triangle_radius_jumps
        );
        repaired_radius_gaps->addFaceScalarQuantity(
            "adaptive radius-jump threshold",
            radius_filter.
                topology_repaired_triangle_radius_jump_thresholds
        );
        repaired_radius_gaps->setEnabled(false);
    }

    if (!radius_filter.flagged_triangles.empty()) {
        auto* flagged_radius_triangles = polyscope::registerSurfaceMesh(
            "radius-discontinuity weights",
            to_glm_points(supported_complex.vertices),
            radius_filter.flagged_triangles
        );
        flagged_radius_triangles->setSurfaceColor(
            glm::vec3(0.82f, 0.12f, 0.12f)
        );
        flagged_radius_triangles->setEdgeColor(
            glm::vec3(0.35f, 0.02f, 0.02f)
        );
        flagged_radius_triangles->setEdgeWidth(0.9);
        flagged_radius_triangles->setBackFacePolicy(
            polyscope::BackFacePolicy::Identical
        );
        flagged_radius_triangles->addFaceScalarQuantity(
            "relative radius jump",
            radius_filter.flagged_triangle_radius_jumps
        );
        flagged_radius_triangles->addFaceScalarQuantity(
            "adaptive radius-jump threshold",
            radius_filter.flagged_triangle_radius_jump_thresholds
        );
        flagged_radius_triangles->setEnabled(false);
    }

    if (!component_filter.removed_triangles.empty()) {
        auto* removed_components = polyscope::registerSurfaceMesh(
            "removed unstable sheet components",
            to_glm_points(supported_complex.vertices),
            component_filter.removed_triangles
        );
        removed_components->setSurfaceColor(
            glm::vec3(0.48f, 0.48f, 0.50f)
        );
        removed_components->setEdgeColor(glm::vec3(0.20f, 0.20f, 0.22f));
        removed_components->setEdgeWidth(0.8);
        removed_components->setBackFacePolicy(
            polyscope::BackFacePolicy::Identical
        );
        removed_components->addFaceScalarQuantity(
            "removed confidence",
            component_filter.removed_triangle_confidences
        );
        removed_components->addFaceScalarQuantity(
            "relative radius jump",
            component_filter.removed_triangle_radius_jumps
        );
        removed_components->addFaceScalarQuantity(
            "adaptive radius-jump threshold",
            component_filter.removed_triangle_radius_jump_thresholds
        );
        std::vector<double> removed_component_ids;
        removed_component_ids.reserve(
            component_filter.removed_triangle_source_components.size()
        );
        for (std::size_t component :
             component_filter.removed_triangle_source_components) {
            removed_component_ids.push_back(
                static_cast<double>(component)
            );
        }
        removed_components->addFaceScalarQuantity(
            "original component",
            removed_component_ids
        );
        removed_components->setEnabled(false);
    }

    if (!validated_complex.junction_vertices.empty()) {
        std::vector<Vec3> junction_positions;
        junction_positions.reserve(
            validated_complex.junction_vertices.size()
        );
        for (std::size_t vertex : validated_complex.junction_vertices) {
            junction_positions.push_back(validated_complex.vertices[vertex]);
        }
        auto* junctions = polyscope::registerPointCloud(
            "medial sheet junctions",
            to_glm_points(junction_positions)
        );
        junctions->setPointColor(glm::vec3(0.95f, 0.05f, 0.65f));
        junctions->setPointRadius(0.011, false);
    }

    if (!validated_complex.rejected_face_centers.empty()) {
        auto* rejected = polyscope::registerPointCloud(
            "rejected sheet candidates",
            to_glm_points(validated_complex.rejected_face_centers)
        );
        rejected->setPointColor(glm::vec3(0.48f, 0.48f, 0.48f));
        rejected->setPointRadius(0.005, false);
        rejected->addScalarQuantity(
            "candidate confidence",
            validated_complex.rejected_face_confidences
        );
        rejected->setEnabled(false);
    }

    if (!medial_points.empty()) {
        auto* candidates =
            polyscope::registerPointCloud("medial candidates", medial_points);
        candidates->setPointColor(glm::vec3(0.95f, 0.22f, 0.12f));
        candidates->setPointRadius(0.008, false);
        candidates->setEnabled(false);

        if (!medial.edges.empty()) {
            auto* axis = polyscope::registerCurveNetwork(
                "medial axis approximation",
                medial_points,
                medial.edges
            );
            axis->setColor(glm::vec3(0.90f, 0.08f, 0.08f));
            axis->setRadius(0.004, false);
            axis->setEnabled(false);
        }
    }

    std::vector<Vec3> all_pole_positions;
    std::vector<Vec3> validated_pole_positions;
    std::vector<double> validated_radii;
    std::vector<double> validated_confidences;
    std::vector<Vec3> contact_positions;
    std::vector<Vec3> spoke_nodes;
    std::vector<std::array<std::size_t, 2>> spoke_edges;
    for (const MedialCandidate& pole : pole_selection.poles) {
        all_pole_positions.push_back(pole.center);
        if (!pole.validated) {
            continue;
        }

        validated_pole_positions.push_back(pole.center);
        validated_radii.push_back(pole.radius);
        validated_confidences.push_back(pole.confidence);
        const std::size_t center_index = spoke_nodes.size();
        spoke_nodes.push_back(pole.center);
        for (const SurfaceContact& contact : pole.contacts) {
            const std::size_t contact_index = spoke_nodes.size();
            spoke_nodes.push_back(contact.position);
            spoke_edges.push_back({center_index, contact_index});
            contact_positions.push_back(contact.position);
        }
    }

    if (!all_pole_positions.empty()) {
        auto* all_poles = polyscope::registerPointCloud(
            "all inward poles",
            to_glm_points(all_pole_positions)
        );
        all_poles->setPointColor(glm::vec3(0.65f, 0.25f, 0.80f));
        all_poles->setPointRadius(0.008, false);
        all_poles->setEnabled(false);
    }

    if (!validated_pole_positions.empty()) {
        const auto validated_glm = to_glm_points(validated_pole_positions);
        auto* validated_poles = polyscope::registerPointCloud(
            "validated poles",
            validated_glm
        );
        validated_poles->setPointColor(glm::vec3(0.95f, 0.72f, 0.08f));
        validated_poles->setPointRadius(0.012, false);
        validated_poles->addScalarQuantity(
            "pole radius",
            validated_radii
        );
        validated_poles->addScalarQuantity(
            "pole confidence",
            validated_confidences
        );

        auto* medial_balls = polyscope::registerPointCloud(
            "validated medial balls",
            validated_glm
        );
        medial_balls->setPointRenderMode(polyscope::PointRenderMode::Sphere);
        medial_balls->setPointColor(glm::vec3(0.95f, 0.55f, 0.10f));
        auto* radius_quantity = medial_balls->addScalarQuantity(
            "ball radius",
            validated_radii
        );
        medial_balls->setPointRadiusQuantity(radius_quantity, false);
        std::vector<double> ball_transparency(
            validated_pole_positions.size(),
            0.18
        );
        auto* ball_transparency_quantity =
            medial_balls->addScalarQuantity(
                "ball transparency",
                ball_transparency
            );
        ball_transparency_quantity->setMapRange({0.0, 1.0});
        medial_balls->setTransparencyQuantity(ball_transparency_quantity);
        medial_balls->setEnabled(
            validated_pole_positions.size() <= 8 &&
            validated_complex.triangles.empty()
        );
    }

    if (!contact_positions.empty()) {
        auto* contacts = polyscope::registerPointCloud(
            "medial ball contacts",
            to_glm_points(contact_positions)
        );
        contacts->setPointColor(glm::vec3(0.10f, 0.75f, 0.25f));
        contacts->setPointRadius(0.008, false);
        contacts->setEnabled(validated_pole_positions.size() <= 8);
    }
    if (!spoke_edges.empty()) {
        auto* spokes = polyscope::registerCurveNetwork(
            "contact spokes",
            to_glm_points(spoke_nodes),
            spoke_edges
        );
        spokes->setColor(glm::vec3(0.95f, 0.78f, 0.12f));
        spokes->setRadius(0.0025, false);
        spokes->setEnabled(validated_pole_positions.size() <= 8);
    }

    if (!screenshot_path.empty()) {
        polyscope::ScreenshotOptions screenshot_options;
        screenshot_options.transparentBackground = false;
        screenshot_options.includeUI = false;
        polyscope::screenshot(screenshot_path.string(), screenshot_options);
        std::cout << "Saved screenshot to " << screenshot_path.string() << ".\n";
        return 0;
    }

    polyscope::show();
    return 0;
}
