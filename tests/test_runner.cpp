#include "../src/core/delaunay3.h"
#include "../src/core/cross_resolution_stability.h"
#include "../src/core/filters.h"
#include "../src/core/medial_axis_approx.h"
#include "../src/core/medial_candidate.h"
#include "../src/core/medial_complex.h"
#include "../src/core/surface_sample.h"
#include "../src/core/surface_query_cpu.h"
#include "../src/core/voronoi_dual.h"
#include "../src/io/tetgen_io.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <sstream>
#include <vector>

#undef assert
#define assert(condition)                                                        \
    do {                                                                         \
        if (!(condition)) {                                                      \
            std::cerr << "Test assertion failed at " << __FILE__ << ':'          \
                      << __LINE__ << ": " #condition << '\n';                    \
            return 1;                                                            \
        }                                                                        \
    } while (false)

bool is_locally_delaunay(const medial_axis_3d::Delaunay3& delaunay) {
    const auto& points = delaunay.points();
    const auto& tetrahedra = delaunay.tetrahedra();

    auto opposite_vertex = [&](const medial_axis_3d::Tetrahedron& tet, const std::array<int, 3>& face) {
        for (int vertex : tet.vertices) {
            if (std::find(face.begin(), face.end(), vertex) == face.end()) {
                return vertex;
            }
        }
        return -1;
    };

    for (std::size_t a = 0; a < tetrahedra.size(); ++a) {
        for (std::size_t b = a + 1; b < tetrahedra.size(); ++b) {
            std::array<int, 3> shared_face{};
            std::size_t shared_count = 0;
            for (int vertex : tetrahedra[a].vertices) {
                if (std::find(tetrahedra[b].vertices.begin(), tetrahedra[b].vertices.end(), vertex) != tetrahedra[b].vertices.end()) {
                    shared_face[shared_count++] = vertex;
                }
            }
            if (shared_count != 3) {
                continue;
            }

            const int opposite_a = opposite_vertex(tetrahedra[a], shared_face);
            const int opposite_b = opposite_vertex(tetrahedra[b], shared_face);
            if (opposite_a < 0 || opposite_b < 0) {
                return false;
            }

            if (delaunay.circumsphere_contains(a, points[opposite_b]) ||
                delaunay.circumsphere_contains(b, points[opposite_a])) {
                return false;
            }
        }
    }
    return true;
}

bool has_unique_vertices(const medial_axis_3d::Delaunay3& delaunay) {
    const auto& tetrahedra = delaunay.tetrahedra();
    for (const auto& tetrahedron : tetrahedra) {
        std::array<int, 4> vertices = tetrahedron.vertices;
        std::sort(vertices.begin(), vertices.end());
        if (vertices[0] == vertices[1] || vertices[1] == vertices[2] || vertices[2] == vertices[3]) {
            return false;
        }
    }
    return true;
}

int main() {
    using namespace medial_axis_3d;

    Delaunay3 delaunay;
    assert(delaunay.point_count() == 0);
    assert(delaunay.tetrahedron_count() == 0);

    const std::vector<Vec3> samples = {
        Vec3(0.0, 0.0, 0.0),
        Vec3(1.0, 0.0, 0.0),
        Vec3(0.0, 1.0, 0.0),
        Vec3(0.0, 0.0, 1.0),
        Vec3(0.2, 0.2, 0.2)
    };

    for (const auto& sample : samples) {
        assert(delaunay.insert(sample));
    }

    Delaunay3 initial_tetra;
    for (const auto& sample : std::vector<Vec3>{
        Vec3(0.0, 0.0, 0.0),
        Vec3(1.0, 0.0, 0.0),
        Vec3(0.0, 1.0, 0.0),
        Vec3(0.0, 0.0, 1.0)
    }) {
        initial_tetra.insert(sample);
    }

    const std::array<int, 4> expected_initial_vertices{0, 1, 2, 3};
    assert(initial_tetra.tetrahedra()[0].vertices == expected_initial_vertices);

    assert(delaunay.point_count() == 5);
    assert(delaunay.tetrahedron_count() == 4);
    assert(delaunay.tetrahedra()[0].vertices[0] == 0);
    assert(delaunay.tetrahedra()[3].vertices[3] == 4);

    Delaunay3 cavity_test;
    for (const auto& sample : std::vector<Vec3>{
        Vec3(0.0, 0.0, 0.0),
        Vec3(1.0, 0.0, 0.0),
        Vec3(0.0, 1.0, 0.0),
        Vec3(0.0, 0.0, 1.0)
    }) {
        cavity_test.insert(sample);
    }

    const auto cavity = cavity_test.cavity_indices(Vec3(0.2, 0.2, 0.2));
    assert(!cavity.empty());

    // Check adjacency: each tetrahedron should have neighbors recorded (or -1 for boundary)
    bool any_neighbor_present = false;
    for (std::size_t t = 0; t < delaunay.tetrahedron_count(); ++t) {
        const auto& tet = delaunay.tetrahedra()[t];
        for (int n = 0; n < 4; ++n) {
            const int neigh = tet.neighbors[n];
            if (neigh != -1) {
                any_neighbor_present = true;
                assert(neigh >= 0 && neigh < static_cast<int>(delaunay.tetrahedron_count()));
            }
        }
    }
    assert(any_neighbor_present);
    assert(is_locally_delaunay(delaunay));

    Delaunay3 repair_test;
    for (const auto& sample : std::vector<Vec3>{
        Vec3(0.0, 0.0, 0.0),
        Vec3(1.0, 0.0, 0.0),
        Vec3(0.0, 1.0, 0.0),
        Vec3(0.0, 0.0, 1.0),
        Vec3(0.2, 0.2, 0.2),
        Vec3(0.3, 0.2, 0.05)
    }) {
        assert(repair_test.insert(sample));
    }
    assert(repair_test.tetrahedron_count() > 0);
    assert(repair_test.is_well_formed());
    assert(has_unique_vertices(repair_test));
    assert(is_locally_delaunay(repair_test));

    Delaunay3 geometry_test;
    for (const auto& sample : std::vector<Vec3>{
        Vec3(0.0, 0.0, 0.0),
        Vec3(1.0, 0.0, 0.0),
        Vec3(0.0, 1.0, 0.0),
        Vec3(0.0, 0.0, 1.0),
        Vec3(0.15, 0.15, 0.15),
        Vec3(0.8, 0.2, 0.2),
        Vec3(0.2, 0.8, 0.2),
        Vec3(0.2, 0.2, 0.8)
    }) {
        assert(geometry_test.insert(sample));
    }
    assert(geometry_test.tetrahedron_count() > 0);
    assert(geometry_test.is_well_formed());
    assert(has_unique_vertices(geometry_test));
    assert(is_locally_delaunay(geometry_test));

    const VoronoiDual dual = build_voronoi_dual(delaunay);
    assert(dual.cells.size() == delaunay.point_count());
    assert(!dual.cells.empty());
    assert(dual.cells[0].owner_point == 0);
    assert(!dual.cells[0].vertices.empty());

    std::istringstream node_input(
        "# five 3D points\n"
        "5 3 1 1\n"
        "1 0 0 0 10.0 1\n"
        "2 1 0 0 20.0 1\n"
        "3 0 1 0 30.0 1\n"
        "4 0 0 1 40.0 1\n"
        "5 0.2 0.2 0.2 50.0 0\n"
    );
    TetGenNodeData node_data;
    std::string io_error;
    assert(read_tetgen_node(node_input, node_data, io_error));
    assert(node_data.points.size() == 5);
    assert(node_data.ids.front() == 1);
    assert(node_data.attributes[4][0] == 50.0);
    assert(node_data.boundary_markers[4] == 0);

    std::istringstream face_input(
        "# tetrahedron boundary\n"
        "4 1\n"
        "1 1 3 2 7\n"
        "2 1 2 4 7\n"
        "3 2 3 4 7\n"
        "4 3 1 4 7\n"
    );
    TetGenFaceData face_data;
    assert(read_tetgen_face(face_input, face_data, io_error));
    assert(face_data.node_ids.size() == 4);
    assert(face_data.boundary_markers[0] == 7);

    Mesh tetrahedron_surface;
    assert(build_surface_mesh(
        node_data,
        face_data,
        tetrahedron_surface,
        io_error
    ));
    assert(tetrahedron_surface.faces.size() == 4);
    assert(tetrahedron_surface.face_normals.size() == 4);
    assert(tetrahedron_surface.vertex_normals.size() == 5);
    assert(point_inside_mesh(tetrahedron_surface, Vec3(0.1, 0.1, 0.1)));
    assert(!point_inside_mesh(tetrahedron_surface, Vec3(2.0, 2.0, 2.0)));

    CpuSurfaceQueryBackend surface_queries(tetrahedron_surface);
    assert(std::string(surface_queries.name()) == "cpu-reference");

    const std::vector<Vec3> containment_queries = {
        Vec3(0.1, 0.1, 0.1),
        Vec3(2.0, 2.0, 2.0)
    };
    std::vector<SurfaceQueryClassification> containment_results;
    surface_queries.points_inside(
        containment_queries,
        containment_results
    );
    assert(containment_results.size() == containment_queries.size());
    assert(containment_results[0] == surface_query_true);
    assert(containment_results[1] == surface_query_false);

    const std::vector<SurfaceSegmentQuery> intersection_queries = {
        {Vec3(0.1, 0.1, 0.1), Vec3(2.0, 2.0, 2.0)},
        {Vec3(0.1, 0.1, 0.1), Vec3(0.2, 0.1, 0.1)}
    };
    std::vector<SurfaceQueryClassification> intersection_results;
    surface_queries.segments_intersect(
        intersection_queries,
        intersection_results
    );
    assert(intersection_results.size() == intersection_queries.size());
    assert(intersection_results[0] == surface_query_true);
    assert(intersection_results[1] == surface_query_false);

    const std::vector<NearestSurfaceQuery> contact_queries = {
        {Vec3(0.1, 0.1, 0.1), 0.0},
        {Vec3(2.0, 2.0, 2.0), 0.05}
    };
    std::vector<std::vector<SurfaceContact>> contact_results;
    surface_queries.nearest_surface_contacts_batch(
        contact_queries,
        contact_results
    );
    assert(contact_results.size() == contact_queries.size());
    for (std::size_t query = 0; query < contact_queries.size(); ++query) {
        const auto expected_contacts = nearest_surface_contacts(
            tetrahedron_surface,
            contact_queries[query].point,
            contact_queries[query].relative_tolerance
        );
        assert(contact_results[query].size() == expected_contacts.size());
        for (std::size_t contact = 0;
             contact < expected_contacts.size();
             ++contact) {
            assert(contact_results[query][contact].position ==
                   expected_contacts[contact].position);
            assert(contact_results[query][contact].normal ==
                   expected_contacts[contact].normal);
            assert(contact_results[query][contact].triangle ==
                   expected_contacts[contact].triangle);
            assert(contact_results[query][contact].distance ==
                   expected_contacts[contact].distance);
        }
    }

    std::vector<SurfaceQueryClassification> empty_classifications = {
        surface_query_true
    };
    surface_queries.points_inside({}, empty_classifications);
    assert(empty_classifications.empty());

    SurfaceSamplingOptions sampling_options;
    sampling_options.target_sample_count = 20;
    const auto sampled_surface =
        sample_surface(tetrahedron_surface, sampling_options);
    const auto sampled_surface_repeat =
        sample_surface(tetrahedron_surface, sampling_options);
    assert(sampled_surface.size() == 20);
    assert(sampled_surface_repeat.size() == sampled_surface.size());
    for (std::size_t i = 0; i < sampled_surface.size(); ++i) {
        const SurfaceSample& sample = sampled_surface[i];
        assert(sample.position == sampled_surface_repeat[i].position);
        assert(sample.source_triangle < tetrahedron_surface.faces.size());
        assert(std::fabs(
            sample.barycentric[0] +
            sample.barycentric[1] +
            sample.barycentric[2] - 1.0
        ) < 1e-12);
        assert(sample.barycentric[0] >= 0.0);
        assert(sample.barycentric[1] >= 0.0);
        assert(sample.barycentric[2] >= 0.0);
        assert(norm(sample.normal) > 0.99);
        assert(sample.boundary_marker == 7);
    }
    SurfaceSamplingOptions reduced_sampling_options;
    reduced_sampling_options.target_sample_count = 4;
    reduced_sampling_options.include_mesh_vertices = false;
    const auto reduced_samples =
        sample_surface(tetrahedron_surface, reduced_sampling_options);
    assert(reduced_samples.size() == 4);
    assert(std::none_of(
        reduced_samples.begin(),
        reduced_samples.end(),
        [](const SurfaceSample& sample) {
            return sample.source_vertex >= 0;
        }
    ));

    SurfaceSamplingOptions weighted_sampling_options;
    weighted_sampling_options.target_sample_count = 20;
    weighted_sampling_options.triangle_importance_weights.assign(
        tetrahedron_surface.faces.size(),
        0.0
    );
    weighted_sampling_options.triangle_importance_weights[0] = 1.0;
    const auto weighted_samples =
        sample_surface(tetrahedron_surface, weighted_sampling_options);
    assert(weighted_samples.size() == 20);
    assert(std::all_of(
        weighted_samples.begin(),
        weighted_samples.end(),
        [](const SurfaceSample& sample) {
            return sample.source_vertex >= 0 ||
                   sample.source_triangle == 0;
        }
    ));

    SurfaceFeatureField adaptive_test_features;
    adaptive_test_features.sampling_densities.assign(
        tetrahedron_surface.vertices.size(),
        0.5
    );
    const int adaptive_focus_vertex =
        tetrahedron_surface.faces[0].vertices[0];
    adaptive_test_features.sampling_densities[
        static_cast<std::size_t>(adaptive_focus_vertex)
    ] = 2.0;
    const std::vector<double> adaptive_importance =
        lfs_adaptive_triangle_importance(
            tetrahedron_surface,
            adaptive_test_features
        );
    assert(adaptive_importance.size() ==
           tetrahedron_surface.faces.size());
    assert(adaptive_importance[0] == 4.0);
    assert(std::count(
        adaptive_importance.begin(),
        adaptive_importance.end(),
        0.25
    ) == 1);

    std::vector<Vec3> sampled_positions;
    std::vector<int> sampled_markers;
    for (const SurfaceSample& sample : sampled_surface) {
        sampled_positions.push_back(sample.position);
        sampled_markers.push_back(sample.boundary_marker);
    }
    std::ostringstream sampled_node_output;
    assert(write_tetgen_node(
        sampled_node_output,
        sampled_positions,
        sampled_markers,
        1,
        io_error
    ));
    std::istringstream sampled_node_input(sampled_node_output.str());
    TetGenNodeData sampled_node_data;
    assert(read_tetgen_node(sampled_node_input, sampled_node_data, io_error));
    assert(sampled_node_data.points.size() == sampled_surface.size());
    assert(sampled_node_data.boundary_markers == sampled_markers);

    TetGenFaceData open_face_data = face_data;
    open_face_data.node_ids.pop_back();
    open_face_data.ids.pop_back();
    open_face_data.boundary_markers.pop_back();
    Mesh invalid_surface;
    assert(!build_surface_mesh(
        node_data,
        open_face_data,
        invalid_surface,
        io_error
    ));

    Mesh dented_cube_surface;
    dented_cube_surface.vertices = {
        Vec3(-1.0, -1.0, -1.0),
        Vec3( 1.0, -1.0, -1.0),
        Vec3( 1.0,  1.0, -1.0),
        Vec3(-1.0,  1.0, -1.0),
        Vec3(-1.0, -1.0,  1.0),
        Vec3( 1.0, -1.0,  1.0),
        Vec3( 1.0,  1.0,  1.0),
        Vec3(-1.0,  1.0,  1.0),
        Vec3( 0.0,  0.0,  0.2)
    };
    dented_cube_surface.faces = {
        Triangle{{0, 2, 1}},
        Triangle{{0, 3, 2}},
        Triangle{{0, 1, 5}},
        Triangle{{0, 5, 4}},
        Triangle{{1, 2, 6}},
        Triangle{{1, 6, 5}},
        Triangle{{2, 3, 7}},
        Triangle{{2, 7, 6}},
        Triangle{{3, 0, 4}},
        Triangle{{3, 4, 7}},
        Triangle{{4, 5, 8}},
        Triangle{{5, 6, 8}},
        Triangle{{6, 7, 8}},
        Triangle{{7, 4, 8}}
    };
    assert(orient_and_analyze_closed_mesh(dented_cube_surface, io_error));
    assert(point_inside_mesh(dented_cube_surface, Vec3(0.0, 0.0, 0.0)));
    assert(!point_inside_mesh(dented_cube_surface, Vec3(0.0, 0.0, 0.8)));
    const Vec3 left_of_dent(-0.8, 0.0, 0.7);
    const Vec3 right_of_dent(0.8, 0.0, 0.7);
    assert(point_inside_mesh(dented_cube_surface, left_of_dent));
    assert(point_inside_mesh(dented_cube_surface, right_of_dent));
    assert(segment_intersects_mesh_surface(
        dented_cube_surface,
        left_of_dent,
        right_of_dent
    ));
    assert(!segment_inside_mesh(
        dented_cube_surface,
        left_of_dent,
        right_of_dent
    ));
    assert(segment_inside_mesh(
        dented_cube_surface,
        Vec3(-0.8, 0.0, 0.0),
        Vec3(0.8, 0.0, 0.0)
    ));

    Delaunay3 batch_delaunay;
    assert(batch_delaunay.build(node_data.points));
    assert(batch_delaunay.point_count() == node_data.points.size());
    assert(batch_delaunay.tetrahedron_count() == 4);
    assert(batch_delaunay.is_well_formed());

    const BoundarySurface boundary = extract_boundary_surface(batch_delaunay);
    assert(boundary.faces.size() == 4);
    const MedialAxisApproximation medial =
        build_medial_axis_approximation(batch_delaunay);
    assert(medial.centers.size() == medial.radii.size());
    assert(medial.centers.size() == medial.source_tetrahedra.size());

    std::ostringstream ele_output;
    assert(write_tetgen_ele(
        ele_output,
        batch_delaunay,
        node_data.ids,
        node_data.index_base,
        io_error
    ));
    assert(ele_output.str().find("4 4 0\n") == 0);

    const std::vector<Vec3> ellipsoid_octahedron_points{
        Vec3(1.5, 0.0, 0.0),
        Vec3(-1.5, 0.0, 0.0),
        Vec3(0.0, 1.0, 0.0),
        Vec3(0.0, -1.0, 0.0),
        Vec3(0.0, 0.0, 0.75),
        Vec3(0.0, 0.0, -0.75)
    };
    Delaunay3 ellipsoid_octahedron;
    assert(ellipsoid_octahedron.build(ellipsoid_octahedron_points));
    Mesh ellipsoid_octahedron_surface;
    ellipsoid_octahedron_surface.vertices = ellipsoid_octahedron_points;
    ellipsoid_octahedron_surface.faces = {
        Triangle{{4, 0, 2}},
        Triangle{{4, 2, 1}},
        Triangle{{4, 1, 3}},
        Triangle{{4, 3, 0}},
        Triangle{{5, 2, 0}},
        Triangle{{5, 1, 2}},
        Triangle{{5, 3, 1}},
        Triangle{{5, 0, 3}}
    };
    assert(orient_and_analyze_closed_mesh(
        ellipsoid_octahedron_surface,
        io_error
    ));
    const auto voronoi_faces =
        build_interior_voronoi_faces(
            ellipsoid_octahedron,
            &ellipsoid_octahedron_surface
        );
    assert(!voronoi_faces.empty());
    for (const VoronoiFaceCandidate& face : voronoi_faces) {
        assert(face.vertices.size() == face.source_tetrahedra.size());
        for (std::size_t source : face.source_tetrahedra) {
            assert(source < ellipsoid_octahedron.tetrahedron_count());
        }
    }
    const auto medial_sheets =
        build_medial_sheet_approximation(
            ellipsoid_octahedron,
            &ellipsoid_octahedron_surface
        );
    assert(medial_sheets.polygon_count > 0);
    assert(!medial_sheets.triangles.empty());
    assert(medial_sheets.triangle_contact_angles.size() ==
           medial_sheets.triangles.size());

    const Vec3 closest = closest_point_on_triangle(
        Vec3(0.25, 0.25, 1.0),
        Vec3(0.0, 0.0, 0.0),
        Vec3(1.0, 0.0, 0.0),
        Vec3(0.0, 1.0, 0.0)
    );
    assert(norm(closest - Vec3(0.25, 0.25, 0.0)) < 1e-12);

    SurfaceSamplingOptions pole_sample_options;
    const auto octahedron_samples = sample_surface(
        ellipsoid_octahedron_surface,
        pole_sample_options
    );
    assert(octahedron_samples.size() == ellipsoid_octahedron.point_count());
    const PoleSelectionResult octahedron_poles = select_inward_poles(
        ellipsoid_octahedron,
        octahedron_samples,
        ellipsoid_octahedron_surface
    );
    assert(!octahedron_poles.poles.empty());
    assert(octahedron_poles.validated_count() == 0);
    for (const MedialCandidate& pole : octahedron_poles.poles) {
        assert(pole.radius >= 0.0);
        assert(pole.confidence >= 0.0 && pole.confidence <= 1.0);
        if (pole.validated) {
            assert(pole.contacts.size() >= 2);
            assert(pole.maximum_contact_angle_degrees >= 90.0);
        }
    }

    SurfaceSamplingOptions dense_pole_sample_options;
    dense_pole_sample_options.target_sample_count = 40;
    const auto dense_octahedron_samples = sample_surface(
        ellipsoid_octahedron_surface,
        dense_pole_sample_options
    );
    std::vector<Vec3> dense_octahedron_points;
    std::vector<Vec3> dense_octahedron_normals;
    for (const SurfaceSample& sample : dense_octahedron_samples) {
        dense_octahedron_points.push_back(sample.position);
        dense_octahedron_normals.push_back(sample.normal);
    }
    Delaunay3 dense_octahedron;
    assert(dense_octahedron.build(dense_octahedron_points));
    const PoleSelectionResult dense_octahedron_poles = select_inward_poles(
        dense_octahedron,
        dense_octahedron_samples,
        ellipsoid_octahedron_surface
    );
    assert(dense_octahedron_poles.validated_count() > 0);
    const SurfaceFeatureField dense_surface_features =
        estimate_surface_feature_field(
            dense_octahedron_points,
            dense_octahedron_normals,
            dense_octahedron_poles
        );
    assert(dense_surface_features.resolutions.size() ==
           dense_octahedron_points.size());
    assert(dense_surface_features.local_feature_sizes.size() ==
           dense_octahedron_points.size());
    assert(dense_surface_features.sampling_densities.size() ==
           dense_octahedron_points.size());
    assert(std::count(
        dense_surface_features.pole_anchored.begin(),
        dense_surface_features.pole_anchored.end(),
        true
    ) > 0);
    for (std::size_t sample = 0;
         sample < dense_surface_features.resolutions.size();
         ++sample) {
        assert(dense_surface_features.resolutions[sample] > 0.0);
        assert(dense_surface_features.local_feature_sizes[sample] >=
               0.25 * dense_surface_features.resolutions[sample] - 1e-12);
        assert(dense_surface_features.sampling_densities[sample] > 0.0);
    }
    const MedialComplex dense_medial_complex =
        build_validated_medial_complex(
            dense_octahedron,
            &ellipsoid_octahedron_surface,
            dense_octahedron_poles,
            {},
            &dense_surface_features
        );
    assert(dense_medial_complex.accepted_polygon_count > 0);
    assert(!dense_medial_complex.triangles.empty());
    assert(dense_medial_complex.vertex_radii.size() ==
           dense_medial_complex.vertices.size());
    assert(dense_medial_complex.vertex_surface_resolutions.size() ==
           dense_medial_complex.vertices.size());
    assert(dense_medial_complex.vertex_local_feature_sizes.size() ==
           dense_medial_complex.vertices.size());
    assert(dense_medial_complex.vertex_sampling_densities.size() ==
           dense_medial_complex.vertices.size());
    assert(dense_medial_complex.vertex_source_tetrahedra.size() ==
           dense_medial_complex.vertices.size());
    assert(dense_medial_complex.component_count > 0);
    assert(dense_medial_complex.triangle_confidences.size() ==
           dense_medial_complex.triangles.size());
    assert(dense_medial_complex.triangle_support_weights.size() ==
           dense_medial_complex.triangles.size());
    assert(dense_medial_complex.triangle_contact_angles.size() ==
           dense_medial_complex.triangles.size());
    assert(dense_medial_complex.triangle_radius_jumps.size() ==
           dense_medial_complex.triangles.size());
    assert(dense_medial_complex.triangle_components.size() ==
           dense_medial_complex.triangles.size());
    for (std::size_t triangle = 0;
         triangle < dense_medial_complex.triangles.size();
         ++triangle) {
        for (std::size_t vertex :
             dense_medial_complex.triangles[triangle]) {
            assert(vertex < dense_medial_complex.vertices.size());
        }
        assert(dense_medial_complex.triangle_confidences[triangle] >= 0.0);
        assert(dense_medial_complex.triangle_confidences[triangle] <= 1.0);
        assert(dense_medial_complex.triangle_radius_jumps[triangle] >= 0.0);
        assert(dense_medial_complex.triangle_radius_jumps[triangle] <= 1.0);
        assert(dense_medial_complex.triangle_components[triangle] <
               dense_medial_complex.component_count);
    }
    for (const auto& edge : dense_medial_complex.boundary_edges) {
        assert(edge[0] < edge[1]);
        assert(edge[1] < dense_medial_complex.vertices.size());
    }
    for (const auto& edge : dense_medial_complex.seam_edges) {
        assert(edge[0] < edge[1]);
        assert(edge[1] < dense_medial_complex.vertices.size());
    }
    for (std::size_t junction :
         dense_medial_complex.junction_vertices) {
        assert(junction < dense_medial_complex.vertices.size());
    }
    const RadiusContinuityFilterResult dense_radius_filtered =
        filter_medial_complex_radius_continuity(dense_medial_complex);
    assert(dense_radius_filtered.retained.triangles.size() +
               dense_radius_filtered.removed_triangles.size() ==
           dense_medial_complex.triangles.size());
    assert(dense_radius_filtered.removed_triangle_radius_jumps.size() ==
           dense_radius_filtered.removed_triangles.size());
    assert(dense_radius_filtered.
               removed_triangle_radius_jump_thresholds.size() ==
           dense_radius_filtered.removed_triangles.size());
    assert(dense_radius_filtered.retained.
               triangle_radius_jump_thresholds.size() ==
           dense_radius_filtered.retained.triangles.size());
    assert(dense_radius_filtered.minimum_applied_threshold >= 0.35);
    assert(dense_radius_filtered.maximum_applied_threshold <= 0.85);
    assert(!dense_radius_filtered.retained.triangles.empty());
    const MedialComplexFilterResult dense_filtered_complex =
        filter_medial_complex_components(
            dense_radius_filtered.retained,
            dense_octahedron_poles
        );
    assert(dense_filtered_complex.source_component_metrics.size() ==
           dense_radius_filtered.retained.component_count);
    assert(dense_filtered_complex.retained.triangles.size() +
               dense_filtered_complex.removed_triangles.size() ==
           dense_radius_filtered.retained.triangles.size());
    assert(!dense_filtered_complex.retained.triangles.empty());
    assert(dense_filtered_complex.retained_component_metrics.size() ==
           dense_filtered_complex.retained.component_count);

    MedialComplex synthetic_complex;
    synthetic_complex.vertices = {
        Vec3(0.0, 0.0, 0.0),
        Vec3(1.0, 0.0, 0.0),
        Vec3(1.0, 1.0, 0.0),
        Vec3(0.0, 1.0, 0.0),
        Vec3(2.0, 0.0, 0.0),
        Vec3(2.01, 0.0, 0.0),
        Vec3(2.0, 0.01, 0.0)
    };
    synthetic_complex.vertex_source_tetrahedra = {
        {0}, {0}, {0}, {0}, {1}, {1}, {1}
    };
    synthetic_complex.vertex_radii = {
        1.0, 1.1, 1.0, 1.05, 1.0, 1.0, 1.0
    };
    synthetic_complex.vertex_surface_resolutions.assign(7, 0.05);
    synthetic_complex.vertex_local_feature_sizes.assign(7, 1.0);
    synthetic_complex.vertex_sampling_densities.assign(7, 0.05);
    synthetic_complex.triangles = {
        {0, 1, 2},
        {0, 2, 3},
        {4, 5, 6}
    };
    synthetic_complex.triangle_confidences = {0.8, 0.8, 0.1};
    synthetic_complex.triangle_contact_angles = {140.0, 140.0, 91.0};
    detail::orient_and_label_sheet_components(synthetic_complex);
    detail::classify_complex_topology(synthetic_complex);
    for (const auto& triangle : synthetic_complex.triangles) {
        synthetic_complex.triangle_radius_jumps.push_back(
            detail::triangle_max_radius_jump(
                synthetic_complex,
                triangle
            )
        );
    }

    PoleSelectionResult synthetic_poles;
    MedialCandidate synthetic_pole;
    synthetic_pole.validated = true;
    synthetic_pole.source_tetrahedron = 0;
    synthetic_pole.confidence = 0.8;
    synthetic_poles.poles.push_back(synthetic_pole);

    MedialComplex synthetic_radius_test = synthetic_complex;
    synthetic_radius_test.vertex_radii[5] = 10.0;
    synthetic_radius_test.vertex_surface_resolutions[4] = 2.0;
    synthetic_radius_test.vertex_surface_resolutions[5] = 2.0;
    synthetic_radius_test.vertex_surface_resolutions[6] = 2.0;
    synthetic_radius_test.vertex_sampling_densities[4] = 2.0;
    synthetic_radius_test.vertex_sampling_densities[5] = 2.0;
    synthetic_radius_test.vertex_sampling_densities[6] = 2.0;
    synthetic_radius_test.triangle_radius_jumps.clear();
    for (const auto& triangle : synthetic_radius_test.triangles) {
        synthetic_radius_test.triangle_radius_jumps.push_back(
            detail::triangle_max_radius_jump(
                synthetic_radius_test,
                triangle
            )
        );
    }
    const RadiusContinuityFilterResult synthetic_radius_filtered =
        filter_medial_complex_radius_continuity(synthetic_radius_test);
    assert(synthetic_radius_filtered.retained.triangles.size() == 3);
    assert(synthetic_radius_filtered.removed_triangles.empty());
    assert(synthetic_radius_filtered.flagged_triangles.size() == 1);
    assert(synthetic_radius_filtered.flagged_triangle_radius_jumps[0] >
           0.6);
    assert(synthetic_radius_filtered.minimum_applied_threshold <
           synthetic_radius_filtered.maximum_applied_threshold);

    MedialComponentFilterOptions synthetic_filter_options;
    synthetic_filter_options.minimum_triangle_count = 2;
    synthetic_filter_options.minimum_area_fraction = 0.01;
    synthetic_filter_options.minimum_mean_confidence = 0.5;
    synthetic_filter_options.minimum_supporting_poles = 1;
    synthetic_filter_options.require_valid_sheet_boundaries = false;
    const MedialComplexFilterResult synthetic_filtered =
        filter_medial_complex_components(
            synthetic_complex,
            synthetic_poles,
            synthetic_filter_options
        );
    assert(synthetic_filtered.source_component_metrics.size() == 2);
    assert(synthetic_filtered.retained.triangles.size() == 2);
    assert(synthetic_filtered.removed_triangles.size() == 1);
    assert(synthetic_filtered.retained.component_count == 1);
    assert(synthetic_filtered.source_component_metrics[0].retained);
    assert(!synthetic_filtered.source_component_metrics[1].retained);
    assert(synthetic_filtered.source_component_metrics[0].
               supporting_pole_count == 1);
    assert(synthetic_filtered.source_component_metrics[1].
               supporting_pole_count == 0);

    // Propagated support is already represented by the component confidence.
    // A supported neighboring stratum must survive the default filter even
    // when it does not contain a pole's source tetrahedron.
    MedialComplex propagated_component = synthetic_complex;
    propagated_component.triangles = {{4, 5, 6}};
    propagated_component.triangle_confidences = {0.8};
    propagated_component.triangle_support_weights = {0.7};
    propagated_component.triangle_contact_angles = {140.0};
    propagated_component.triangle_radius_jumps = {0.02};
    propagated_component.triangle_components.clear();
    propagated_component.permitted_boundary_vertices = {4, 5, 6};
    detail::orient_and_label_sheet_components(propagated_component);
    detail::classify_complex_topology(propagated_component);
    const MedialComplexFilterResult propagated_filtered =
        filter_medial_complex_components(
            propagated_component,
            synthetic_poles
        );
    assert(propagated_filtered.retained.triangles.size() == 1);
    assert(propagated_filtered.removed_triangles.empty());

    ResolutionMedialComplex stability_base;
    stability_base.sample_count = 40;
    stability_base.complex = synthetic_complex;
    stability_base.component_metrics = analyze_medial_components(
        stability_base.complex,
        synthetic_poles
    );

    ResolutionMedialComplex stability_80;
    stability_80.sample_count = 80;
    stability_80.complex = synthetic_complex;
    stability_80.complex.triangles.resize(2);
    stability_80.complex.triangle_confidences.resize(2);
    stability_80.complex.triangle_contact_angles.resize(2);
    stability_80.complex.triangle_radius_jumps.resize(2);
    stability_80.complex.triangle_components.clear();
    stability_80.complex.boundary_edges.clear();
    stability_80.complex.seam_edges.clear();
    stability_80.complex.junction_vertices.clear();
    detail::orient_and_label_sheet_components(stability_80.complex);
    detail::classify_complex_topology(stability_80.complex);
    stability_80.component_metrics = analyze_medial_components(
        stability_80.complex,
        synthetic_poles
    );

    ResolutionMedialComplex stability_160 = stability_80;
    stability_160.sample_count = 160;
    for (std::size_t vertex = 0; vertex < 4; ++vertex) {
        stability_160.complex.vertices[vertex].z += 0.05;
    }
    stability_160.component_metrics = analyze_medial_components(
        stability_160.complex,
        synthetic_poles
    );

    CrossResolutionStabilityOptions synthetic_stability_options;
    synthetic_stability_options.remove_unstable_components = true;
    const CrossResolutionStabilityResult synthetic_stability =
        analyze_cross_resolution_stability(
            {
                stability_base,
                stability_80,
                stability_160
            },
            0,
            synthetic_stability_options
        );
    assert(synthetic_stability.compared_sample_counts.size() == 3);
    assert(synthetic_stability.retained.triangles.size() == 2);
    assert(synthetic_stability.removed_triangles.size() == 1);
    assert(synthetic_stability.retained_triangle_stabilities.size() == 2);
    assert(synthetic_stability.removed_triangle_stabilities.size() == 1);
    assert(synthetic_stability.retained_triangle_stabilities[0] == 1.0);
    assert(synthetic_stability.removed_triangle_stabilities[0] <
           0.5);
    assert(synthetic_stability.source_component_stabilities.size() == 2);
    assert(synthetic_stability.source_component_stabilities[0] == 1.0);
    assert(synthetic_stability.source_component_stabilities[1] <
           0.5);

    const CrossResolutionStabilityResult diagnostic_stability =
        analyze_cross_resolution_stability({
            stability_base,
            stability_80,
            stability_160
        });
    assert(diagnostic_stability.retained.triangles.size() ==
           stability_base.complex.triangles.size());
    assert(diagnostic_stability.removed_triangles.empty());
    assert(diagnostic_stability.retained_triangle_stabilities.size() ==
           stability_base.complex.triangles.size());

    MedialComplex enclosed_hole_complex;
    enclosed_hole_complex.vertices = {
        Vec3(0.0, 0.0, 0.0),
        Vec3(1.0, 0.0, 0.0),
        Vec3(0.0, 1.0, 0.0),
        Vec3(0.0, 0.0, 1.0)
    };
    enclosed_hole_complex.triangles = {
        {0, 2, 1},
        {0, 1, 3},
        {1, 2, 3},
        {2, 0, 3}
    };
    detail::orient_and_label_sheet_components(
        enclosed_hole_complex
    );
    detail::classify_complex_topology(enclosed_hole_complex);
    std::vector<double> enclosed_stabilities{
        1.0 / 3.0,
        2.0 / 3.0,
        2.0 / 3.0,
        2.0 / 3.0
    };
    std::vector<bool> enclosed_keep{false, true, true, true};
    CrossResolutionStabilityOptions repair_options;
    repair_options.maximum_repaired_hole_area_fraction = 0.50;
    const detail::StabilityTopologyRepair enclosed_repair =
        detail::repair_stability_holes(
            enclosed_hole_complex,
            enclosed_stabilities,
            {1.0},
            enclosed_keep,
            repair_options
        );
    assert(enclosed_repair.repaired_patch_count == 1);
    assert(enclosed_repair.repaired[0]);
    assert(enclosed_keep[0]);

    MedialComplex boundary_gap_complex;
    boundary_gap_complex.vertices.assign(
        synthetic_complex.vertices.begin(),
        synthetic_complex.vertices.begin() + 4
    );
    boundary_gap_complex.triangles = {
        {0, 1, 2},
        {0, 2, 3}
    };
    detail::orient_and_label_sheet_components(boundary_gap_complex);
    detail::classify_complex_topology(boundary_gap_complex);
    std::vector<double> boundary_stabilities{
        1.0 / 3.0,
        2.0 / 3.0
    };
    std::vector<bool> boundary_keep{false, true};
    const detail::StabilityTopologyRepair boundary_repair =
        detail::repair_stability_holes(
            boundary_gap_complex,
            boundary_stabilities,
            {1.0 / 3.0},
            boundary_keep,
            repair_options
        );
    assert(boundary_repair.repaired_patch_count == 0);
    assert(!boundary_repair.repaired[0]);
    assert(!boundary_keep[0]);

    MedialComplex radius_hole_complex = enclosed_hole_complex;
    radius_hole_complex.vertex_radii.assign(4, 1.0);
    radius_hole_complex.triangle_radius_jumps = {
        0.9, 0.1, 0.1, 0.1
    };
    RadiusContinuityFilterOptions radius_hole_options;
    radius_hole_options.use_adaptive_lfs_thresholds = false;
    radius_hole_options.maximum_relative_radius_jump = 0.5;
    radius_hole_options.maximum_topology_repair_area_fraction = 0.50;
    const RadiusContinuityFilterResult radius_hole_repair =
        filter_medial_complex_radius_continuity(
            radius_hole_complex,
            radius_hole_options
        );
    assert(radius_hole_repair.removed_triangles.empty());
    assert(radius_hole_repair.retained.triangles.size() == 4);
    assert(radius_hole_repair.flagged_triangles.size() == 1);
    assert(radius_hole_repair.topology_repaired_triangles.empty());
    assert(radius_hole_repair.topology_repaired_patch_count == 0);

    MedialComplex candidate_gap_reference = enclosed_hole_complex;
    candidate_gap_reference.vertex_radii.assign(4, 1.0);
    candidate_gap_reference.triangle_radius_jumps.assign(4, 0.1);
    std::vector<bool> candidate_gap_keep{false, true, true, true};
    MedialComplexOptions candidate_gap_options;
    candidate_gap_options.minimum_gap_normal_alignment = 0.0;
    candidate_gap_options.maximum_gap_area_fraction = 0.50;
    const detail::CandidateGapRepair candidate_gap_repair =
        detail::restore_enclosed_candidate_gaps(
            candidate_gap_reference,
            candidate_gap_keep,
            candidate_gap_options
        );
    assert(candidate_gap_repair.repaired_patch_count == 1);
    assert(candidate_gap_repair.repaired[0]);
    assert(candidate_gap_keep[0]);

    std::vector<VoronoiFaceCandidate> polygon_gap_faces(4);
    polygon_gap_faces[0].source_tetrahedra = {0, 1, 2};
    polygon_gap_faces[1].source_tetrahedra = {1, 0, 3};
    polygon_gap_faces[2].source_tetrahedra = {2, 1, 4};
    polygon_gap_faces[3].source_tetrahedra = {0, 2, 5};
    const std::vector<bool> polygon_gap_eligible(4, true);
    const std::vector<bool> polygon_gap_inside(4, true);
    std::vector<bool> polygon_gap_keep{false, true, true, true};
    const detail::CandidatePolygonGapRepair polygon_gap_repair =
        detail::restore_enclosed_polygon_gaps(
            polygon_gap_faces,
            polygon_gap_eligible,
            polygon_gap_inside,
            polygon_gap_keep
        );
    assert(polygon_gap_repair.repaired_patch_count == 1);
    assert(polygon_gap_repair.repaired[0]);
    assert(polygon_gap_keep[0]);

    polygon_gap_keep = {false, true, true, false};
    const detail::CandidatePolygonGapRepair open_polygon_gap =
        detail::restore_enclosed_polygon_gaps(
            polygon_gap_faces,
            polygon_gap_eligible,
            polygon_gap_inside,
            polygon_gap_keep
        );
    assert(open_polygon_gap.repaired_patch_count == 0);
    assert(!polygon_gap_keep[0]);

    std::vector<VoronoiFaceCandidate> stratum_faces(3);
    stratum_faces[0].source_tetrahedra = {0, 1, 2};
    stratum_faces[1].source_tetrahedra = {1, 0, 3};
    stratum_faces[2].source_tetrahedra = {0, 3, 4};
    const std::vector<bool> stratum_eligible(3, true);
    const std::vector<bool> stratum_inside(3, true);
    std::vector<bool> stratum_keep{true, false, false};
    const detail::CandidatePolygonStratumCompletion stratum_completion =
        detail::complete_seeded_polygon_strata(
            stratum_faces,
            stratum_eligible,
            stratum_inside,
            stratum_keep
        );
    assert(stratum_completion.completed_stratum_count == 1);
    assert(!stratum_completion.repaired[0]);
    assert(stratum_completion.repaired[1]);
    assert(stratum_completion.repaired[2]);
    assert(std::all_of(
        stratum_keep.begin(),
        stratum_keep.end(),
        [](bool kept) { return kept; }
    ));

    std::vector<VoronoiFaceCandidate> seam_faces(3);
    seam_faces[0].source_tetrahedra = {0, 1, 2};
    seam_faces[1].source_tetrahedra = {1, 0, 3};
    seam_faces[2].source_tetrahedra = {0, 1, 4};
    std::vector<bool> seam_keep{true, false, false};
    const detail::CandidatePolygonStratumCompletion seam_completion =
        detail::complete_seeded_polygon_strata(
            seam_faces,
            stratum_eligible,
            stratum_inside,
            seam_keep
        );
    assert(seam_completion.completed_stratum_count == 0);
    assert(seam_keep[0]);
    assert(!seam_keep[1]);
    assert(!seam_keep[2]);

    MedialComplex termination_loop_complex;
    termination_loop_complex.vertices.assign(
        synthetic_complex.vertices.begin(),
        synthetic_complex.vertices.begin() + 4
    );
    termination_loop_complex.vertex_radii.assign(4, 1.0);
    termination_loop_complex.vertex_local_feature_sizes.assign(4, 1.0);
    termination_loop_complex.triangles = {
        {0, 1, 2},
        {0, 2, 3}
    };
    termination_loop_complex.triangle_contact_angles = {95.0, 95.0};
    detail::orient_and_label_sheet_components(termination_loop_complex);
    detail::classify_complex_topology(termination_loop_complex, 105.0);
    assert(termination_loop_complex.boundary_loops.size() == 1);
    assert(termination_loop_complex.boundary_loops[0].closed);
    assert(termination_loop_complex.boundary_loops[0].allowed);
    assert(termination_loop_complex.artificial_boundary_edges.empty());

    termination_loop_complex.triangle_contact_angles = {140.0, 140.0};
    detail::classify_complex_topology(termination_loop_complex, 105.0);
    assert(termination_loop_complex.boundary_loops.size() == 1);
    assert(!termination_loop_complex.boundary_loops[0].allowed);
    assert(termination_loop_complex.artificial_boundary_edges.size() == 4);

    PoleSelectionResult boundary_test_poles;
    MedialCandidate boundary_test_pole;
    boundary_test_pole.validated = true;
    boundary_test_pole.source_tetrahedron = 0;
    boundary_test_poles.poles.push_back(boundary_test_pole);
    termination_loop_complex.vertex_source_tetrahedra.assign(4, {0});
    termination_loop_complex.triangle_confidences.assign(2, 0.8);
    termination_loop_complex.triangle_radius_jumps.assign(2, 0.1);
    const auto artificial_boundary_metrics = analyze_medial_components(
        termination_loop_complex,
        boundary_test_poles
    );
    assert(artificial_boundary_metrics.size() == 1);
    assert(artificial_boundary_metrics[0].artificial_boundary_fraction >
           0.99);
    const MedialComplexFilterResult diagnostic_boundary_filter =
        filter_medial_complex_components(
            termination_loop_complex,
            boundary_test_poles
        );
    assert(diagnostic_boundary_filter.removed_triangles.empty());
    MedialComponentFilterOptions destructive_boundary_options;
    destructive_boundary_options.require_valid_sheet_boundaries = true;
    const MedialComplexFilterResult artificial_boundary_filtered =
        filter_medial_complex_components(
            termination_loop_complex,
            boundary_test_poles,
            destructive_boundary_options
        );
    assert(artificial_boundary_filtered.retained.triangles.empty());

    MedialComplex triangular_hole_reference = enclosed_hole_complex;
    triangular_hole_reference.vertex_radii.assign(4, 1.0);
    triangular_hole_reference.triangle_confidences.assign(4, 0.8);
    triangular_hole_reference.triangle_contact_angles.assign(4, 140.0);
    triangular_hole_reference.triangle_radius_jumps.assign(4, 0.1);
    triangular_hole_reference.
        triangle_radius_jump_thresholds.assign(4, 0.5);
    MedialComplex triangular_hole_input = triangular_hole_reference;
    triangular_hole_input.triangles.erase(
        triangular_hole_input.triangles.begin()
    );
    triangular_hole_input.triangle_confidences.erase(
        triangular_hole_input.triangle_confidences.begin()
    );
    triangular_hole_input.triangle_contact_angles.erase(
        triangular_hole_input.triangle_contact_angles.begin()
    );
    triangular_hole_input.triangle_radius_jumps.erase(
        triangular_hole_input.triangle_radius_jumps.begin()
    );
    triangular_hole_input.triangle_radius_jump_thresholds.erase(
        triangular_hole_input.
            triangle_radius_jump_thresholds.begin()
    );
    triangular_hole_input.triangle_components.clear();
    triangular_hole_input.boundary_edges.clear();
    triangular_hole_input.seam_edges.clear();
    triangular_hole_input.junction_vertices.clear();
    detail::orient_and_label_sheet_components(triangular_hole_input);
    detail::classify_complex_topology(triangular_hole_input);
    const TriangularHoleSealResult triangular_hole_result =
        seal_isolated_triangular_holes(
            triangular_hole_reference,
            triangular_hole_input
        );
    assert(triangular_hole_result.restored_triangles.size() == 1);
    assert(triangular_hole_result.sealed.triangles.size() == 4);
    assert(triangular_hole_result.sealed.boundary_edges.empty());

    std::cout << "Delaunay scaffold tests passed.\n";
    return 0;
}
