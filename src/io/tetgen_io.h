#pragma once

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <istream>
#include <limits>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "../core/delaunay3.h"
#include "../core/mesh.h"
#include "../core/vec3.h"

namespace medial_axis_3d {

struct TetGenNodeData {
    std::vector<int> ids;
    std::vector<Vec3> points;
    std::vector<std::vector<double>> attributes;
    std::vector<int> boundary_markers;
    int index_base{1};
};

struct TetGenFaceData {
    std::vector<int> ids;
    std::vector<std::array<int, 3>> node_ids;
    std::vector<int> boundary_markers;
    int index_base{1};
};

namespace detail {

inline bool next_data_line(std::istream& input, std::string& line) {
    while (std::getline(input, line)) {
        const std::size_t comment = line.find('#');
        if (comment != std::string::npos) {
            line.erase(comment);
        }
        if (line.find_first_not_of(" \t\r\n") != std::string::npos) {
            return true;
        }
    }
    return false;
}

}  // namespace detail

inline bool read_tetgen_node(std::istream& input, TetGenNodeData& data, std::string& error) {
    data = TetGenNodeData{};
    error.clear();

    std::string line;
    if (!detail::next_data_line(input, line)) {
        error = "The .node file is empty.";
        return false;
    }

    std::size_t point_count = 0;
    int dimension = 0;
    int attribute_count = 0;
    int marker_count = 0;
    {
        std::istringstream header(line);
        if (!(header >> point_count >> dimension >> attribute_count >> marker_count)) {
            error = "Invalid .node header; expected: <point count> <dimension> <attribute count> <marker count>.";
            return false;
        }
    }

    if (dimension != 3) {
        error = "This project requires a 3D TetGen .node file (header dimension must be 3).";
        return false;
    }
    if (attribute_count < 0 || (marker_count != 0 && marker_count != 1)) {
        error = "Invalid .node header attribute or boundary-marker count.";
        return false;
    }
    if (point_count < 4) {
        error = "At least four non-coplanar points are required for a 3D tetrahedralization.";
        return false;
    }

    data.ids.reserve(point_count);
    data.points.reserve(point_count);
    data.attributes.reserve(point_count);
    if (marker_count == 1) {
        data.boundary_markers.reserve(point_count);
    }

    std::set<int> seen_ids;
    for (std::size_t row = 0; row < point_count; ++row) {
        if (!detail::next_data_line(input, line)) {
            error = "The .node file ended before all declared points were read.";
            return false;
        }

        std::istringstream record(line);
        int id = 0;
        Vec3 point{};
        if (!(record >> id >> point.x >> point.y >> point.z)) {
            error = "Invalid point record at .node data row " + std::to_string(row + 1) + ".";
            return false;
        }
        if (!seen_ids.insert(id).second) {
            error = "Duplicate point id " + std::to_string(id) + " in .node file.";
            return false;
        }

        std::vector<double> attributes(static_cast<std::size_t>(attribute_count));
        for (double& attribute : attributes) {
            if (!(record >> attribute)) {
                error = "Missing attribute in .node point id " + std::to_string(id) + ".";
                return false;
            }
        }

        int boundary_marker = 0;
        if (marker_count == 1 && !(record >> boundary_marker)) {
            error = "Missing boundary marker in .node point id " + std::to_string(id) + ".";
            return false;
        }

        std::string extra;
        if (record >> extra) {
            error = "Unexpected extra field in .node point id " + std::to_string(id) + ".";
            return false;
        }

        data.ids.push_back(id);
        data.points.push_back(point);
        data.attributes.push_back(std::move(attributes));
        if (marker_count == 1) {
            data.boundary_markers.push_back(boundary_marker);
        }
    }

    int min_id = data.ids.front();
    int max_id = data.ids.front();
    for (int id : data.ids) {
        min_id = std::min(min_id, id);
        max_id = std::max(max_id, id);
    }
    if ((min_id != 0 && min_id != 1) ||
        max_id != min_id + static_cast<int>(point_count) - 1) {
        error = "TetGen point ids must be consecutive and start at 0 or 1.";
        return false;
    }
    data.index_base = min_id;

    return true;
}

inline bool read_tetgen_node(const std::filesystem::path& path, TetGenNodeData& data, std::string& error) {
    std::ifstream input(path);
    if (!input) {
        error = "Could not open .node file: " + path.string();
        return false;
    }
    return read_tetgen_node(input, data, error);
}

inline bool read_tetgen_face(std::istream& input, TetGenFaceData& data, std::string& error) {
    data = TetGenFaceData{};
    error.clear();

    std::string line;
    if (!detail::next_data_line(input, line)) {
        error = "The .face file is empty.";
        return false;
    }

    std::size_t face_count = 0;
    int marker_count = 0;
    {
        std::istringstream header(line);
        if (!(header >> face_count >> marker_count)) {
            error = "Invalid .face header; expected: <face count> <marker count>.";
            return false;
        }
    }
    if (face_count < 4) {
        error = "A closed 3D surface requires at least four faces.";
        return false;
    }
    if (marker_count != 0 && marker_count != 1) {
        error = "The .face boundary-marker count must be 0 or 1.";
        return false;
    }

    data.ids.reserve(face_count);
    data.node_ids.reserve(face_count);
    if (marker_count == 1) {
        data.boundary_markers.reserve(face_count);
    }

    std::set<int> seen_ids;
    for (std::size_t row = 0; row < face_count; ++row) {
        if (!detail::next_data_line(input, line)) {
            error = "The .face file ended before all declared faces were read.";
            return false;
        }

        std::istringstream record(line);
        int id = 0;
        std::array<int, 3> nodes{};
        if (!(record >> id >> nodes[0] >> nodes[1] >> nodes[2])) {
            error = "Invalid triangle record at .face data row " +
                    std::to_string(row + 1) + ".";
            return false;
        }
        if (!seen_ids.insert(id).second) {
            error = "Duplicate face id " + std::to_string(id) + " in .face file.";
            return false;
        }
        if (nodes[0] == nodes[1] ||
            nodes[1] == nodes[2] ||
            nodes[2] == nodes[0]) {
            error = "Face id " + std::to_string(id) + " repeats a node id.";
            return false;
        }

        int boundary_marker = 0;
        if (marker_count == 1 && !(record >> boundary_marker)) {
            error = "Missing boundary marker in .face id " +
                    std::to_string(id) + ".";
            return false;
        }
        std::string extra;
        if (record >> extra) {
            error = "Unexpected extra field in .face id " +
                    std::to_string(id) + ".";
            return false;
        }

        data.ids.push_back(id);
        data.node_ids.push_back(nodes);
        if (marker_count == 1) {
            data.boundary_markers.push_back(boundary_marker);
        }
    }

    const auto [minimum, maximum] =
        std::minmax_element(data.ids.begin(), data.ids.end());
    if ((*minimum != 0 && *minimum != 1) ||
        *maximum != *minimum + static_cast<int>(face_count) - 1) {
        error = "TetGen face ids must be consecutive and start at 0 or 1.";
        return false;
    }
    data.index_base = *minimum;
    return true;
}

inline bool read_tetgen_face(const std::filesystem::path& path,
                             TetGenFaceData& data,
                             std::string& error) {
    std::ifstream input(path);
    if (!input) {
        error = "Could not open .face file: " + path.string();
        return false;
    }
    return read_tetgen_face(input, data, error);
}

inline bool build_surface_mesh(const TetGenNodeData& nodes,
                               const TetGenFaceData& faces,
                               Mesh& mesh,
                               std::string& error) {
    mesh = Mesh{};
    error.clear();

    std::map<int, int> vertex_for_node_id;
    for (std::size_t i = 0; i < nodes.ids.size(); ++i) {
        vertex_for_node_id[nodes.ids[i]] = static_cast<int>(i);
    }

    mesh.vertices = nodes.points;
    mesh.faces.reserve(faces.node_ids.size());
    mesh.face_boundary_markers.reserve(faces.node_ids.size());
    for (std::size_t face_index = 0;
         face_index < faces.node_ids.size();
         ++face_index) {
        Triangle triangle;
        for (int corner = 0; corner < 3; ++corner) {
            const int node_id =
                faces.node_ids[face_index][static_cast<std::size_t>(corner)];
            const auto vertex = vertex_for_node_id.find(node_id);
            if (vertex == vertex_for_node_id.end()) {
                error = "Face id " + std::to_string(faces.ids[face_index]) +
                        " references missing node id " + std::to_string(node_id) + ".";
                return false;
            }
            triangle.vertices[static_cast<std::size_t>(corner)] = vertex->second;
        }
        mesh.faces.push_back(triangle);
        mesh.face_boundary_markers.push_back(
            face_index < faces.boundary_markers.size()
                ? faces.boundary_markers[face_index]
                : 0
        );
    }

    return orient_and_analyze_closed_mesh(mesh, error);
}

inline bool write_tetgen_node(std::ostream& output,
                              const std::vector<Vec3>& points,
                              const std::vector<int>& boundary_markers,
                              int index_base,
                              std::string& error) {
    error.clear();
    const bool write_markers = boundary_markers.size() == points.size();
    if (!boundary_markers.empty() && !write_markers) {
        error = "Cannot write .node: boundary-marker count does not match point count.";
        return false;
    }

    output << std::setprecision(std::numeric_limits<double>::max_digits10);
    output << points.size() << " 3 0 " << (write_markers ? 1 : 0) << '\n';
    for (std::size_t i = 0; i < points.size(); ++i) {
        const Vec3& point = points[i];
        output << (index_base + static_cast<int>(i)) << ' '
               << point.x << ' ' << point.y << ' ' << point.z;
        if (write_markers) {
            output << ' ' << boundary_markers[i];
        }
        output << '\n';
    }
    if (!output) {
        error = "Failed while writing .node output.";
        return false;
    }
    return true;
}

inline bool write_tetgen_node(const std::filesystem::path& path,
                              const std::vector<Vec3>& points,
                              const std::vector<int>& boundary_markers,
                              int index_base,
                              std::string& error) {
    std::ofstream output(path);
    if (!output) {
        error = "Could not open .node output file: " + path.string();
        return false;
    }
    return write_tetgen_node(
        output,
        points,
        boundary_markers,
        index_base,
        error
    );
}

inline bool write_tetgen_ele(std::ostream& output,
                             const Delaunay3& delaunay,
                             const std::vector<int>& node_ids,
                             int element_index_base,
                             std::string& error) {
    error.clear();
    if (node_ids.size() != delaunay.point_count()) {
        error = "Cannot write .ele: node id count does not match the tetrahedralization point count.";
        return false;
    }

    output << delaunay.tetrahedron_count() << " 4 0\n";
    for (std::size_t i = 0; i < delaunay.tetrahedron_count(); ++i) {
        const auto& tetrahedron = delaunay.tetrahedra()[i];
        output << (element_index_base + static_cast<int>(i));
        for (int vertex : tetrahedron.vertices) {
            if (vertex < 0 || vertex >= static_cast<int>(node_ids.size())) {
                error = "Cannot write .ele: tetrahedron contains an invalid vertex index.";
                return false;
            }
            output << ' ' << node_ids[static_cast<std::size_t>(vertex)];
        }
        output << '\n';
    }

    if (!output) {
        error = "Failed while writing .ele output.";
        return false;
    }
    return true;
}

inline bool write_tetgen_ele(const std::filesystem::path& path,
                             const Delaunay3& delaunay,
                             const std::vector<int>& node_ids,
                             int element_index_base,
                             std::string& error) {
    std::ofstream output(path);
    if (!output) {
        error = "Could not open .ele output file: " + path.string();
        return false;
    }
    return write_tetgen_ele(output, delaunay, node_ids, element_index_base, error);
}

}  // namespace medial_axis_3d
