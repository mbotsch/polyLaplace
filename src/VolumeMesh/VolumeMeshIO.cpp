//=============================================================================
// Copyright 2021 Hendrik Meyer, Astrid Bunge, Mario Botsch
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================

#include "VolumeMeshIO.h"
#include <OpenVolumeMesh/FileManager/FileManager.hh>
#include <clocale>

//=============================================================================

typedef OpenVolumeMesh::FaceHandle FHandle;
typedef std::vector<VHandle> FaceTuple;
typedef std::map<FaceTuple, FHandle> FaceMap;

//=============================================================================

bool VolumeMeshIO::read(VolumeMesh &mesh) {
    std::setlocale(LC_NUMERIC, "C");

    // clear mesh before reading from file what does it to a volume mesh
    //mesh.clear();
    for (auto v: mesh.vertices()) {
        mesh.delete_vertex(v);
    }
    mesh.collect_garbage();

    // extract file extension
    std::string::size_type dot(filename_.rfind('.'));
    if (dot == std::string::npos)
        return false;
    std::string ext = filename_.substr(dot + 1, filename_.length() - dot - 1);
    std::transform(ext.begin(), ext.end(), ext.begin(), tolower);

    // extension determines reader
    if (ext == "ovm") {
        return read_ovm(mesh);
    } else if (ext == "mesh") {
        std::cout << "test" << std::endl;
        return read_mesh(mesh);
    } else if (ext == "hybrid") {
        return read_hybrid(mesh);
    }

    return false;
}
//-----------------------------------------------------------------------------

bool VolumeMeshIO::write_histrogram() {

    std::setlocale(LC_NUMERIC, "C");

    // extract file extension
    std::string::size_type dot(filename_.rfind('.'));
    if (dot == std::string::npos)
        return false;
    std::string ext = filename_.substr(dot + 1, filename_.length() - dot - 1);
    std::transform(ext.begin(), ext.end(), ext.begin(), tolower);

    // extension determines reader
    if (ext == "hybrid") {

        typedef std::map<int, int> MapType;
        MapType hist;

        auto iff = std::ifstream(filename_, std::ios::in);

        if (!iff.good()) {
            std::cerr << "Could not open file " << filename_ << " for reading!"
                      << std::endl;
            return false;
        }

        std::string next;
        int c = 0;
        int nv, nf, np;
//
        while (iff.good()) {
            if (c == 0) {
                iff >> nv >> nf >> np;
                np = 2 * np / 3;
                std::cout << nv << ", " << nf << ", " << np << std::endl;
            }

            if (c > nv + nf && c <= nv + nf + np) {
                int val;
                iff >> val;
                std::cout << val << std::endl;
                if ((c - nv + nf) % 2 == 0) {
                    // Histogramm of cell valences
                    int k = val;
                    int v = 1;

                    auto lb = hist.lower_bound(k);

                    if (lb != hist.end() && !(hist.key_comp()(k, lb->first))) {
                        // key already exists
                        // update lb->second if you care to
                        v = lb->second;
                        lb->second = v + 1;
                    } else {
                        // the key does not exist in the map
                        // add it to the map
                        hist.insert(lb, MapType::value_type(k, v));    // Use lb as a hint to insert,
                        // so it can avoid another lookup
                    }
                }
            }
            c++;
        }
        iff.close();

        std::string::size_type slash(filename_.rfind('/'));
        if (slash == std::string::npos)
            return false;
        std::string hist_name = filename_.substr(slash + 1, dot - 1) + "_histogram.txt";
        std::ofstream file_(hist_name);
        MapType::iterator it;
        for (it = hist.begin(); it != hist.end(); it++) {
            file_ << it->first    // string (key)
                  << ':'
                  << it->second   // string's value
                  << std::endl;
        }
        file_.close();
    }
    return true;
}
//-----------------------------------------------------------------------------

bool VolumeMeshIO::write(const VolumeMesh &mesh) {
    // extract file extension
    std::string::size_type dot(filename_.rfind('.'));
    if (dot == std::string::npos)
        return false;
    std::string ext = filename_.substr(dot + 1, filename_.length() - dot - 1);
    std::transform(ext.begin(), ext.end(), ext.begin(), tolower);

    // extension determines reader
    if (ext == "ovm") {
        return write_ovm(mesh);
    } else if (ext == "mesh") {
//        return write_mesh(mesh);
        std::cout << "Write .mesh currently not supported " << std::endl;
        return false;
    } else if (ext == "hybrid") {
        std::cout << "Wrote Hyprid file " << std::endl;
        return write_hybrid(mesh);
    }

    return false;
}

//-----------------------------------------------------------------------------
bool VolumeMeshIO::read_hybrid(VolumeMesh &mesh) {
    auto iff = std::ifstream(filename_, std::ios::in);
    if (!iff.good()) {
        std::cerr << "Could not open file " << filename_ << " for reading!"
                  << std::endl;
        return false;
    }

    std::string next;
    int c = 0;
    int cell_c = 0;
    int nv, nf, np;
    Vec3d vec;
    std::vector<FHandle> faces;
    std::vector<std::vector<FHandle>> cells;
    std::vector<std::vector<int>> cell_orientation;

    while (iff.good()) {
        if (c == 0) {
            iff >> nv >> nf >> np;
            std::cout << "Nr polyhedra: " << 2 * np / 3 << std::endl;
            np = 2 * np / 3;
        } else if (c <= nv) {
            iff >> vec[0] >> vec[1] >> vec[2];
            mesh.add_vertex(vec);
        } else if (c > nv && c <= nv + nf) {
            int val;
            iff >> val;
            std::vector<VHandle> v_vec;
            v_vec.resize(val);
            for (int i = 0; i < val; i++) {
                iff >> v_vec[i];
            }
            FHandle fh = mesh.add_face(v_vec);
            faces.emplace_back(fh);
        } else if (c > nv + nf && c <= nv + nf + np) {
            int val;
            iff >> val;
            std::vector<FHandle> cell_faces;
            std::vector<int> orientation;
            std::vector<HFHandle> cell_hfaces;
            if ((c - nv + nf) % 2 == 1) {
                int idx;
                for (int i = 0; i < val; i++) {
                    iff >> idx;
                    if (idx > (int) faces.size()) {
                        std::cout << "faces out of range" << std::endl;
                    }
                    if (idx == -1) {
                        std::cout << "faces idx not valid" << std::endl;
                    }
                    cell_faces.emplace_back(faces[idx]);
                }
                cells.emplace_back(cell_faces);
            } else {
                int sign;
                if (cell_c > (int) cells.size()) {
                    std::cout << "cells out of range" << std::endl;
                }
                cell_faces = cells[cell_c];
                for (int i = 0; i < val; i++) {
                    iff >> sign;
                    if (!mesh.face_halffaces(cell_faces[i])[sign].is_valid()) {
                        std::cout << "halfface idx is a problem" << std::endl;
                    }
                    cell_hfaces.emplace_back(
                            mesh.face_halffaces(cell_faces[i])[sign]);
                }
                CHandle ch = mesh.add_cell(cell_hfaces);
                if (!ch.is_valid()) {
                    std::cout << "cell handle problem" << std::endl;
                }
                cell_c++;
                //                std::cout << "valence: " << val << std::endl;
                std::cout << "cell nr: " << cell_c << std::endl;
                std::cout << "--------------------" << std::endl;
            }
        } else {
            std::cout << "c: " << c << " nv+nf+np: " << nv + nf + np
                      << std::endl;
            //            continue;
            break;
        }
        c++;
    }
    iff.close();

    std::cout << "Converted " << mesh.n_vertices() << " vertices," << std::endl;
    std::cout << "\t  " << mesh.n_edges() << " edges," << std::endl;
    std::cout << "\t  " << mesh.n_faces() << " faces," << std::endl;
    std::cout << "\t  " << mesh.n_cells() << " cells!" << std::endl;

    return true;
}

//-----------------------------------------------------------------------------

bool VolumeMeshIO::read_ovm(VolumeMesh &mesh) {
    OpenVolumeMesh::IO::FileManager file_manager;
    return file_manager.readFile(filename_, mesh);
}

//-----------------------------------------------------------------------------

bool VolumeMeshIO::write_ovm(const VolumeMesh &mesh) {
    OpenVolumeMesh::IO::FileManager file_manager;
    return file_manager.writeFile(filename_, mesh);
}

//-----------------------------------------------------------------------------

void read_mesh_vertices(std::istream &_istream, VolumeMesh &_mesh) {
    int n_vertices;
    _istream >> n_vertices;

    std::string coord;
    Vec3d vec;
    double x;

    for (int i = 0; i < n_vertices; i++) {
        _istream >> vec[0] >> vec[1] >> vec[2] >> x;

        _mesh.add_vertex(vec);
    }
}

//-----------------------------------------------------------------------------

void read_mesh_tetrahedra(std::istream &_istream, VolumeMesh &_mesh) {
    int n_tetrahedra;
    _istream >> n_tetrahedra;

    std::vector<int> idx;
    idx.resize(4);
    float x;

    std::vector<VHandle> c_vertices_;
    std::vector<FaceTuple> tuples;
    std::vector<HFHandle> cell_halffaces;

    FaceMap face_map;

    for (int i = 0; i < n_tetrahedra; i++) {
        idx.clear();
        _istream >> idx[0] >> idx[1] >> idx[2] >> idx[3] >> x;

        c_vertices_.clear();
        c_vertices_.emplace_back(VHandle(idx[0] - 1));
        c_vertices_.emplace_back(VHandle(idx[1] - 1));
        c_vertices_.emplace_back(VHandle(idx[2] - 1));
        c_vertices_.emplace_back(VHandle(idx[3] - 1));

        Vec3d midP(0.0, 0.0, 0.0);
        double valence = 0.0;
        for (auto c_vertex: c_vertices_) {
            midP += _mesh.vertex(c_vertex);
            valence += 1.0;
        }
        midP /= valence;

        std::sort(c_vertices_.begin(), c_vertices_.end());

        tuples.clear();
        tuples.emplace_back(
                FaceTuple{c_vertices_[0], c_vertices_[1], c_vertices_[2]});
        tuples.emplace_back(
                FaceTuple{c_vertices_[1], c_vertices_[2], c_vertices_[3]});
        tuples.emplace_back(
                FaceTuple{c_vertices_[0], c_vertices_[2], c_vertices_[3]});
        tuples.emplace_back(
                FaceTuple{c_vertices_[0], c_vertices_[1], c_vertices_[3]});

        cell_halffaces.clear();

        for (const auto &it: tuples) {

            // Check if face exists for current tuple
            auto f = face_map.find(it);
            if (f == face_map.end()) {
                // Face does not exist, create it

                // Find right orientation, s.t. normal
                // points inside the cell

                Vec3d e1 = _mesh.vertex(it[1]) - _mesh.vertex(it[0]);
                Vec3d e2 = _mesh.vertex(it[2]) - _mesh.vertex(it[1]);

                // Get face normal (cross product)
                Vec3d n = (e1 % e2).normalize();

                std::vector<VHandle> v_vec;
                v_vec.push_back(it[0]);
                v_vec.push_back(it[1]);
                v_vec.push_back(it[2]);
                FHandle fh = _mesh.add_face(v_vec);

                // Add face to face map
                face_map[it] = fh;

                // Check whether normal points inside cell
                if (((midP - _mesh.vertex(it[0])) | n) > 0.0) {

                    // Normal points inside cell, just add half-face 0
                    // Add corresponding half-face to cell definition
                    cell_halffaces.push_back(VolumeMesh::halfface_handle(fh, 0));
                } else {

                    // Normal points outside cell, just add half-face 1
                    // Add corresponding half-face to cell definition
                    cell_halffaces.push_back(VolumeMesh::halfface_handle(fh, 1));
                }
            } else {

                // Face exists, find right orientation
                FHandle fh = f->second;

                std::vector<OpenVolumeMesh::HalfEdgeHandle> hes =
                        _mesh.face(fh).halfedges();

                assert(hes.size() == 3);

                Vec3d e1 = _mesh.vertex(_mesh.halfedge(hes[0]).to_vertex()) -
                           _mesh.vertex(_mesh.halfedge(hes[0]).from_vertex());
                Vec3d e2 = _mesh.vertex(_mesh.halfedge(hes[1]).to_vertex()) -
                           _mesh.vertex(_mesh.halfedge(hes[1]).from_vertex());

                Vec3d n = (e1 % e2).normalize();

                if (((midP -
                      _mesh.vertex(_mesh.halfedge(hes[0]).from_vertex())) |
                     n) > 0.0) {
                    // Normal points inside cell
                    cell_halffaces.push_back(VolumeMesh::halfface_handle(fh, 0));
                } else {
                    // Normal points outisde cell
                    cell_halffaces.push_back(VolumeMesh::halfface_handle(fh, 1));
                }
            }
        }
        _mesh.add_cell(cell_halffaces);
    }
}

void read_mesh_hexahedra(std::istream &_istream, VolumeMesh &_mesh) {
    int n_hexahedra;
    _istream >> n_hexahedra;

    std::vector<int> idx;
    idx.resize(8);
    float x;

    std::vector<VHandle> c_vertices_;
    std::vector<FaceTuple> tuples;
    std::vector<HFHandle> cell_halffaces;

    FaceMap face_map;

    for (int i = 0; i < n_hexahedra; i++) {
        idx.clear();
        _istream >> idx[0] >> idx[1] >> idx[2] >> idx[3] >> idx[4] >> idx[5] >>
                 idx[6] >> idx[7] >> x;

        c_vertices_.clear();
        for (int j = 0; j < 8; j++) {
            c_vertices_.emplace_back(VHandle(idx[j] - 1));
        }

        Vec3d midP(0.0, 0.0, 0.0);
        double valence = 0.0;
        for (auto c_vertex: c_vertices_) {
            midP += _mesh.vertex(c_vertex);
            valence += 1.0;
        }
        midP /= valence;

        int hexahedron_indices_[6][4] = {{0, 1, 2, 3},
                                         {5, 4, 7, 6},
                                         {4, 0, 3, 7},
                                         {1, 5, 6, 2},
                                         {4, 5, 1, 0},
                                         {6, 7, 3, 2}};

        tuples.clear();
        FaceTuple tuple;
        for (auto face: hexahedron_indices_) {
            tuple.clear();
            for (size_t k = 0; k < 4; ++k) {
                tuple.emplace_back(c_vertices_[face[k]]);
            }
            tuples.emplace_back(tuple);
        }

        cell_halffaces.clear();

        for (const auto &it: tuples) {
            std::vector<VHandle> key(it);
            std::sort(key.begin(), key.end());
            auto f = face_map.find(key);
            if (f == face_map.end()) {

                Vec3d e1 = _mesh.vertex(it[1]) - _mesh.vertex(it[0]);
                Vec3d e2 = _mesh.vertex(it[2]) - _mesh.vertex(it[1]);

                Vec3d n = (e1 % e2).normalize();

                FHandle fh = _mesh.add_face(it);

                face_map[key] = fh;

                if (((midP - _mesh.vertex(it[0])) | n) > 0.0) {
                    cell_halffaces.push_back(VolumeMesh::halfface_handle(fh, 0));
                } else {
                    cell_halffaces.push_back(VolumeMesh::halfface_handle(fh, 1));
                }
            } else {
                FHandle fh = f->second;

                std::vector<OpenVolumeMesh::HalfEdgeHandle> hes =
                        _mesh.face(fh).halfedges();

                assert(hes.size() == 4);

                Vec3d e1 = _mesh.vertex(_mesh.halfedge(hes[0]).to_vertex()) -
                           _mesh.vertex(_mesh.halfedge(hes[0]).from_vertex());
                Vec3d e2 = _mesh.vertex(_mesh.halfedge(hes[1]).to_vertex()) -
                           _mesh.vertex(_mesh.halfedge(hes[1]).from_vertex());

                Vec3d n = (e1 % e2).normalize();

                if (((midP -
                      _mesh.vertex(_mesh.halfedge(hes[0]).from_vertex())) |
                     n) > 0.0) {
                    // Normal points inside cell
                    cell_halffaces.push_back(VolumeMesh::halfface_handle(fh, 0));
                } else {
                    // Normal points outisde cell
                    cell_halffaces.push_back(VolumeMesh::halfface_handle(fh, 1));
                }
            }
        }
        _mesh.add_cell(cell_halffaces);
    }
}

//-----------------------------------------------------------------------------

bool VolumeMeshIO::read_mesh(VolumeMesh &mesh) {
    auto iff = std::ifstream(filename_, std::ios::in);

    if (!iff.good()) {
        std::cerr << "Could not open file " << filename_ << " for reading!"
                  << std::endl;
        return false;
    }

    std::string next;
    float value;
    while (iff.good()) {
        iff >> next;
        if (next == "MeshVersionFormatted") {
            iff >> value;
        } else if (next == "Dimension") {
            iff >> value;
        } else if (next[0] == '#') {
            getline(iff, next);
        } else if (next == "Vertices") {
            read_mesh_vertices(iff, mesh);
        } else if (next == "Tetrahedra") {
            read_mesh_tetrahedra(iff, mesh);
        } else if (next == "Hexahedra") {
            read_mesh_hexahedra(iff, mesh);
        } else if (next == "End") {
            break;
        } else {
            continue;
        }
    }
    iff.close();

    std::cout << "Converted " << mesh.n_vertices() << " vertices," << std::endl;
    std::cout << "\t  " << mesh.n_edges() << " edges," << std::endl;
    std::cout << "\t  " << mesh.n_faces() << " faces," << std::endl;
    std::cout << "\t  " << mesh.n_cells() << " cells!" << std::endl;

    return true;
}

//-----------------------------------------------------------------------------

//bool VolumeMeshIO::write_mesh(const VolumeMesh &mesh)
//{
//    std::cout << "Currently not supported." << std::endl;
//    return false;
//}

bool VolumeMeshIO::write_hybrid(const VolumeMesh &_mesh) {
    std::ofstream _ostream(filename_.c_str(), std::ios::out);
    if (!_ostream.good()) {
        return false;
    }
    // Write header

    uint64_t n_vertices(_mesh.n_vertices());
    uint64_t n_faces(_mesh.n_faces());
    uint64_t n_cells(_mesh.n_cells());
    _ostream << n_vertices << " " << n_faces << " " << 3 * n_cells << std::endl;

    typedef typename VolumeMesh::PointT Point;

    // write vertices
    for (auto v_it = _mesh.v_iter(); v_it; ++v_it) {

        Point v = _mesh.vertex(*v_it);
        _ostream << v[0] << " " << v[1] << " " << v[2] << std::endl;
    }

    // write faces
    for (auto f: _mesh.faces()) {

        _ostream << static_cast<uint64_t>(_mesh.valence(f)) << " ";
        for (auto v: _mesh.face_vertices(f)) {
            _ostream << v.idx() << " ";
        }
        _ostream << std::endl;
    }

    for (auto c: _mesh.cells()) {

        _ostream << static_cast<uint64_t>(_mesh.valence(c)) << " ";
        std::vector<int> orientation;
        std::vector<FHandle> faces;
        for (auto f: _mesh.cell_faces(c)) {
            faces.emplace_back(f);
        }
        int ctr = 0;
        for (auto h: _mesh.cell_halffaces(c)) {
            FHandle f = faces[ctr];
            _ostream << f.idx() << " ";
            if (h == _mesh.face_halffaces(f)[0]) {
                orientation.emplace_back(0);
            } else {
                orientation.emplace_back(1);
            }
            ctr++;
        }
        _ostream << std::endl;
        _ostream << static_cast<uint64_t>(_mesh.valence(c)) << " ";
        for (int o: orientation) {
            _ostream << o << " ";
        }
        _ostream << std::endl;
    }
    for (auto c: _mesh.cells()) {
        if (_mesh.valence(c) == 6) {
            _ostream << static_cast<uint64_t>(1) << std::endl;
        } else {
            _ostream << static_cast<uint64_t>(0) << std::endl;
        }

    }

    return _ostream.good();

}
