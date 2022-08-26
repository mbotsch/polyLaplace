//=============================================================================
// Copyright 2021 Hendrik Meyer, Astrid Bunge, Mario Botsch
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================

#include "VolumeSubdivision.h"
#include "../Volume/diffgeo_3D.h"

//=============================================================================

VolumeSubdivision::VolumeSubdivision(VolumeMesh &mesh) : mesh_(mesh) {}

void VolumeSubdivision::tetrahedra(unsigned int face_point,
                                   unsigned int cell_point)
{

    auto fp_prop =
        mesh_.request_face_property<VolumeMesh::PointT>("face points");
    auto cp_prop =
        mesh_.request_cell_property<VolumeMesh::PointT>("cell points");

    compute_3D_virtual_points_and_weights(mesh_, face_point, cell_point);
    // save current faces
    // we need to delete them later on
    auto f_it_pair = mesh_.faces();

    std::vector<VHandle> c_bary_vert;

    c_bary_vert.reserve(mesh_.n_vertices());

    std::vector<HFHandle> c_halffaces;
    std::vector<VHandle> v_handles;

    c_halffaces.reserve(4);
    v_handles.reserve(3);

    // save the new triangles that lie within the plane of the halfface
    auto hf_prop = mesh_.request_halfface_property<std::vector<HFHandle>>();
    // save barycenter vertex for faces and cells
    auto f_prop = mesh_.request_face_property<VHandle>();
    auto c_prop = mesh_.request_cell_property<VHandle>();

    // we precompute them so they are ordered
    for (auto f_it = mesh_.f_iter(); f_it.valid(); ++f_it)
    {
        f_prop[*f_it] = mesh_.add_vertex(fp_prop[*f_it]);
        //        f_prop[*f_it] = mesh_.add_vertex(mesh_.barycenter(*f_it));
    }

    for (auto c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
    {
//        Vec3d p(0, 0, 0);
//        unsigned int valence = 0;
//        for (auto cv_it = mesh_.cv_iter(*c_it); cv_it.valid();
//             ++cv_it, ++valence)
//        {
//            p += mesh_.vertex(*cv_it);
//        }
//        for (auto cf_it = mesh_.cf_iter(*c_it); cf_it.valid();
//             ++cf_it, ++valence)
//        {
//            p += mesh_.vertex(f_prop[*cf_it]);
//        }
//        p /= valence;
        //        c_prop[*c_it] = mesh_.add_vertex(p);
        c_prop[*c_it] = mesh_.add_vertex(cp_prop[*c_it]);
    }

    auto c_it_pair = mesh_.cells();
    for (auto c_it = c_it_pair.first; c_it < c_it_pair.second; ++c_it)
    {
        // get cell barycenter vertex
        auto c_bary = c_prop[*c_it];

        // save halfface in between halfegde and cell barycenter
        auto he_prop = mesh_.request_halfedge_property<HFHandle>();

        // create tetrahedrons for all halffaces of a cell
        for (auto c_hf_it = mesh_.chf_iter(*c_it); c_hf_it.valid(); ++c_hf_it)
        {

            // get face barycenter vertex
            auto f_bary = f_prop[mesh_.face_handle(*c_hf_it)];

            // save halffaces in between the vertex, face barycenter and cell barycenter
            auto v_prop = mesh_.request_vertex_property<HFHandle>();

            auto hf_hes_it_pair = mesh_.halfface_halfedges(*c_hf_it);
            OVM::FaceHandle face_handle;

            // halfface triangles for the current halfface
            std::vector<HFHandle> hf_hfs;

            // check if the halffacewas already split into the smaller triangles
            if (hf_prop[*c_hf_it].empty())
            {
                // create triangles constructed from one halfedge and the halfface barycenter
                for (auto hf_he_it = hf_hes_it_pair.first;
                     hf_he_it != hf_hes_it_pair.second; ++hf_he_it)
                {
                    auto he = mesh_.halfedge(*hf_he_it);
                    v_handles.clear();
                    v_handles.emplace_back(f_bary);
                    v_handles.emplace_back(he.from_vertex());
                    v_handles.emplace_back(he.to_vertex());
                    face_handle = mesh_.add_face(v_handles);
                    hf_hfs.emplace_back(mesh_.halfface_handle(face_handle, 0));
                    // add them to the opposite halfface
                    hf_prop[mesh_.opposite_halfface_handle(*c_hf_it)]
                        .emplace_back(mesh_.halfface_handle(face_handle, 1));
                }
                // reversing the vector because the edges of the opposite halface have a reversed order
                std::reverse(
                    hf_prop[mesh_.opposite_halfface_handle(*c_hf_it)].begin(),
                    hf_prop[mesh_.opposite_halfface_handle(*c_hf_it)].end());
            }
            else
            {
                // get triangles for this halfface
                hf_hfs = hf_prop[*c_hf_it];
            }

            int counter = 0;
            for (auto hf_he_it = hf_hes_it_pair.first;
                 hf_he_it != hf_hes_it_pair.second; ++hf_he_it, ++counter)
            {
                c_halffaces.clear();
                auto he = mesh_.halfedge(*hf_he_it);

                // check if halfface from halfedge to cell barycenter already exists
                if (!he_prop[*hf_he_it].is_valid())
                {
                    v_handles.clear();
                    v_handles.emplace_back(c_bary);
                    v_handles.emplace_back(he.from_vertex());
                    v_handles.emplace_back(he.to_vertex());
                    face_handle = mesh_.add_face(v_handles);

                    // add created halffaces to the corresponding halfedges
                    he_prop[*hf_he_it] = mesh_.halfface_handle(face_handle, 1);
                    he_prop[mesh_.opposite_halfedge_handle(*hf_he_it)] =
                        mesh_.halfface_handle(face_handle, 0);
                }
                c_halffaces.push_back(he_prop[*hf_he_it]);

                // check if halfface between from_vertex, halfface barycenter
                // and cell barycenter already exists
                if (!v_prop[he.from_vertex()].is_valid())
                {
                    v_handles.clear();
                    v_handles.emplace_back(f_bary);
                    v_handles.emplace_back(he.from_vertex());
                    v_handles.emplace_back(c_bary);
                    face_handle = mesh_.add_face(v_handles);

                    // add created halffaces to the corresponding vertices
                    c_halffaces.push_back(
                        mesh_.halfface_handle(face_handle, 1));
                    v_prop[he.from_vertex()] =
                        mesh_.halfface_handle(face_handle, 0);
                }
                else
                {
                    // get previously calculated halfface
                    c_halffaces.push_back(v_prop[he.from_vertex()]);
                }

                // check if halfface between from_vertex, halfface barycenter
                // and cell barycenter already exists
                if (!v_prop[he.to_vertex()].is_valid())
                {
                    v_handles.clear();
                    v_handles.emplace_back(f_bary);
                    v_handles.emplace_back(he.to_vertex());
                    v_handles.emplace_back(c_bary);
                    face_handle = mesh_.add_face(v_handles);

                    // add created halffaces to the corresponding vertices
                    c_halffaces.push_back(
                        mesh_.halfface_handle(face_handle, 0));
                    v_prop[he.to_vertex()] =
                        mesh_.halfface_handle(face_handle, 1);
                }
                else
                {
                    // get previously calculated halfface
                    c_halffaces.push_back(v_prop[he.to_vertex()]);
                }

                // add halfface between halfedge and face barycenter
                c_halffaces.emplace_back(hf_hfs[counter]);

                // create tetrahedron
                mesh_.add_cell(c_halffaces);
            }
        }
    }

    // derefence all faces
    // we can not directly delete them because then the indices would not match
    // the ones we had before adding cells
    for (auto f_it = f_it_pair.first; f_it < f_it_pair.second; ++f_it)
    {
        mesh_.delete_face(*f_it);
    }

    // delete all derferenced entities
    mesh_.collect_garbage();
}

void VolumeSubdivision::irregular_mesh(int n)
{
    auto divide_cell = mesh_.request_cell_property<bool>();

    int min = 0;
    int max = mesh_.n_cells();
    for (int i = 0; i < n && i < max; ++i)
    {
        int random = min + (rand() % static_cast<int>(max - min + 1));
        divide_cell[CHandle(random)] = true;
    }

    auto c_it_pair = mesh_.cells();
    for (auto c_it = c_it_pair.first; c_it < c_it_pair.second; ++c_it)
    {
        if (divide_cell[*c_it])
        {
            auto v = mesh_.add_vertex(mesh_.barycenter(*c_it));
            std::vector<OpenVolumeMesh::VertexHandle> vertices;
            std::vector<OpenVolumeMesh::HalfFaceHandle> halffaces;

            auto he_prop = mesh_.request_halfedge_property<HFHandle>();

            for (auto hf : mesh_.cell_halffaces(*c_it))
            {
                halffaces.clear();

                halffaces.emplace_back(hf);

                for (auto he : mesh_.halfface_halfedges(hf))
                {
                    if (!he_prop[he].is_valid())
                    {
                        vertices.clear();

                        auto e = mesh_.halfedge(he);
                        vertices.emplace_back(v);
                        vertices.emplace_back(e.from_vertex());
                        vertices.emplace_back(e.to_vertex());

                        auto f = mesh_.add_face(vertices);

                        halffaces.emplace_back(mesh_.halfface_handle(f, 1));
                        he_prop[mesh_.opposite_halfedge_handle(he)] =
                            mesh_.halfface_handle(f, 0);
                    }
                    else
                    {
                        halffaces.emplace_back(he_prop[he]);
                    }
                }

                mesh_.add_cell(halffaces);
            }
            mesh_.delete_cell(*c_it);
        }
    }
}

void VolumeSubdivision::full_truncation()
{

    auto edges = mesh_.edges();

    auto e_prop = mesh_.request_edge_property<VHandle>();

    // precompute edge middle vertices
    for (auto e : mesh_.edges())
    {
        e_prop[e] = mesh_.add_vertex(mesh_.barycenter(e));
    }

    std::vector<HFHandle> c_halffaces;
    std::vector<HFHandle> v_halffaces;
    std::vector<VHandle> c_vertices;
    std::vector<VHandle> v_vertices;

    OpenVolumeMesh::FaceHandle f;

    auto hf_prop = mesh_.request_halfface_property<HFHandle>();
    auto hf_prop2 = mesh_.request_halfface_property<std::vector<HFHandle>>();

    auto c_it_pair = mesh_.cells();
    for (auto c_it = c_it_pair.first; c_it < c_it_pair.second; ++c_it)
    {
        c_halffaces.clear();

        // create center cell
        for (auto hf : mesh_.cell_halffaces(*c_it))
        {
            if (!hf_prop[hf].is_valid())
            {
                c_vertices.clear();
                for (auto he : mesh_.halfface_halfedges(hf))
                {
                    c_vertices.push_back(e_prop[mesh_.edge_handle(he)]);
                }
                f = mesh_.add_face(c_vertices);
                hf_prop[hf] = mesh_.halfface_handle(f, 0);
                hf_prop[mesh_.opposite_halfface_handle(hf)] =
                    mesh_.halfface_handle(f, 1);
            }
            c_halffaces.push_back(hf_prop[hf]);
        }

        // create cell for each vertex in cell
        for (auto v : mesh_.cell_vertices(*c_it))
        {
            v_halffaces.clear();
            c_vertices.clear();

            HFHandle hf;
            OpenVolumeMesh::HalfEdgeHandle first_he;

            // find any halfedge that is incident to the cell and vertex
            for (auto v_he : mesh_.outgoing_halfedges(v))
            {
                for (auto c_he : mesh_.cell_halfedges(*c_it))
                {
                    if (v_he.idx() == c_he.idx())
                    {
                        first_he = v_he;
                    }
                }
            }

            // find any halfface that is incident to halfedge and cell
            for (auto he_hf : mesh_.halfedge_halffaces(first_he))
            {
                for (auto c_hf : mesh_.cell_halffaces(*c_it))
                {
                    if (he_hf.idx() == c_hf.idx())
                    {
                        hf = c_hf;
                    }
                }
            }

            auto cur_he = first_he;

            // iterate over faces incident to vertex and cell in order
            do
            {
                c_vertices.push_back(e_prop[mesh_.edge_handle(cur_he)]);
                v_vertices.clear();
                v_vertices.push_back(v);
                v_vertices.push_back(e_prop[mesh_.edge_handle(cur_he)]);

                for (auto he : mesh_.halfface_halfedges(hf))
                {
                    for (auto hee : mesh_.incoming_halfedges(v))
                    {
                        if (he.idx() == hee.idx() && he.idx() != cur_he.idx())
                        {
                            cur_he = mesh_.opposite_halfedge_handle(he);
                            break;
                        }
                    }
                }
                v_vertices.push_back(e_prop[mesh_.edge_handle(cur_he)]);

                auto hff = mesh_.halfface(v_vertices);
                if (!hff.is_valid())
                {
                    f = mesh_.add_face(v_vertices);
                    v_halffaces.push_back(mesh_.halfface_handle(f, 0));
                }
                else
                {
                    v_halffaces.push_back(hff);
                }

                hf = mesh_.adjacent_halfface_in_cell(hf, cur_he);
            } while (cur_he != first_he);

            f = mesh_.add_face(c_vertices);

            // TODO choose the corect side, but it seems to work fine for now
            c_halffaces.push_back(mesh_.halfface_handle(f, 0));
            v_halffaces.push_back(mesh_.halfface_handle(f, 1));

            mesh_.add_cell(v_halffaces);
        }

        mesh_.add_cell(c_halffaces);
    }

    for (auto e_it = edges.first; e_it < edges.second; ++e_it)
    {
        mesh_.delete_edge(*e_it);
    }

    mesh_.collect_garbage();
}

void VolumeSubdivision::quads()
{
    auto c_vertex = mesh_.request_cell_property<VHandle>();
    auto f_vertex = mesh_.request_face_property<VHandle>();
    auto e_vertex = mesh_.request_edge_property<VHandle>();

    auto vertices = mesh_.vertices();

    Vec3d v_pos;
    for (auto c : mesh_.cells())
    {
        v_pos = mesh_.barycenter(c);
        c_vertex[c] = mesh_.add_vertex(v_pos);
    }

    for (auto f : mesh_.faces())
    {
        v_pos = mesh_.barycenter(f);
        f_vertex[f] = mesh_.add_vertex(v_pos);
    }

    for (auto e : mesh_.edges())
    {
        v_pos = mesh_.barycenter(e);
        e_vertex[e] = mesh_.add_vertex(v_pos);
    }

    auto edges = mesh_.edges();

    auto cells = mesh_.cells();
    for (auto c_it = cells.first; c_it < cells.second; ++c_it)
    {
        for (auto v : mesh_.cell_vertices(*c_it))
        {

            HFHandle hf;
            OpenVolumeMesh::HalfEdgeHandle he;

            // find any halfedge that is incident to the cell and vertex
            for (auto v_he : mesh_.outgoing_halfedges(v))
            {
                for (auto c_he : mesh_.cell_halfedges(*c_it))
                {
                    if (v_he.idx() == c_he.idx())
                    {
                        he = v_he;
                    }
                }
            }

            // find any halfface that is incident to halfedge and cell
            for (auto he_hf : mesh_.halfedge_halffaces(he))
            {
                for (auto c_hf : mesh_.cell_halffaces(*c_it))
                {
                    if (he_hf.idx() == c_hf.idx())
                    {
                        hf = c_hf;
                    }
                }
            }

            auto cur_he = he;
            auto cur_hf = hf;

            OpenVolumeMesh::HalfEdgeHandle next_he;
            HFHandle next_hf;

            std::vector<VHandle> vertices;
            std::vector<HFHandle> halffaces;

            // iterate over faces ad edges incident to vertex and cell in order
            do
            {
                for (auto hee : mesh_.halfface_halfedges(cur_hf))
                {
                    for (auto heee : mesh_.incoming_halfedges(v))
                    {
                        if (hee.idx() == heee.idx() &&
                            hee.idx() != cur_he.idx())
                        {
                            next_he = mesh_.opposite_halfedge_handle(hee);
                            goto found;
                        }
                    }
                }

            found:

                vertices.clear();
                vertices.push_back(e_vertex[mesh_.edge_handle(cur_he)]);
                vertices.push_back(f_vertex[mesh_.face_handle(cur_hf)]);
                vertices.push_back(e_vertex[mesh_.edge_handle(next_he)]);
                vertices.push_back(v);

                auto hff = mesh_.halfface(vertices);
                if (!hff.is_valid())
                {
                    auto f = mesh_.add_face(vertices);
                    halffaces.push_back(mesh_.halfface_handle(f, 0));
                }
                else
                {
                    halffaces.push_back(hff);
                }

                next_hf = mesh_.adjacent_halfface_in_cell(cur_hf, next_he);

                vertices.clear();
                vertices.push_back(f_vertex[mesh_.face_handle(next_hf)]);
                vertices.push_back(e_vertex[mesh_.edge_handle(next_he)]);
                vertices.push_back(f_vertex[mesh_.face_handle(cur_hf)]);
                vertices.push_back(c_vertex[*c_it]);

                hff = mesh_.halfface(vertices);
                if (!hff.is_valid())
                {
                    auto f = mesh_.add_face(vertices);
                    halffaces.push_back(mesh_.halfface_handle(f, 0));
                }
                else
                {
                    halffaces.push_back(hff);
                }

                cur_he = next_he;
                cur_hf = next_hf;
            } while (cur_he != he);

            mesh_.add_cell(halffaces);
        }
    }

    for (auto e_it = edges.first; e_it < edges.second; ++e_it)
    {
        mesh_.delete_edge(*e_it);
    }

    mesh_.collect_garbage();
}

void VolumeSubdivision::linear_subdivision()
{
    auto c_vertex = mesh_.request_cell_property<VHandle>();
    auto e_vertex = mesh_.request_edge_property<VHandle>();

    auto vertices = mesh_.vertices();

    Vec3d v_pos;
    for (auto c : mesh_.cells())
    {
        v_pos = mesh_.barycenter(c);
        c_vertex[c] = mesh_.add_vertex(v_pos);
    }

    for (auto e : mesh_.edges())
    {
        v_pos = mesh_.barycenter(e);
        e_vertex[e] = mesh_.add_vertex(v_pos);
    }

    auto edges = mesh_.edges();

    auto cells = mesh_.cells();
    for (auto c_it = cells.first; c_it < cells.second; ++c_it)
    {
        for (auto v : mesh_.cell_vertices(*c_it))
        {

            HFHandle hf;
            OpenVolumeMesh::HalfEdgeHandle he;

            // find any halfedge that is incident to the cell and vertex
            for (auto v_he : mesh_.outgoing_halfedges(v))
            {
                for (auto c_he : mesh_.cell_halfedges(*c_it))
                {
                    if (v_he.idx() == c_he.idx())
                    {
                        he = v_he;
                    }
                }
            }

            // find any halfface that is incident to halfedge and cell
            for (auto he_hf : mesh_.halfedge_halffaces(he))
            {
                for (auto c_hf : mesh_.cell_halffaces(*c_it))
                {
                    if (he_hf.idx() == c_hf.idx())
                    {
                        hf = c_hf;
                    }
                }
            }

            auto cur_he = he;
            auto cur_hf = hf;

            OpenVolumeMesh::HalfEdgeHandle next_he;
            HFHandle next_hf;

            std::vector<VHandle> vertices;
            std::vector<HFHandle> halffaces;

            // iterate over faces ad edges incident to vertex and cell in order
            do
            {
                for (auto hee : mesh_.halfface_halfedges(cur_hf))
                {
                    for (auto heee : mesh_.incoming_halfedges(v))
                    {
                        if (hee.idx() == heee.idx() &&
                            hee.idx() != cur_he.idx())
                        {
                            next_he = mesh_.opposite_halfedge_handle(hee);
                            goto found;
                        }
                    }
                }

            found:

                vertices.clear();
                vertices.push_back(e_vertex[mesh_.edge_handle(cur_he)]);
                vertices.push_back(e_vertex[mesh_.edge_handle(next_he)]);
                vertices.push_back(v);

                auto hff = mesh_.halfface(vertices);
                if (!hff.is_valid())
                {
                    auto f = mesh_.add_face(vertices);
                    halffaces.push_back(mesh_.halfface_handle(f, 0));
                }
                else
                {
                    halffaces.push_back(hff);
                }

                vertices.clear();
                vertices.push_back(e_vertex[mesh_.edge_handle(cur_he)]);
                vertices.push_back(c_vertex[*c_it]);
                vertices.push_back(e_vertex[mesh_.edge_handle(next_he)]);

                next_hf = mesh_.adjacent_halfface_in_cell(cur_hf, next_he);

                hff = mesh_.halfface(vertices);
                if (!hff.is_valid())
                {
                    auto f = mesh_.add_face(vertices);
                    halffaces.push_back(mesh_.halfface_handle(f, 0));
                }
                else
                {
                    halffaces.push_back(hff);
                }

                cur_he = next_he;
                cur_hf = next_hf;
            } while (cur_he != he);

            mesh_.add_cell(halffaces);
        }

        // add cell for all faces of the original cell

        std::vector<VHandle> vertices;
        std::vector<VHandle> hf_vertices;
        std::vector<HFHandle> halffaces;

        for (auto hf : mesh_.cell_halffaces(*c_it))
        {
            halffaces.clear();
            hf_vertices.clear();

            auto hf_hes = mesh_.halfface_halfedges(hf);
            auto first_he = *hf_hes.first;
            auto cur_he = first_he;
            auto prev_he = mesh_.prev_halfedge_in_halfface(cur_he, hf);

            do
            {
                vertices.clear();
                vertices.push_back(e_vertex[mesh_.edge_handle(prev_he)]);
                vertices.push_back(c_vertex[*c_it]);
                vertices.push_back(e_vertex[mesh_.edge_handle(cur_he)]);

                auto hff = mesh_.halfface(vertices);
                if (!hff.is_valid())
                {
                    auto f = mesh_.add_face(vertices);
                    halffaces.push_back(mesh_.halfface_handle(f, 0));
                }
                else
                {
                    halffaces.push_back(hff);
                }

                hf_vertices.push_back(e_vertex[mesh_.edge_handle(cur_he)]);

                prev_he = cur_he;
                cur_he = mesh_.next_halfedge_in_halfface(cur_he, hf);
            } while (cur_he.idx() != first_he.idx());

            auto hff = mesh_.halfface(hf_vertices);
            if (!hff.is_valid())
            {
                auto f = mesh_.add_face(hf_vertices);
                halffaces.push_back(mesh_.halfface_handle(f, 0));
            }
            else
            {
                halffaces.push_back(hff);
            }

            mesh_.add_cell(halffaces);
        }
    }

    for (auto e_it = edges.first; e_it < edges.second; ++e_it)
    {
        mesh_.delete_edge(*e_it);
    }

    mesh_.collect_garbage();
}


