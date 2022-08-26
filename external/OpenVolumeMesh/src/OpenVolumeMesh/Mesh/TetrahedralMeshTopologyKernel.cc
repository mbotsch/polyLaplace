/*===========================================================================*\
 *                                                                           *
 *                            OpenVolumeMesh                                 *
 *        Copyright (C) 2011 by Computer Graphics Group, RWTH Aachen         *
 *                        www.openvolumemesh.org                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *  This file is part of OpenVolumeMesh.                                     *
 *                                                                           *
 *  OpenVolumeMesh is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU Lesser General Public License as           *
 *  published by the Free Software Foundation, either version 3 of           *
 *  the License, or (at your option) any later version with the              *
 *  following exceptions:                                                    *
 *                                                                           *
 *  If other files instantiate templates or use macros                       *
 *  or inline functions from this file, or you compile this file and         *
 *  link it with other files to produce an executable, this file does        *
 *  not by itself cause the resulting executable to be covered by the        *
 *  GNU Lesser General Public License. This exception does not however       *
 *  invalidate any other reasons why the executable file might be            *
 *  covered by the GNU Lesser General Public License.                        *
 *                                                                           *
 *  OpenVolumeMesh is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *  GNU Lesser General Public License for more details.                      *
 *                                                                           *
 *  You should have received a copy of the GNU LesserGeneral Public          *
 *  License along with OpenVolumeMesh.  If not,                              *
 *  see <http://www.gnu.org/licenses/>.                                      *
 *                                                                           *
\*===========================================================================*/

#include "TetrahedralMeshTopologyKernel.hh"

#include <iostream>

namespace OpenVolumeMesh {

FaceHandle TetrahedralMeshTopologyKernel::add_face(const std::vector<HalfEdgeHandle>& _halfedges, bool _topologyCheck) {

    if(_halfedges.size() != 3) {
#ifndef NDEBUG
        std::cerr << "TetrahedralMeshTopologyKernel::add_face(): Face valence is not three! Returning" << std::endl;
        std::cerr << "invalid handle." << std::endl;
#endif
        return TopologyKernel::InvalidFaceHandle;
    }

    return TopologyKernel::add_face(_halfedges, _topologyCheck);
}

//========================================================================================


FaceHandle
TetrahedralMeshTopologyKernel::add_face(const std::vector<VertexHandle>& _vertices) {

    if(_vertices.size() != 3) {
#ifndef NDEBUG
        std::cerr << "TetrahedralMeshTopologyKernel::add_face(): Face valence is not three! Returning" << std::endl;
        std::cerr << "invalid handle." << std::endl;
#endif
        return TopologyKernel::InvalidFaceHandle;
    }

    return TopologyKernel::add_face(_vertices);
}

//========================================================================================


CellHandle
TetrahedralMeshTopologyKernel::add_cell(const std::vector<HalfFaceHandle>& _halffaces, bool _topologyCheck) {

    if(_halffaces.size() != 4) {
// To make this consistent with add_face
#ifndef NDEBUG
        std::cerr << "Cell valence is not four! Aborting." << std::endl;
#endif
        return TopologyKernel::InvalidCellHandle;
    }
    for(std::vector<HalfFaceHandle>::const_iterator it = _halffaces.begin();
            it != _halffaces.end(); ++it) {
        if(TopologyKernel::halfface(*it).halfedges().size() != 3) {
#ifndef NDEBUG
            std::cerr << "Incident face does not have valence three! Aborting." << std::endl;
#endif
            return TopologyKernel::InvalidCellHandle;
        }
    }

    return TopologyKernel::add_cell(_halffaces, _topologyCheck);
}


HalfEdgeHandle TetrahedralMeshTopologyKernel::add_halfedge(const VertexHandle& _fromVertex, const VertexHandle& _toVertex)
{
    HalfEdgeHandle he = halfedge(_fromVertex, _toVertex);
    if (he != InvalidHalfEdgeHandle)
        return he;
    else
        return halfedge_handle(add_edge(_fromVertex, _toVertex), 0);
}

HalfFaceHandle TetrahedralMeshTopologyKernel::add_halfface(const std::vector<HalfEdgeHandle>& _halfedges, bool _topologyCheck)
{
    HalfFaceHandle hf = halfface(_halfedges);
    if (hf != InvalidHalfFaceHandle)
        return hf;
    else
        return halfface_handle(add_face(_halfedges, _topologyCheck), 0);
}

HalfFaceHandle TetrahedralMeshTopologyKernel::add_halfface(VertexHandle _vh0, VertexHandle _vh1, VertexHandle _vh2, bool _topologyCheck)
{
    std::vector<HalfEdgeHandle> halfedges;
    halfedges.push_back(add_halfedge(_vh0, _vh1));
    halfedges.push_back(add_halfedge(_vh1, _vh2));
    halfedges.push_back(add_halfedge(_vh2, _vh0));
    return add_halfface(halfedges, _topologyCheck);
}

/*void TetrahedralMeshTopologyKernel::replaceHalfFace(CellHandle ch, HalfFaceHandle hf_del, HalfFaceHandle hf_ins)
{
    Cell& c = cells_[ch.idx()];
    std::vector<HalfFaceHandle> hfs;
    for (unsigned int i = 0; i < c.halffaces().size(); ++i)
        if (c.halffaces()[i] != hf_del)
            hfs.push_back(c.halffaces()[i]);
        else
            hfs.push_back(hf_ins);
    c.set_halffaces(hfs);
}

void TetrahedralMeshTopologyKernel::replaceHalfEdge(HalfFaceHandle hfh, HalfEdgeHandle he_del, HalfEdgeHandle he_ins)
{
    FaceHandle fh = face_handle(hfh);
    unsigned char oppF = hfh.idx() - halfface_handle(fh, 0);
    if (oppF == 1)
    {
        he_del = opposite_halfedge_handle(he_del);
        he_ins = opposite_halfedge_handle(he_ins);
    }
    Face& f = faces_[fh.idx()];
    std::vector<HalfEdgeHandle> hes;
    for (unsigned int i = 0; i < f.halfedges().size(); ++i)
        if (f.halfedges()[i] != he_del)
            hes.push_back(f.halfedges()[i]);
        else
            hes.push_back(he_ins);
    f.set_halfedges(hes);

}*/

/*
void TetrahedralMeshTopologyKernel::collapse_edge(HalfEdgeHandle _heh)
{

    std::vector<bool> deleteTagFaces(faces_.size(), false);
    std::vector<bool> deleteTagEdges(edges_.size(), false);
    std::vector<bool> deleteTagCells(cells_.size(), false);

    for (HalfEdgeHalfFaceIter hehf_it = hehf_iter(_heh); hehf_it.valid(); ++hehf_it)
    {
        CellHandle ch = incident_cell(*hehf_it);
        if (ch.is_valid())
        {
            HalfFaceHandle hf134 = *hehf_it;
            HalfFaceHandle hf143 = opposite_halfface_handle(hf143);
            HalfEdgeHandle he13  = prev_halfedge_in_halfface(_heh, hf134);
            HalfEdgeHandle he41  = next_halfedge_in_halfface(_heh, hf134);
            HalfFaceHandle hf123 = adjacent_halfface_in_cell(hf134, he13);
            HalfFaceHandle hf142 = adjacent_halfface_in_cell(hf134, he41);
            HalfFaceHandle hf132 = opposite_halfface_handle(hf123);
            CellHandle ch0123 = incident_cell(hf132);
            HalfEdgeHandle he32  = next_halfedge_in_halfface(he13, hf132);
            HalfEdgeHandle he23  = opposite_halfedge_handle(he32);
            HalfEdgeHandle he14  = opposite_halfedge_handle(he41);
            HalfEdgeHandle he31  = opposite_halfedge_handle(he13);
            HalfEdgeHandle he42  = next_halfedge_in_halfface(he14, hf142);
            HalfEdgeHandle he24  = opposite_halfedge_handle(he42);
            HalfFaceHandle hf243 = adjacent_halfface_in_cell(hf134, he41);
            HalfFaceHandle hf234 = opposite_halfface_handle(hf234);
            HalfEdgeHandle he12  = next_halfedge_in_halfface(he42, hf142);
            HalfEdgeHandle he21  = opposite_halfedge_handle(he12);

            if (ch0123.is_valid())
            {
                HalfFaceHandle hf031 = adjacent_halfface_in_cell(hf132, he13);
                HalfFaceHandle hf023 = adjacent_halfface_in_cell(hf132, he32);

                replaceHalfEdge(hf031, he31, he41);
                replaceHalfEdge(hf023, he23, he24);
                replaceHalfFace(ch0123, hf132, hf142);
            }

            //copyHalfFaceAndHalfEdgeProperties(hf132, 142);

            incident_cell_per_hf_[hf142.idx()] = ch0123;


            deleteTagCells[ch.idx()] = true;
            deleteTagFaces[face_handle(hf132).idx()] = true;
            deleteTagFaces[face_handle(*hehf_it).idx()] = true;
            deleteTagEdges[edge_handle(he13).idx()] = true;
            deleteTagEdges[edge_handle(he32).idx()] = true;

            std::set<HalfFaceHandle> excludeFaces;
            excludeFaces.insert(hf134);
            excludeFaces.insert(hf143);
            excludeFaces.insert(hf123);
            excludeFaces.insert(hf132);
            excludeFaces.insert(hf243);
            excludeFaces.insert(hf234);

            std::vector<std::pair<HalfEdgeHandle, HalfEdgeHandle> > joinpartners;
            joinpartners.push_back(std::make_pair(he41, he31));
            joinpartners.push_back(std::make_pair(he14, he13));
            joinpartners.push_back(std::make_pair(he42, he32));
            joinpartners.push_back(std::make_pair(he24, he23));

            for (unsigned int i = 0; i < joinpartners.size(); ++i)
            {
                HalfEdgeHandle target = joinpartners[i].first;
                HalfEdgeHandle source = joinpartners[i].second;
                std::vector<HalfFaceHandle> incidentHfs;
                for (unsigned int j = 0; j < incident_hfs_per_he_[target.idx()].size(); ++j)
                {
                    HalfFaceHandle cur_hf = incident_hfs_per_he_[target.idx()][j];
                    if ((excludeFaces.find(cur_hf) == excludeFaces.end()) && !deleteTagFaces[face_handle(cur_hf).idx()])
                        incidentHfs.push_back(cur_hf);
                }
                for (unsigned int i = 0; i < incident_hfs_per_he_[source.idx()].size(); ++i)
                {
                    HalfFaceHandle cur_hf = incident_hfs_per_he_[source.idx()][i];
                    if ((excludeFaces.find(cur_hf) == excludeFaces.end()) && !deleteTagFaces[face_handle(cur_hf).idx()])
                        incidentHfs.push_back(cur_hf);
                }

                std::swap(incident_hfs_per_he_[target], incidentHfs);
            }

            std::vector<HalfFaceHandle>& vec = incident_hfs_per_he_[he21];
            vec.erase(std::remove(vec.begin(), vec.end(), hf132), vec.end());
            std::vector<HalfFaceHandle>& vec2 = incident_hfs_per_he_[he12];
            vec2.erase(std::remove(vec2.begin(), vec2.end(), hf123), vec2.end());

        }
        else
        {
            deleteTagFaces[face_handle(*hehf_it).idx()] = true;
        }
    }


    VertexHandle from_vh = halfedge(_heh).from_vertex();
    VertexHandle to_vh = halfedge(_heh).to_vertex();
    for (VertexOHalfEdgeIter voh_it = voh_iter(from_vh); voh_it.valid(); ++voh_it )
    {
        Edge he = halfedge(*voh_it);
        if (he.to_vertex() == to_vh)
        {
            std::vector<HalfEdgeHandle>& vec = outgoing_hes_per_vertex_[to_vh];
            vec.erase(std::remove(vec.begin(), vec.end(), opposite_halfedge_handle(*voh_it)), vec.end());
        }
        EdgeHandle eh = edge_handle(*voh_it);
        if (!deleteTagEdges[eh.idx()])
        {
            std::vector<HalfEdgeHandle>& vec = outgoing_hes_per_vertex_[to_vh];
            vec.push_back(opposite_halfedge_handle(*voh_it));

            Edge& e = edges_[eh.idx()];
            if (e.from_vertex() == from_vh)
                e.set_from_vertex(to_vh);
            if (e.to_vertex() == from_vh)
                e.set_to_vertex(to_vh);

        }
    }

    outgoing_hes_per_vertex_[from_vh].clear();

    deleteTagEdges[edge_handle(_heh).idx()] = true;

    delete_multiple_cells(deleteTagCells);
    delete_multiple_faces(deleteTagFaces);
    delete_multiple_edges(deleteTagEdges);
    delete_vertex(from_vh);
}
*/


//void TetrahedralMeshTopologyKernel::swapCellProperties(CellHandle source, CellHandle destination)
//{
//    swapPropertyElements(cell_props_begin(), cell_props_end(), source, destination);
//}

//void TetrahedralMeshTopologyKernel::swapHalfFaceProperties(HalfFaceHandle source, HalfFaceHandle destination)
//{
//    swapPropertyElements(halfface_props_begin(), halfface_props_end(), source, destination);
//}

//void TetrahedralMeshTopologyKernel::swapHalfEdgeProperties(HalfEdgeHandle source, HalfEdgeHandle destination)
//{
//    swapPropertyElements(halfedge_props_begin(), halfedge_props_end(), source, destination);
//}

// cppcheck-suppress unusedFunction ; public interface
VertexHandle TetrahedralMeshTopologyKernel::collapse_edge(HalfEdgeHandle _heh)
{
    bool deferred_deletion_tmp = deferred_deletion_enabled();

    if (!deferred_deletion_tmp)
        enable_deferred_deletion(true);

    VertexHandle from_vh = halfedge(_heh).from_vertex();
    VertexHandle to_vh   = halfedge(_heh).to_vertex();


    // find cells that will collapse, i.e. are incident to the collapsing halfedge
    std::set<CellHandle> collapsingCells;
    for (HalfEdgeHalfFaceIter hehf_it = hehf_iter(_heh); hehf_it.valid(); ++hehf_it)
    {
        HalfFaceHandle hfh = *hehf_it;
        CellHandle ch = incident_cell(hfh);
        if (ch.is_valid())
            collapsingCells.insert(ch);
    }

    std::vector<CellHandle> incidentCells;
    for (VertexCellIter vc_it = vc_iter(from_vh); vc_it.valid(); ++vc_it)
        incidentCells.push_back(*vc_it);

    for (const CellHandle &ch: incidentCells)
    {
        if (collapsingCells.find(ch) != collapsingCells.end())
            continue;

        Cell c = cell(ch);

        std::vector<HalfFaceHandle> newHalffaces;

        for (unsigned int hf_idx = 0; hf_idx < 4; ++hf_idx)
        {
            Face hf = halfface(c.halffaces()[hf_idx]);
            std::vector<HalfEdgeHandle> newHalfedges;

            for (unsigned int j = 0; j < 3; ++j)
            {
                Edge e = halfedge(hf.halfedges()[j]);
                VertexHandle newStart = (e.from_vertex() == from_vh) ? to_vh: e.from_vertex();
                VertexHandle newEnd   = (e.to_vertex()   == from_vh) ? to_vh : e.to_vertex();

                HalfEdgeHandle heh = add_halfedge(newStart, newEnd);
                newHalfedges.push_back(heh);
                swap_halfedge_properties(hf.halfedges()[j], heh);
            }

            HalfFaceHandle hfh = add_halfface(newHalfedges);
            newHalffaces.push_back(hfh);
            swap_halfface_properties(c.halffaces()[hf_idx], hfh);
        }

        delete_cell(ch);

        CellHandle newCell = add_cell(newHalffaces);

        swap_cell_properties(ch, newCell);

    }


    VertexHandle survivingVertex = to_vh;

    if (!deferred_deletion_tmp)
    {
        if (fast_deletion_enabled())
        {
            // from_vh is swapped with last vertex and then deleted
            if (to_vh.idx() == (int)n_vertices() - 1)
                survivingVertex = from_vh;
        }
        else
        {
            // from_vh is deleted and every vertex id larger than from_vh is reduced by one
            if (from_vh.idx() < to_vh.idx())
                survivingVertex = VertexHandle(to_vh.idx() - 1);
        }
    }

    delete_vertex(from_vh);

    enable_deferred_deletion(deferred_deletion_tmp);

    return survivingVertex;

}

// cppcheck-suppress unusedFunction ; public interface
void TetrahedralMeshTopologyKernel::split_edge(HalfEdgeHandle _heh, VertexHandle _vh)
{
    bool deferred_deletion_tmp = deferred_deletion_enabled();

    if (!deferred_deletion_tmp)
        enable_deferred_deletion(true);

    std::vector<HalfFaceHandle> incident_halffaces_with_cells;
    for (HalfEdgeHalfFaceIter hehf_it = hehf_iter(_heh); hehf_it.valid(); ++hehf_it)
    {
        CellHandle ch = incident_cell(*hehf_it);
        if (ch.is_valid())
            incident_halffaces_with_cells.push_back(*hehf_it);
    }

    for (auto hfh : incident_halffaces_with_cells)
    {
        CellHandle ch = incident_cell(hfh);

        std::vector<VertexHandle> vertices = get_cell_vertices(hfh, _heh);

        delete_cell(ch);

        add_cell(vertices[0], _vh, vertices[2], vertices[3]);
        add_cell(_vh, vertices[1], vertices[2], vertices[3]);
    }

    delete_edge(edge_handle(_heh));

    enable_deferred_deletion(deferred_deletion_tmp);

}

// cppcheck-suppress unusedFunction ; public interface
void TetrahedralMeshTopologyKernel::split_face(FaceHandle _fh, VertexHandle _vh)
{
    bool deferred_deletion_tmp = deferred_deletion_enabled();

    if (!deferred_deletion_tmp)
        enable_deferred_deletion(true);

    for (char i = 0; i < 2; ++i)
    {
        HalfFaceHandle hfh = halfface_handle(_fh, i);
        CellHandle ch = incident_cell(hfh);
        if (ch.is_valid())
        {
            std::vector<VertexHandle> vertices = get_cell_vertices(hfh);

            delete_cell(ch);

            add_cell(vertices[0], vertices[1], _vh, vertices[3]);
            add_cell(vertices[0], _vh, vertices[2], vertices[3]);
            add_cell(_vh, vertices[1], vertices[2], vertices[3]);
        }
    }

    delete_face(_fh);

    enable_deferred_deletion(deferred_deletion_tmp);

}


std::vector<VertexHandle> TetrahedralMeshTopologyKernel::get_cell_vertices(CellHandle ch) const
{
    return get_cell_vertices(cell(ch).halffaces().front());
}

std::vector<VertexHandle> TetrahedralMeshTopologyKernel::get_cell_vertices(CellHandle ch, VertexHandle vh) const
{
    HalfFaceHandle hfh = cell(ch).halffaces()[0];
    Face f = halfface(hfh);
    HalfEdgeHandle heh;
    for (unsigned int i = 0; i < 3; ++i)
    {
        Edge e = halfedge(f.halfedges()[i]);
        if (e.from_vertex() == vh)
        {
            heh = f.halfedges()[i];
            break;
        }
    }
    if (!heh.is_valid())
    {
        hfh = adjacent_halfface_in_cell(hfh, f.halfedges()[0]);
        heh = prev_halfedge_in_halfface(opposite_halfedge_handle(f.halfedges()[0]), hfh);
    }

    return get_cell_vertices(hfh,heh);

}

std::vector<VertexHandle> TetrahedralMeshTopologyKernel::get_cell_vertices(HalfFaceHandle hfh) const
{
    return get_cell_vertices(hfh, halfface(hfh).halfedges().front());
}

std::vector<VertexHandle> TetrahedralMeshTopologyKernel::get_cell_vertices(HalfFaceHandle hfh, HalfEdgeHandle heh) const
{
    std::vector<VertexHandle> vertices;

    // add vertices of halfface
    for (unsigned int i = 0; i < 3; ++i)
    {
        Edge e = halfedge(heh);
        vertices.push_back(e.from_vertex());
        heh = next_halfedge_in_halfface(heh, hfh);
    }

    Cell c = cell(incident_cell(hfh));
    HalfFaceHandle otherHfh = c.halffaces()[0];
    if (otherHfh == hfh)
        otherHfh = c.halffaces()[1];

    Face otherF = halfface(otherHfh);

    for (unsigned int i = 0; i < otherF.halfedges().size(); ++i)
    {
        HalfEdgeHandle he = otherF.halfedges()[i];
        Edge e = halfedge(he);
        if (std::find(vertices.begin(), vertices.end(), e.to_vertex()) == vertices.end())
        {
            vertices.push_back(e.to_vertex());
            return vertices;
        }
    }

    return vertices;
}

std::vector<VertexHandle> TetrahedralMeshTopologyKernel::get_halfface_vertices(HalfFaceHandle hfh) const
{
    return get_halfface_vertices(hfh, halfface(hfh).halfedges().front());
}

std::vector<VertexHandle> TetrahedralMeshTopologyKernel::get_halfface_vertices(HalfFaceHandle hfh, VertexHandle vh) const
{
    Face hf = halfface(hfh);
    for (unsigned int i = 0; i < 3; ++i)
        if (halfedge(hf.halfedges()[i]).from_vertex() == vh)
            return get_halfface_vertices(hfh, hf.halfedges()[i]);

    return std::vector<VertexHandle>();
}

// cppcheck-suppress unusedFunction ; public interface
std::vector<VertexHandle> TetrahedralMeshTopologyKernel::get_halfface_vertices(HalfFaceHandle hfh, HalfEdgeHandle heh) const
{
    std::vector<VertexHandle> vertices;

    // add vertices of halfface
    for (unsigned int i = 0; i < 3; ++i)
    {
        Edge e = halfedge(heh);
        vertices.push_back(e.from_vertex());
        heh = next_halfedge_in_halfface(heh, hfh);
    }

    return vertices;
}

VertexHandle TetrahedralMeshTopologyKernel::halfface_opposite_vertex(HalfFaceHandle hfh) const
{
    if (is_boundary(hfh)) {
        return InvalidVertexHandle;
    }

    const std::vector<VertexHandle> base = get_halfface_vertices(hfh);
    for (CellVertexIter it = cv_iter(incident_cell(hfh)); it.valid(); ++it) {
        const VertexHandle vh = *it;
        if (vh != base[0] && vh != base[1] && vh != base[2]) {
            return vh;
        }
    }

    return InvalidVertexHandle;
}


//========================================================================================

CellHandle
TetrahedralMeshTopologyKernel::add_cell(const std::vector<VertexHandle>& _vertices, bool _topologyCheck) {

    // debug mode checks
    assert(TopologyKernel::has_full_bottom_up_incidences());
    assert(_vertices.size() == 4);

    // release mode checks
    if(!TopologyKernel::has_full_bottom_up_incidences()) {
        return CellHandle(-1);
    }

    if(_vertices.size() != 4) {
        return CellHandle(-1);
    }

    HalfFaceHandle hf0, hf1, hf2, hf3;

    std::vector<VertexHandle> vs;

    vs.push_back(_vertices[0]);
    vs.push_back(_vertices[1]);
    vs.push_back(_vertices[2]);
    hf0 = TopologyKernel::halfface(vs);
    if(!hf0.is_valid()) {
        FaceHandle fh = TopologyKernel::add_face(vs);
        hf0 = halfface_handle(fh, 0);
    }
    vs.clear();

    vs.push_back(_vertices[0]);
    vs.push_back(_vertices[2]);
    vs.push_back(_vertices[3]);
    hf1 = TopologyKernel::halfface(vs);
    if(!hf1.is_valid()) {
        FaceHandle fh = TopologyKernel::add_face(vs);
        hf1 = halfface_handle(fh, 0);
    }
    vs.clear();

    vs.push_back(_vertices[0]);
    vs.push_back(_vertices[3]);
    vs.push_back(_vertices[1]);
    hf2 = TopologyKernel::halfface(vs);
    if(!hf2.is_valid()) {
        FaceHandle fh = TopologyKernel::add_face(vs);
        hf2 = halfface_handle(fh, 0);
    }
    vs.clear();

    vs.push_back(_vertices[1]);
    vs.push_back(_vertices[3]);
    vs.push_back(_vertices[2]);
    hf3 = TopologyKernel::halfface(vs);
    if(!hf3.is_valid()) {
        FaceHandle fh = TopologyKernel::add_face(vs);
        hf3 = halfface_handle(fh, 0);
    }
    vs.clear();

    assert(hf0.is_valid());
    assert(hf1.is_valid());
    assert(hf2.is_valid());
    assert(hf3.is_valid());


    std::vector<HalfFaceHandle> hfs;
    hfs.push_back(hf0);
    hfs.push_back(hf1);
    hfs.push_back(hf2);
    hfs.push_back(hf3);

    if (_topologyCheck) {
        /*
        * Test if all halffaces are connected and form a two-manifold
        * => Cell is closed
        *
        * This test is simple: The number of involved half-edges has to be
        * exactly twice the number of involved edges.
        */

        std::set<HalfEdgeHandle> incidentHalfedges;
        std::set<EdgeHandle>     incidentEdges;

        for(std::vector<HalfFaceHandle>::const_iterator it = hfs.begin(),
                end = hfs.end(); it != end; ++it) {

            OpenVolumeMeshFace hface = halfface(*it);
            for(std::vector<HalfEdgeHandle>::const_iterator he_it = hface.halfedges().begin(),
                    he_end = hface.halfedges().end(); he_it != he_end; ++he_it) {
                incidentHalfedges.insert(*he_it);
                incidentEdges.insert(edge_handle(*he_it));
            }
        }

        if(incidentHalfedges.size() != (incidentEdges.size() * 2u)) {
#ifndef NDEBUG
            std::cerr << "The specified halffaces are not connected!" << std::endl;
#endif
            return InvalidCellHandle;
        }
        // The halffaces are now guaranteed to form a two-manifold

        if(has_face_bottom_up_incidences()) {

            for(std::vector<HalfFaceHandle>::const_iterator it = hfs.begin(),
                    end = hfs.end(); it != end; ++it) {
                if(incident_cell(*it) != InvalidCellHandle) {
#ifndef NDEBUG
                    std::cerr << "Warning: One of the specified half-faces is already incident to another cell!" << std::endl;
#endif
                    return InvalidCellHandle;
                }
            }

        }

    }

    return TopologyKernel::add_cell(hfs, false);
}

CellHandle TetrahedralMeshTopologyKernel::add_cell(VertexHandle _vh0, VertexHandle _vh1, VertexHandle _vh2, VertexHandle _vh3, bool _topologyCheck)
{
    std::vector<HalfFaceHandle> halffaces;
    halffaces.push_back(add_halfface(_vh0, _vh1, _vh2));
    halffaces.push_back(add_halfface(_vh0, _vh2, _vh3));
    halffaces.push_back(add_halfface(_vh0, _vh3, _vh1));
    halffaces.push_back(add_halfface(_vh1, _vh3, _vh2));
    return add_cell(halffaces, _topologyCheck);
}

//========================================================================================

} // Namespace OpenVolumeMesh
