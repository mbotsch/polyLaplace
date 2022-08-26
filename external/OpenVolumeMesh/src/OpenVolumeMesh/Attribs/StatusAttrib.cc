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

/*===========================================================================*\
 *                                                                           *
 *   $Revision$                                                         *
 *   $Date$                    *
 *   $LastChangedBy$                                                *
 *                                                                           *
\*===========================================================================*/

#ifndef NDEBUG
#include <iostream>
#endif

#include "../Core/TopologyKernel.hh"
#include "StatusAttrib.hh"

namespace OpenVolumeMesh {

StatusAttrib::StatusAttrib(TopologyKernel& _kernel) :
kernel_(_kernel),
v_status_(_kernel.request_vertex_property<OpenVolumeMeshStatus>("vertex_status")),
e_status_(_kernel.request_edge_property<OpenVolumeMeshStatus>("edge_status")),
he_status_(_kernel.request_halfedge_property<OpenVolumeMeshStatus>("halfedge_status")),
f_status_(_kernel.request_face_property<OpenVolumeMeshStatus>("face_status")),
hf_status_(_kernel.request_halfface_property<OpenVolumeMeshStatus>("halfface_status")),
c_status_(_kernel.request_cell_property<OpenVolumeMeshStatus>("cell_status")),
m_status_(_kernel.request_mesh_property<OpenVolumeMeshStatus>("mesh_status")) {

}

//========================================================================================

StatusAttrib::~StatusAttrib() {

}

//========================================================================================

void StatusAttrib::mark_higher_dim_entities() {

    // Edges
    if(kernel_.has_vertex_bottom_up_incidences()) {

        for(VertexIter v_it = kernel_.vertices_begin(); v_it != kernel_.vertices_end(); ++v_it) {
            if(v_status_[*v_it].deleted()) {
                for(VertexOHalfEdgeIter voh_it = kernel_.voh_iter(*v_it);
                        voh_it.valid(); ++voh_it) {
                    e_status_[kernel_.edge_handle(*voh_it)].set_deleted(true);
                }
            }
        }
    } else {

        for(EdgeIter e_it = kernel_.edges_begin(); e_it != kernel_.edges_end(); ++e_it) {
            if(v_status_[kernel_.edge(*e_it).from_vertex()].deleted() ||
                    v_status_[kernel_.edge(*e_it).to_vertex()].deleted()) {
                e_status_[*e_it].set_deleted(true);
            }
        }
    }

    // Faces
    if(kernel_.has_edge_bottom_up_incidences()) {

        for(EdgeIter e_it = kernel_.edges_begin(); e_it != kernel_.edges_end(); ++e_it) {
            if(e_status_[*e_it].deleted()) {
                for(HalfEdgeHalfFaceIter hehf_it = kernel_.hehf_iter(kernel_.halfedge_handle(*e_it, 0));
                        hehf_it.valid(); ++hehf_it) {
                    f_status_[kernel_.face_handle(*hehf_it)].set_deleted(true);
                }
            }
        }
    } else {

        for(FaceIter f_it = kernel_.faces_begin(); f_it != kernel_.faces_end(); ++f_it) {

            const std::vector<HalfEdgeHandle>& hes = kernel_.face(*f_it).halfedges();
            for(std::vector<HalfEdgeHandle>::const_iterator he_it = hes.begin(),
                    he_end = hes.end(); he_it != he_end; ++he_it) {
                if(e_status_[kernel_.edge_handle(*he_it)].deleted()) {
                    f_status_[*f_it].set_deleted(true);
                    break;
                }
            }
        }
    }

    // Cells
    if(kernel_.has_face_bottom_up_incidences()) {

        for(FaceIter f_it = kernel_.faces_begin(); f_it != kernel_.faces_end(); ++f_it) {
            if(f_status_[*f_it].deleted()) {
                CellHandle c0 = kernel_.incident_cell(kernel_.halfface_handle(*f_it, 0));
                CellHandle c1 = kernel_.incident_cell(kernel_.halfface_handle(*f_it, 1));
                if(c0.is_valid()) {
                    c_status_[c0].set_deleted(true);
                }
                if(c1.is_valid()) {
                    c_status_[c1].set_deleted(true);
                }
            }
        }
    } else {

        for(CellIter c_it = kernel_.cells_begin(); c_it != kernel_.cells_end(); ++c_it) {

            const std::vector<HalfFaceHandle>& hfs = kernel_.cell(*c_it).halffaces();
            for(std::vector<HalfFaceHandle>::const_iterator hf_it = hfs.begin(),
                    hf_end = hfs.end(); hf_it != hf_end; ++hf_it) {
                if(f_status_[kernel_.face_handle(*hf_it)].deleted()) {
                    c_status_[*c_it].set_deleted(true);
                    break;
                }
            }
        }
    }
}

//========================================================================================

void StatusAttrib::garbage_collection(bool _preserveManifoldness) {
    std::vector<VertexHandle*> vh_empty;
    std::vector<HalfEdgeHandle*> hh_empty;
    std::vector<HalfFaceHandle*> hfh_empty;
    std::vector<CellHandle*> ch_empty;

    garbage_collection(vh_empty, hh_empty, hfh_empty, ch_empty, _preserveManifoldness);
}


} // Namespace OpenVolumeMesh
