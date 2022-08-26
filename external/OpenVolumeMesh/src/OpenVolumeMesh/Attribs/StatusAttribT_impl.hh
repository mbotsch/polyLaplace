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
 *   $Revision: 36 $                                                         *
 *   $Date: 2012-01-10 18:00:06 +0100 (Di, 10 Jan 2012) $                    *
 *   $LastChangedBy: kremer $                                                *
 *                                                                           *
\*===========================================================================*/

#define STATUSATTRIBT_CC

#include "StatusAttrib.hh"

#include "../Core/TopologyKernel.hh"
#include "../Core/PropertyDefines.hh"

#include <map>

namespace OpenVolumeMesh {
//========================================================================================

template<typename std_API_Container_VHandlePointer,
         typename std_API_Container_HHandlePointer,
         typename std_API_Container_HFHandlePointer,
         typename std_API_Container_CHandlePointer>
void StatusAttrib::garbage_collection(std_API_Container_VHandlePointer &vh_to_update,
                                      std_API_Container_HHandlePointer &hh_to_update,
                                      std_API_Container_HFHandlePointer &hfh_to_update,
                                      std_API_Container_CHandlePointer &ch_to_update,
                                      bool _preserveManifoldness) {

    /*
     * This is not a real garbage collection in its conventional
     * sense. What happens in this routine are the following steps:
     *
     * 1. If an entity of dimension n is marked to be deleted,
     *    also mark all incident entities of dimension n + 1
     *    for deletion. Do this in a bottom-up fashion.
     * 2. Then delete all entities in top-down manner, so that
     *    no invalid incident higher-dimensional entity is generated.
     * 3. If desired, search for all isolated entities and mark
     *    them deleted in a top-down manner.
     * 4. Delete all entities marked deleted in step 4 in order
     *    to prevent manifoldness.
     */

    // setup tracking so we can update the given handles
    bool track_vh = !vh_to_update.empty();
    bool track_hh = !hh_to_update.empty();
    bool track_hfh = !hfh_to_update.empty();
    bool track_ch = !ch_to_update.empty();
    int offset_vh = 0;
    int offset_hh = 0;
    int offset_hfh = 0;
    int offset_ch = 0;

    std::map<int,int> vh_map;
    std::map<int,int> hh_map;
    std::map<int,int> hfh_map;
    std::map<int,int> ch_map;

    // initialise the maps
    if (track_vh) {
        typename std_API_Container_VHandlePointer::iterator it = vh_to_update.begin();
        typename std_API_Container_VHandlePointer::iterator end = vh_to_update.end();

        for (; it != end; ++it) {
            vh_map[(*it)->idx()] = (*it)->idx();
        }
    }
    if (track_hh) {
        typename std_API_Container_HHandlePointer::iterator it = hh_to_update.begin();
        typename std_API_Container_HHandlePointer::iterator end = hh_to_update.end();

        for (; it != end; ++it) {
            hh_map[(*it)->idx()] = (*it)->idx();
        }
    }
    if (track_hfh) {
        typename std_API_Container_HFHandlePointer::iterator it = hfh_to_update.begin();
        typename std_API_Container_HFHandlePointer::iterator end = hfh_to_update.end();

        for (; it != end; ++it) {
            hfh_map[(*it)->idx()] = (*it)->idx();
        }
    }
    if (track_ch) {
        typename std_API_Container_CHandlePointer::iterator it = ch_to_update.begin();
        typename std_API_Container_CHandlePointer::iterator end = ch_to_update.end();

        for (; it != end; ++it) {
            ch_map[(*it)->idx()] = (*it)->idx();
        }
    }

    // Mark all higher-dimensional entities incident to
    // entities marked as deleted from bottom to top
    mark_higher_dim_entities();

    std::vector<int> vertexIndexMap(kernel_.n_vertices(), -1);

    // Turn off bottom-up incidences
    bool v_bu = kernel_.has_vertex_bottom_up_incidences();
    bool e_bu = kernel_.has_edge_bottom_up_incidences();
    bool f_bu = kernel_.has_face_bottom_up_incidences();

    kernel_.enable_bottom_up_incidences(false);

    std::vector<bool> tags(kernel_.n_cells(), false);
    std::vector<bool>::iterator tag_it = tags.begin();

    for(CellIter c_it = kernel_.cells_begin(); c_it != kernel_.cells_end(); ++c_it, ++tag_it) {
        *tag_it = c_status_[*c_it].deleted();

        if (track_ch) {
            if (c_status_[*c_it].deleted()) {
                ++offset_ch;
                if (ch_map.find(c_it->idx()) != ch_map.end())
                    ch_map[c_it->idx()] = -1;
            } else {
                if (ch_map.find(c_it->idx()) != ch_map.end())
                    ch_map[c_it->idx()] = c_it->idx() - offset_ch;
            }
        }
    }
    kernel_.delete_multiple_cells(tags);

    tags.resize(kernel_.n_faces(), false);
    tag_it = tags.begin();

    for(FaceIter f_it = kernel_.faces_begin(); f_it != kernel_.faces_end(); ++f_it, ++tag_it) {
        *tag_it = f_status_[*f_it].deleted();

        if (track_hfh) {
            int halfface_idx = f_it->idx() * 2;
            if (f_status_[*f_it].deleted()) {
                offset_hfh += 2;
                if (hfh_map.find(halfface_idx) != hfh_map.end()) {
                    hfh_map[halfface_idx] = -1;
                }
                if (hfh_map.find(halfface_idx + 1) != hfh_map.end()) {
                    hfh_map[halfface_idx + 1] = -1;
                }
            } else {
                if (hfh_map.find(halfface_idx) != hfh_map.end()) {
                    hfh_map[halfface_idx] = halfface_idx - offset_hfh;
                }
                if (hfh_map.find(halfface_idx + 1) != hfh_map.end()) {
                    hfh_map[halfface_idx + 1] = halfface_idx + 1 - offset_hfh;
                }
            }
        }
    }
    kernel_.delete_multiple_faces(tags);

    tags.resize(kernel_.n_edges(), false);
    tag_it = tags.begin();

    for(EdgeIter e_it = kernel_.edges_begin(); e_it != kernel_.edges_end(); ++e_it, ++tag_it) {
        *tag_it = e_status_[*e_it].deleted();

        if (track_hh) {
            int halfedge_idx = e_it->idx() * 2;
            if (e_status_[*e_it].deleted()) {
                offset_hh += 2;
                if (hh_map.find(halfedge_idx) != hh_map.end()) {
                    hh_map[halfedge_idx] = -1;
                }
                if (hh_map.find(halfedge_idx + 1) != hh_map.end()) {
                    hh_map[halfedge_idx + 1] = -1;
                }
            } else {
                if (hh_map.find(halfedge_idx) != hh_map.end()) {
                    hh_map[halfedge_idx] = halfedge_idx - offset_hh;
                }
                if (hh_map.find(halfedge_idx + 1) != hh_map.end()) {
                    hh_map[halfedge_idx + 1] = halfedge_idx + 1 - offset_hh;
                }
            }
        }
    }
    kernel_.delete_multiple_edges(tags);

    tags.resize(kernel_.n_vertices(), false);
    tag_it = tags.begin();

    for(VertexIter v_it = kernel_.vertices_begin(); v_it != kernel_.vertices_end(); ++v_it, ++tag_it) {
        *tag_it = v_status_[*v_it].deleted();

        if (track_vh) {
            if (v_status_[*v_it].deleted()) {
                if (vh_map.find(v_it->idx()) != vh_map.end()) {
                    ++offset_vh;
                    vh_map[v_it->idx()] = -1;
                }
            } else {
                if (vh_map.find(v_it->idx()) != vh_map.end()) {
                    vh_map[v_it->idx()] = v_it->idx() - offset_vh;
                }
            }
        }
    }
    kernel_.delete_multiple_vertices(tags);

    // update given handles
    if (track_vh) {
        typename std_API_Container_VHandlePointer::iterator it = vh_to_update.begin();
        typename std_API_Container_VHandlePointer::iterator end = vh_to_update.end();

        for (; it != end; ++it) {
            *(*it) = VertexHandle( vh_map[(*it)->idx()] );
        }
    }
    if (track_hh) {
        typename std_API_Container_HHandlePointer::iterator it = hh_to_update.begin();
        typename std_API_Container_HHandlePointer::iterator end = hh_to_update.end();

        for (; it != end; ++it) {
            *(*it) = HalfEdgeHandle( hh_map[(*it)->idx()] );
        }
    }
    if (track_hfh) {
        typename std_API_Container_HFHandlePointer::iterator it = hfh_to_update.begin();
        typename std_API_Container_HFHandlePointer::iterator end = hfh_to_update.end();

        for (; it != end; ++it) {
            *(*it) = HalfFaceHandle( hfh_map[(*it)->idx()] );
        }
    }
    if (track_ch) {
        typename std_API_Container_CHandlePointer::iterator it = ch_to_update.begin();
        typename std_API_Container_CHandlePointer::iterator end = ch_to_update.end();

        for (; it != end; ++it) {
            *(*it) = CellHandle( ch_map[(*it)->idx()] );
        }
    }

    // Todo: Resize props

    if(v_bu) kernel_.enable_vertex_bottom_up_incidences(true);
    if(e_bu) kernel_.enable_edge_bottom_up_incidences(true);
    if(f_bu) kernel_.enable_face_bottom_up_incidences(true);

    // Step 6
    if(_preserveManifoldness) {
        if(kernel_.has_full_bottom_up_incidences()) {

            // Go over all faces and find those
            // that are not incident to any cell
            for(FaceIter f_it = kernel_.faces_begin(); f_it != kernel_.faces_end(); ++f_it) {

                // Get half-faces
                HalfFaceHandle hf0 = kernel_.halfface_handle(*f_it, 0);
                HalfFaceHandle hf1 = kernel_.halfface_handle(*f_it, 1);

                // If neither of the half-faces is incident to a cell, delete face
                if(kernel_.incident_cell(hf0) == TopologyKernel::InvalidCellHandle &&
                        kernel_.incident_cell(hf1) == TopologyKernel::InvalidCellHandle) {

                    f_status_[*f_it].set_deleted(true);
                }
            }

            // Go over all edges and find those
            // whose half-edges are not incident to any half-face
            for(EdgeIter e_it = kernel_.edges_begin(); e_it != kernel_.edges_end(); ++e_it) {

                // Get half-edges
                HalfEdgeHandle he = kernel_.halfedge_handle(*e_it, 0);

                // If the half-edge isn't incident to a half-face, delete edge
                HalfEdgeHalfFaceIter hehf_it = kernel_.hehf_iter(he);

                if(!hehf_it.valid()) {

                    e_status_[*e_it].set_deleted(true);

                } else {
                    bool validFace = false;
                    for(; hehf_it.valid(); ++hehf_it) {
                        if(!f_status_[kernel_.face_handle(*hehf_it)].deleted()) {
                            validFace = true;
                            break;
                        }
                    }
                    if(!validFace) {
                        e_status_[*e_it].set_deleted(true);
                    }
                }
            }

            // Go over all vertices and find those
            // that are not incident to any edge
            for(VertexIter v_it = kernel_.vertices_begin(); v_it != kernel_.vertices_end(); ++v_it) {

                // If neither of the half-edges is incident to a half-face, delete edge
                VertexOHalfEdgeIter voh_it = kernel_.voh_iter(*v_it);

                if(!voh_it.valid()) {

                    v_status_[*v_it].set_deleted(true);
                } else {

                    bool validEdge = false;
                    for(; voh_it.valid(); ++voh_it) {
                        if(!e_status_[kernel_.edge_handle(*voh_it)].deleted()) {
                            validEdge = true;
                            break;
                        }
                    }
                    if(!validEdge) {
                        v_status_[*v_it].set_deleted(true);
                    }
                }
            }

            // Recursive call
            garbage_collection(vh_to_update, hh_to_update, hfh_to_update, ch_to_update, false);

        } else {
#ifndef NDEBUG
            std::cerr << "Preservation of three-manifoldness in garbage_collection() "
                    << "requires bottom-up incidences!" << std::endl;
#endif
            return;
        }
    }
}
} // Namespace OpenVolumeMesh
