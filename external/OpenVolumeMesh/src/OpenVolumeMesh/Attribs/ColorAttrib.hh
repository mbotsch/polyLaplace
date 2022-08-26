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

#ifndef COLORATTRIB_HH_
#define COLORATTRIB_HH_

#include <cassert>

#include "../Core/OpenVolumeMeshHandle.hh"
#include "OpenVolumeMeshStatus.hh"
#include "../Core/PropertyDefines.hh"
#include "../Core/TopologyKernel.hh"

namespace OpenVolumeMesh {

//== CLASS DEF ================================================================

template <class ColT>
class ColorAttrib {
public:

    ColorAttrib(TopologyKernel& _kernel, const ColT _def = ColT());

    virtual ~ColorAttrib();

    //==================
    // Vertices
    //==================
    const ColT& operator[](const VertexHandle& _h) const {
        assert((unsigned int)_h.idx() < kernel_.n_vertices());
        return vcolor_prop_[_h];
    }

    ColT& operator[](const VertexHandle& _h) {
        assert((unsigned int)_h.idx() < kernel_.n_vertices());
        vertex_colors_available_ = true;
        return vcolor_prop_[_h];
    }

    //==================
    // Edges
    //==================
    const ColT& operator[](const EdgeHandle& _h) const {
        assert((unsigned int)_h.idx() < kernel_.n_edges());
        return ecolor_prop_[_h];
    }

    ColT& operator[](const EdgeHandle& _h) {
        assert((unsigned int)_h.idx() < kernel_.n_edges());
        edge_colors_available_ = true;
        return ecolor_prop_[_h];
    }

    //==================
    // Half-Edges
    //==================
    const ColT& operator[](const HalfEdgeHandle& _h) const {
        assert((unsigned int)_h.idx() < kernel_.n_halfedges());
        return hecolor_prop_[_h];
    }

    ColT& operator[](const HalfEdgeHandle& _h) {
        assert((unsigned int)_h.idx() < kernel_.n_halfedges());
        halfedge_colors_available_ = true;
        return hecolor_prop_[_h];
    }

    //==================
    // Faces
    //==================
    const ColT& operator[](const FaceHandle& _h) const {
        assert((unsigned int)_h.idx() < kernel_.n_faces());
        return fcolor_prop_[_h];
    }

    ColT& operator[](const FaceHandle& _h) {
        assert((unsigned int)_h.idx() < kernel_.n_faces());
        face_colors_available_ = true;
        return fcolor_prop_[_h];
    }

    //==================
    // Half-Faces
    //==================
    const ColT& operator[](const HalfFaceHandle& _h) const {
        assert((unsigned int)_h.idx() < kernel_.n_halffaces());
        return hfcolor_prop_[_h];
    }

    ColT& operator[](const HalfFaceHandle& _h) {
        assert((unsigned int)_h.idx() < kernel_.n_halffaces());
        halfface_colors_available_ = true;
        return hfcolor_prop_[_h];
    }

    //==================
    // Cells
    //==================
    const ColT& operator[](const CellHandle& _h) const {
        assert((unsigned int)_h.idx() < kernel_.n_cells());
        return ccolor_prop_[_h];
    }

    ColT& operator[](const CellHandle& _h) {
        assert((unsigned int)_h.idx() < kernel_.n_cells());
        cell_colors_available_ = true;
        return ccolor_prop_[_h];
    }


    bool vertex_colors_available()   { return vertex_colors_available_;   }
    bool halfedge_colors_available() { return halfedge_colors_available_; }
    bool edge_colors_available()     { return edge_colors_available_;     }
    bool halfface_colors_available() { return halfface_colors_available_; }
    bool face_colors_available()     { return face_colors_available_;     }
    bool cell_colors_available()     { return cell_colors_available_;     }

    void clear_vertex_colors();
    void clear_halfedge_colors();
    void clear_edge_colors();
    void clear_halfface_colors();
    void clear_face_colors();
    void clear_cell_colors();


private:

    VertexPropertyT<ColT> vcolor_prop_;
    EdgePropertyT<ColT> ecolor_prop_;
    HalfEdgePropertyT<ColT> hecolor_prop_;
    FacePropertyT<ColT> fcolor_prop_;
    HalfFacePropertyT<ColT> hfcolor_prop_;
    CellPropertyT<ColT> ccolor_prop_;

    TopologyKernel& kernel_;

    bool vertex_colors_available_;
    bool halfedge_colors_available_;
    bool edge_colors_available_;
    bool halfface_colors_available_;
    bool face_colors_available_;
    bool cell_colors_available_;

    ColT default_color_;

};

} // Namespace OpenVolumeMesh

#if defined(INCLUDE_TEMPLATES) && !defined(COLORATTRIBT_CC)
#include "ColorAttribT_impl.hh"
#endif

#endif /* COLORATTRIB_HH_ */
