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

#define NORMALATTRIBT_CC

#include <set>

#include "NormalAttrib.hh"

#include "../Core/GeometryKernel.hh"

namespace OpenVolumeMesh {

template <class GeomKernelT>
NormalAttrib<GeomKernelT>::NormalAttrib(GeomKernelT& _kernel) :
kernel_(_kernel),
v_normals_(_kernel.template request_vertex_property<typename GeomKernelT::PointT>("vertex_normals", typename GeomKernelT::PointT(0.0))),
f_normals_(_kernel.template request_face_property<typename GeomKernelT::PointT>("face_normals", typename GeomKernelT::PointT(0.0)))
{

}

template <class GeomKernelT>
NormalAttrib<GeomKernelT>::~NormalAttrib() {

}

template <class GeomKernelT>
void NormalAttrib<GeomKernelT>::update_vertex_normals() {

    if(!kernel_.has_face_bottom_up_incidences()) {
        std::cerr << "Error: update_vertex_normals() needs bottom-up incidences!" << std::endl;
        return;
    }

    // Compute face normals
    update_face_normals();

    for(const auto &_vh: kernel_.vertices())
    {
        compute_vertex_normal(_vh);
    }
}

template <class GeomKernelT>
void NormalAttrib<GeomKernelT>::update_face_normals() {

    if(!kernel_.has_face_bottom_up_incidences()) {
        std::cerr << "Error: update_normals() needs bottom-up incidences!" << std::endl;
        return;
    }

    for (const auto &_fh: kernel_.faces()) {
        f_normals_[_fh] = kernel_.normal(kernel_.halfface_handle(_fh, 0));
    }
}

template <class GeomKernelT>
void NormalAttrib<GeomKernelT>::compute_vertex_normal(const VertexHandle& _vh) {

    std::set<HalfFaceHandle> halffaces;
    for(VertexOHalfEdgeIter voh_it = kernel_.voh_iter(_vh);
            voh_it.valid(); ++voh_it) {

        for(HalfEdgeHalfFaceIter hehf_it = kernel_.hehf_iter(*voh_it);
                hehf_it.valid(); ++hehf_it) {
            if(kernel_.is_boundary(*hehf_it)) {
                halffaces.insert(*hehf_it);
            }
        }
    }
    typename GeomKernelT::PointT normal = typename GeomKernelT::PointT(0.0);
    for(std::set<HalfFaceHandle>::const_iterator hf_it = halffaces.begin();
            hf_it != halffaces.end(); ++hf_it) {
        normal += (*this)[*hf_it];
    }

    normal.normalize();

}

} // Namespace OpenVolumeMesh
