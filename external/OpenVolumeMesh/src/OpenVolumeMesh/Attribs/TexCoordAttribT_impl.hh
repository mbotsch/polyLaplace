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

#define TEXCOORDATTRIBT_CC

#include "TexCoordAttrib.hh"

namespace OpenVolumeMesh {

template <class TexCoordT>
TexCoordAttrib<TexCoordT>::TexCoordAttrib(TopologyKernel& _kernel, const TexCoordT _def) :
        vtexcoord_prop_(_kernel.request_vertex_property<TexCoordT>("vertex_texcoord", _def)),
        kernel_(_kernel),
        vertex_texcoords_available_(false),
        default_texcoord_(_def)
{

}

template <class TexCoordT>
TexCoordAttrib<TexCoordT>::~TexCoordAttrib() {

}

template <class TexCoordT>
void TexCoordAttrib<TexCoordT>::clear_vertex_texcoords()
{
    for (VertexIter v_it = kernel_.vertices_begin(); v_it != kernel_.vertices_end(); ++v_it)
        vtexcoord_prop_[v_it->idx()] = default_texcoord_;
    vertex_texcoords_available_   = false;
}


} // Namespace OpenVolumeMesh
