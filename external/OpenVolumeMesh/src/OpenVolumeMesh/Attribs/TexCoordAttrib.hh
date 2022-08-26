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

#ifndef TEXCOORDATTRIB_HH_
#define TEXCOORDATTRIB_HH_

#include <cassert>

#include "../Core/OpenVolumeMeshHandle.hh"
#include "OpenVolumeMeshStatus.hh"
#include "../Core/PropertyDefines.hh"
#include "../Core/TopologyKernel.hh"

namespace OpenVolumeMesh {

//== CLASS DEF ================================================================

template <class TexCoordT>
class TexCoordAttrib {
public:

    TexCoordAttrib(TopologyKernel& _kernel, const TexCoordT _def = TexCoordT());

    virtual ~TexCoordAttrib();

    //==================
    // Vertices
    //==================
    const TexCoordT& operator[](const VertexHandle& _h) const {
        assert((unsigned int)_h.idx() < kernel_.n_vertices());
        return vtexcoord_prop_[_h];
    }

    TexCoordT& operator[](const VertexHandle& _h) {
        assert((unsigned int)_h.idx() < kernel_.n_vertices());
        vertex_texcoords_available_ = true;
        return vtexcoord_prop_[_h];
    }

    bool vertex_texcoords_available() const  { return vertex_texcoords_available_;   }

    void clear_vertex_texcoords();


private:

    VertexPropertyT<TexCoordT> vtexcoord_prop_;

    TopologyKernel& kernel_;

    bool vertex_texcoords_available_;

    TexCoordT default_texcoord_;

};

} // Namespace OpenVolumeMesh

#if defined(INCLUDE_TEMPLATES) && !defined(TEXCOORDATTRIBT_CC)
#include "TexCoordAttribT_impl.hh"
#endif

#endif /* TEXCOORDATTRIB_HH_ */
