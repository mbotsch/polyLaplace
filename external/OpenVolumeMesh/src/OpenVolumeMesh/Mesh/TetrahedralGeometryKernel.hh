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

#ifndef TETRAHEDRALGEOMETRYKERNEL_HH_
#define TETRAHEDRALGEOMETRYKERNEL_HH_

#include <cassert>
#include <iostream>

#include "../Geometry/VectorT.hh"
#include "../Core/GeometryKernel.hh"
#include "TetrahedralMeshTopologyKernel.hh"

namespace OpenVolumeMesh {

template <class VecT, class TopologyKernelT = TetrahedralMeshTopologyKernel>
class TetrahedralGeometryKernel : public GeometryKernel<VecT, TopologyKernelT> {
public:

    typedef VecT PointT;
    typedef TopologyKernelT KernelT;
    typedef GeometryKernel<VecT, TopologyKernelT> ParentT;

    /// Constructor
    TetrahedralGeometryKernel() {}

    /// Destructor
    ~TetrahedralGeometryKernel() {}

    VertexHandle split_edge(HalfEdgeHandle heh, double alpha = 0.5)
    {
        OpenVolumeMeshEdge e = TopologyKernelT::halfedge(heh);
        PointT newPos = alpha*ParentT::vertex(e.from_vertex()) + (1.0-alpha)*ParentT::vertex(e.to_vertex());
        VertexHandle splitVertex = ParentT::add_vertex(newPos);
        TopologyKernelT::split_edge(heh, splitVertex);
        return splitVertex;
    }

    VertexHandle split_edge(EdgeHandle eh)
    {
        return split_edge(TopologyKernelT::halfedge_handle(eh,0));
    }

    VertexHandle split_face(FaceHandle fh, PointT pos)
    {
        VertexHandle splitVertex = ParentT::add_vertex(pos);
        TopologyKernelT::split_face(fh, splitVertex);
        return splitVertex;
    }

    VertexHandle split_face(FaceHandle fh)
    {
        VertexHandle splitVertex = ParentT::add_vertex(ParentT::barycenter(fh));
        TopologyKernelT::split_face(fh, splitVertex);
        return splitVertex;
    }

protected:

};

} // Namespace OpenVolumeMesh

#endif /* TETRAHEDRALGEOMETRYKERNEL_HH_ */
