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

#ifndef BASEENTITIES_HH_
#define BASEENTITIES_HH_

#include <vector>

#include "OpenVolumeMesh/Config/Export.hh"
#include "OpenVolumeMeshHandle.hh"

namespace OpenVolumeMesh {

class OVM_EXPORT OpenVolumeMeshEdge {
friend class TopologyKernel;
public:
    OpenVolumeMeshEdge(const VertexHandle& _fromVertex,
                       const VertexHandle& _toVertex) :
        fromVertex_(_fromVertex),
        toVertex_(_toVertex) {
    }

    const VertexHandle from_vertex() const {
        return fromVertex_;
    }
    const VertexHandle to_vertex() const {
        return toVertex_;
    }

protected:

    void set_from_vertex(const VertexHandle& _vertex) {
        fromVertex_ = _vertex;
    }
    void set_to_vertex(const VertexHandle& _vertex) {
        toVertex_ = _vertex;
    }

private:
    VertexHandle fromVertex_;
    VertexHandle toVertex_;
};

// Stream operator for edges
std::ostream& operator<<(std::ostream& _os, const OpenVolumeMeshEdge& _edge);

//***************************************************************************

class OVM_EXPORT OpenVolumeMeshFace {
friend class TopologyKernel;
public:
    explicit OpenVolumeMeshFace(const std::vector<HalfEdgeHandle>& _halfedges) :
        halfedges_(_halfedges) {
    }

    const std::vector<HalfEdgeHandle>& halfedges() const & {
        return halfedges_;
    }

    const std::vector<HalfEdgeHandle>& halfedges() const && = delete;
    std::vector<HalfEdgeHandle> halfedges() && {
        return std::move(halfedges_);
    }

protected:

    void set_halfedges(const std::vector<HalfEdgeHandle>& _halfedges) {
        halfedges_ = _halfedges;
    }

private:
    std::vector<HalfEdgeHandle> halfedges_;
};

// Stream operator for faces
std::ostream& operator<<(std::ostream& _os, const OpenVolumeMeshFace& _face);

//***************************************************************************

class OVM_EXPORT OpenVolumeMeshCell {
friend class TopologyKernel;
public:
    explicit OpenVolumeMeshCell(const std::vector<HalfFaceHandle>& _halffaces) :
        halffaces_(_halffaces) {
    }

    const std::vector<HalfFaceHandle>& halffaces() const & {
        return halffaces_;
    }

    const std::vector<HalfFaceHandle>& halffaces() const && = delete;
    std::vector<HalfFaceHandle> halffaces() && {
        return std::move(halffaces_);
    }

protected:

    void set_halffaces(const std::vector<HalfFaceHandle>& _halffaces) {
        halffaces_ = _halffaces;
    }

private:
    std::vector<HalfFaceHandle> halffaces_;
};

// Stream operator for cells
std::ostream& operator<<(std::ostream& _os, const OpenVolumeMeshCell& _cell);

} // Namespace OpenVolumeMesh

#endif /* BASEENTITIES_HH_ */
