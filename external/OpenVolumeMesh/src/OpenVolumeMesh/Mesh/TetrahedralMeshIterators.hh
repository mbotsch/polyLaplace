#pragma once
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

#include "../Core/Iterators.hh"
#include "OpenVolumeMesh/Config/Export.hh"

#include <array>

namespace OpenVolumeMesh {

class TetrahedralMeshTopologyKernel;


/** \brief Iterate over all vertices of a hexahedron in a specific order
 *
 * Vertices are addressed in the following order: vertices of one halfface in ccw order, then the remaining vertex
 *
 */

class OVM_EXPORT TetVertexIter : public BaseCirculator<CellHandle,
    VertexHandle> {
private:
    typedef BaseCirculator<CellHandle,
            VertexHandle> BaseIter;
public:
    TetVertexIter(const CellHandle& _ref_h,
                  const TetrahedralMeshTopologyKernel* _mesh,
                  int _max_laps = 1);

    // Post increment/decrement operator
    TetVertexIter operator++(int) {
        TetVertexIter cpy = *this;
        ++(*this);
        return cpy;
    }
    TetVertexIter operator--(int) {
        TetVertexIter cpy = *this;
        --(*this);
        return cpy;
    }
    TetVertexIter operator+(int _n) {
        TetVertexIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            ++cpy;
        }
        return cpy;
    }
    TetVertexIter operator-(int _n) {
        TetVertexIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            --cpy;
        }
        return cpy;
    }
    TetVertexIter& operator+=(int _n) {
        for(int i = 0; i < _n; ++i) {
            ++(*this);
        }
        return *this;
    }
    TetVertexIter& operator-=(int _n) {
        for(int i = 0; i < _n; ++i) {
            --(*this);
        }
        return *this;
    }

    TetVertexIter& operator++();
    TetVertexIter& operator--();

private:
    std::array<VertexHandle, 4> vertices_;
    size_t cur_index_;
};

} // Namespace OpenVolumeMesh
