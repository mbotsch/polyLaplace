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

#ifndef HEXAHEDRALMESHITERATORS_HH
#define HEXAHEDRALMESHITERATORS_HH

#include "../Core/Iterators.hh"
#include "OpenVolumeMesh/Config/Export.hh"

namespace OpenVolumeMesh {

class HexahedralMeshTopologyKernel;


class OVM_EXPORT CellSheetCellIter : public BaseCirculator<CellHandle, CellHandle> {
private:
    typedef BaseCirculator<CellHandle, CellHandle>    BaseIter;
public:
	CellSheetCellIter(const CellHandle& _ref_h, const unsigned char _orthDir,
            const HexahedralMeshTopologyKernel* _mesh, int _max_laps = 1);

	// Post increment/decrement operator
	CellSheetCellIter operator++(int) {
		CellSheetCellIter cpy = *this;
		++(*this);
		return cpy;
	}
	CellSheetCellIter operator--(int) {
		CellSheetCellIter cpy = *this;
		--(*this);
		return cpy;
	}
	CellSheetCellIter operator+(int _n) {
		CellSheetCellIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			++cpy;
		}
		return cpy;
	}
	CellSheetCellIter operator-(int _n) {
		CellSheetCellIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			--cpy;
		}
		return cpy;
	}
	CellSheetCellIter& operator+=(int _n) {
		for(int i = 0; i < _n; ++i) {
			++(*this);
		}
		return *this;
	}
	CellSheetCellIter& operator-=(int _n) {
		for(int i = 0; i < _n; ++i) {
			--(*this);
		}
		return *this;
	}

	CellSheetCellIter& operator++();
	CellSheetCellIter& operator--();

private:
    std::vector<CellHandle> neighb_sheet_cell_hs_;
    size_t cur_index_;
};


class OVM_EXPORT HalfFaceSheetHalfFaceIter : public BaseCirculator<HalfFaceHandle,HalfFaceHandle> {
private:
    typedef BaseCirculator<HalfFaceHandle, HalfFaceHandle> BaseIter;
public:
	HalfFaceSheetHalfFaceIter(const HalfFaceHandle& _ref_h,
            const HexahedralMeshTopologyKernel* _mesh, int _max_laps = 1);

	// Post increment/decrement operator
	HalfFaceSheetHalfFaceIter operator++(int) {
		HalfFaceSheetHalfFaceIter cpy = *this;
		++(*this);
		return cpy;
	}
	HalfFaceSheetHalfFaceIter operator--(int) {
		HalfFaceSheetHalfFaceIter cpy = *this;
		--(*this);
		return cpy;
	}
	HalfFaceSheetHalfFaceIter operator+(int _n) {
		HalfFaceSheetHalfFaceIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			++cpy;
		}
		return cpy;
	}
	HalfFaceSheetHalfFaceIter operator-(int _n) {
		HalfFaceSheetHalfFaceIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			--cpy;
		}
		return cpy;
	}
	HalfFaceSheetHalfFaceIter& operator+=(int _n) {
		for(int i = 0; i < _n; ++i) {
			++(*this);
		}
		return *this;
	}
	HalfFaceSheetHalfFaceIter& operator-=(int _n) {
		for(int i = 0; i < _n; ++i) {
			--(*this);
		}
		return *this;
	}

	HalfFaceSheetHalfFaceIter& operator++();
	HalfFaceSheetHalfFaceIter& operator--();

    const EdgeHandle& common_edge() const { return common_edges_[cur_index_]; }

private:
	std::vector<HalfFaceHandle> adjacent_halffaces_;
    std::vector<EdgeHandle> common_edges_;
    size_t cur_index_;
};

/** \brief Iterate over all vertices of a hexahedron in a specific order
 *
 * Vertices are addressed in the following order:
 *
 *      5-------6
 *     /|      /|
 *    / |     / |
 *   3-------2  |
 *   |  4----|--7
 *   | /     | /
 *   |/      |/
 *   0-------1
 */

class OVM_EXPORT HexVertexIter : public BaseCirculator<CellHandle,
    VertexHandle> {
private:
    typedef BaseCirculator<CellHandle,
            VertexHandle> BaseIter;
public:
    HexVertexIter(const CellHandle& _ref_h,
                  const HexahedralMeshTopologyKernel* _mesh,
                  int _max_laps = 1);

    // Post increment/decrement operator
    HexVertexIter operator++(int) {
        HexVertexIter cpy = *this;
        ++(*this);
        return cpy;
    }
    HexVertexIter operator--(int) {
        HexVertexIter cpy = *this;
        --(*this);
        return cpy;
    }
    HexVertexIter operator+(int _n) {
        HexVertexIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            ++cpy;
        }
        return cpy;
    }
    HexVertexIter operator-(int _n) {
        HexVertexIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            --cpy;
        }
        return cpy;
    }
    HexVertexIter& operator+=(int _n) {
        for(int i = 0; i < _n; ++i) {
            ++(*this);
        }
        return *this;
    }
    HexVertexIter& operator-=(int _n) {
        for(int i = 0; i < _n; ++i) {
            --(*this);
        }
        return *this;
    }

    HexVertexIter& operator++();
    HexVertexIter& operator--();

private:
    std::vector<VertexHandle> vertices_;
    size_t cur_index_;
};

} // Namespace OpenVolumeMesh

#endif /* HEXAHEDRALMESHITERATORS_HH */
