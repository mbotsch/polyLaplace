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

#include <set>

#include "HexahedralMeshIterators.hh"
#include "HexahedralMeshTopologyKernel.hh"

namespace OpenVolumeMesh {

//================================================================================================
// CellSheetCellIter
//================================================================================================


CellSheetCellIter::CellSheetCellIter(const CellHandle& _ref_h,
        const unsigned char _orthDir, const HexahedralMeshTopologyKernel* _mesh, int _max_laps) :
BaseIter(_mesh, _ref_h, _max_laps) {

    if(!_mesh->has_face_bottom_up_incidences()) {
#ifndef NDEBUG
        std::cerr << "This iterator needs bottom-up incidences!" << std::endl;
#endif
        BaseIter::valid(false);
        return;
    }

	// First off, get all surrounding cells
	std::vector<HalfFaceHandle> halffaces = _mesh->cell(_ref_h).halffaces();
	for(std::vector<HalfFaceHandle>::const_iterator hf_it = halffaces.begin();
			hf_it != halffaces.end(); ++hf_it) {
		// Add those, that are perpendicular to the specified _orthDir
		if(_mesh->orientation(*hf_it, _ref_h) != _orthDir &&
				_mesh->orientation(*hf_it, _ref_h) != _mesh->opposite_orientation(_orthDir)) {
			CellHandle ch = _mesh->incident_cell(_mesh->opposite_halfface_handle(*hf_it));
			if(ch != TopologyKernel::InvalidCellHandle) {
                neighb_sheet_cell_hs_.push_back(ch);
			}
		}
	}

    // Remove all duplicate entries
    std::sort(neighb_sheet_cell_hs_.begin(), neighb_sheet_cell_hs_.end());
    neighb_sheet_cell_hs_.resize(std::unique(neighb_sheet_cell_hs_.begin(), neighb_sheet_cell_hs_.end()) - neighb_sheet_cell_hs_.begin());

    cur_index_ = 0;
    BaseIter::valid(neighb_sheet_cell_hs_.size() > 0);
	if(BaseIter::valid()) {
        BaseIter::cur_handle(neighb_sheet_cell_hs_[cur_index_]);
	}
}


CellSheetCellIter& CellSheetCellIter::operator--() {

    if (cur_index_ == 0) {
        cur_index_ = neighb_sheet_cell_hs_.size() - 1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    } else {
        --cur_index_;
    }

    BaseIter::cur_handle(neighb_sheet_cell_hs_[cur_index_]);

    return *this;
}


CellSheetCellIter& CellSheetCellIter::operator++() {

    ++cur_index_;
    if(cur_index_ == neighb_sheet_cell_hs_.size()) {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
	}

    BaseIter::cur_handle(neighb_sheet_cell_hs_[cur_index_]);

    return *this;
}

//================================================================================================
// HalfFaceSheetHalfFaceIter
//================================================================================================


HalfFaceSheetHalfFaceIter::HalfFaceSheetHalfFaceIter(const HalfFaceHandle& _ref_h,
        const HexahedralMeshTopologyKernel* _mesh, int _max_laps) :
BaseIter(_mesh, _ref_h, _max_laps) {

	if(!_mesh->has_face_bottom_up_incidences()) {
#ifndef NDEBUG
        std::cerr << "This iterator needs bottom-up incidences!" << std::endl;
#endif
        BaseIter::valid(false);
        return;
    }

	/*
	 * Each halfface uniquely belongs to either a cell
	 * or the boundary. If the halfface belongs
	 * to a cell, it suffices to determine the local axis
	 * the halfface represents w.r.t. the cell and to
	 * iterate over all neighboring cells orthogonal to
	 * this direction. We have to find those halffaces
	 * of the neighboring cells that contain exactly one
	 * of the initial halfface's opposite halfedges.
	 */

	if(_mesh->is_boundary(_ref_h)) {
#ifndef NDEBUG
		std::cerr << "HalfFaceSheetHalfFaceIter: HalfFace is boundary!" << std::endl;
#endif
		BaseIter::valid(false);
        return;
	}

	CellHandle ch = _mesh->incident_cell(_ref_h);
	unsigned char orientation = _mesh->orientation(_ref_h, ch);
	std::vector<HalfEdgeHandle> hes_v = _mesh->opposite_halfface(_mesh->halfface(_ref_h)).halfedges();
	std::set<HalfEdgeHandle> hes;
	hes.insert(hes_v.begin(), hes_v.end());

	for(CellSheetCellIter csc_it = _mesh->csc_iter(ch, orientation);
			csc_it.valid(); ++csc_it) {

		std::vector<HalfFaceHandle> hfs = _mesh->cell(*csc_it).halffaces();
		for(std::vector<HalfFaceHandle>::const_iterator hf_it = hfs.begin();
				hf_it != hfs.end(); ++hf_it) {

			std::vector<HalfEdgeHandle> hf_hes = _mesh->halfface(*hf_it).halfedges();
			for(std::vector<HalfEdgeHandle>::const_iterator he_it = hf_hes.begin();
					he_it != hf_hes.end(); ++he_it) {

				if(hes.count(*he_it) > 0) {
					// Found halfface that lies on the same sheet
					adjacent_halffaces_.push_back(*hf_it);
					common_edges_.push_back(_mesh->edge_handle(*he_it));
					break;
				}
			}
		}
	}

    cur_index_ = 0;
    BaseIter::valid(adjacent_halffaces_.size() > 0);
    if(BaseIter::valid()) {
        BaseIter::cur_handle(adjacent_halffaces_[cur_index_]);
    }
}


HalfFaceSheetHalfFaceIter& HalfFaceSheetHalfFaceIter::operator--() {

    if (cur_index_ == 0) {
        cur_index_ = adjacent_halffaces_.size() - 1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    } else {
        --cur_index_;
    }

    BaseIter::cur_handle(adjacent_halffaces_[cur_index_]);

    return *this;
}


HalfFaceSheetHalfFaceIter& HalfFaceSheetHalfFaceIter::operator++() {

    ++cur_index_;
    if(cur_index_ == adjacent_halffaces_.size()) {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }

    BaseIter::cur_handle(adjacent_halffaces_[cur_index_]);

    return *this;
}

//================================================================================================
// HexVertexIter
//================================================================================================


HexVertexIter::HexVertexIter(const CellHandle& _ref_h,
        const HexahedralMeshTopologyKernel* _mesh, int _max_laps) :
BaseIter(_mesh, _ref_h, _max_laps) {

    assert(_ref_h.is_valid());

    HexahedralMeshTopologyKernel::Cell cell = _mesh->cell(_ref_h);
    assert(cell.halffaces().size() == 6);

    // Get first half-face
    HalfFaceHandle curHF = *cell.halffaces().begin();
    assert(curHF.is_valid());

    // Get first half-edge
    assert(_mesh->halfface(curHF).halfedges().size() == 4);
    HalfEdgeHandle curHE = *_mesh->halfface(curHF).halfedges().begin();
    assert(curHE.is_valid());

    vertices_.push_back(_mesh->halfedge(curHE).from_vertex());

    curHE = _mesh->prev_halfedge_in_halfface(curHE, curHF);

    vertices_.push_back(_mesh->halfedge(curHE).from_vertex());

    curHE = _mesh->prev_halfedge_in_halfface(curHE, curHF);

    vertices_.push_back(_mesh->halfedge(curHE).from_vertex());

    curHE = _mesh->prev_halfedge_in_halfface(curHE, curHF);

    vertices_.push_back(_mesh->halfedge(curHE).from_vertex());

    curHE = _mesh->prev_halfedge_in_halfface(curHE, curHF);
    curHF = _mesh->adjacent_halfface_in_cell(curHF, curHE);
    curHE = _mesh->opposite_halfedge_handle(curHE);
    curHE = _mesh->next_halfedge_in_halfface(curHE, curHF);
    curHE = _mesh->next_halfedge_in_halfface(curHE, curHF);
    curHF = _mesh->adjacent_halfface_in_cell(curHF, curHE);
    curHE = _mesh->opposite_halfedge_handle(curHE);

    vertices_.push_back(_mesh->halfedge(curHE).to_vertex());

    curHE = _mesh->prev_halfedge_in_halfface(curHE, curHF);

    vertices_.push_back(_mesh->halfedge(curHE).to_vertex());

    curHE = _mesh->prev_halfedge_in_halfface(curHE, curHF);

    vertices_.push_back(_mesh->halfedge(curHE).to_vertex());

    vertices_.push_back(_mesh->halfedge(curHE).from_vertex());

    cur_index_ = 0;
    BaseIter::valid(vertices_.size() > 0);
    if(BaseIter::valid()) {
        BaseIter::cur_handle(vertices_[cur_index_]);
    }
}


HexVertexIter& HexVertexIter::operator--() {

    if (cur_index_ == 0) {
        cur_index_ = vertices_.size() - 1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    } else {
        --cur_index_;
    }

    BaseIter::cur_handle(vertices_[cur_index_]);

    return *this;
}


HexVertexIter& HexVertexIter::operator++() {

    ++cur_index_;
    if(cur_index_ == vertices_.size()) {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }

    BaseIter::cur_handle(vertices_[cur_index_]);

    return *this;
}

} // Namespace OpenVolumeMesh
