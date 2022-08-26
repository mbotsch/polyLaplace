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

#include <algorithm>
#include <set>

#include "Iterators.hh"
#include "TopologyKernel.hh"

namespace OpenVolumeMesh {

//================================================================================================
// VertexOHalfEdgeIter
//================================================================================================


VertexOHalfEdgeIter::VertexOHalfEdgeIter(const VertexHandle& _ref_h,
        const TopologyKernel* _mesh, int _max_laps) :
BaseIter(_mesh, _ref_h, _max_laps),
cur_index_(0) {

  if(!_mesh->has_vertex_bottom_up_incidences()) {
#ifndef NDEBUG
        std::cerr << "This iterator needs bottom-up incidences!" << std::endl;
#endif
        BaseIter::valid(false);
        return;
    }

  if((unsigned int)_ref_h.idx() >= BaseIter::mesh()->outgoing_hes_per_vertex_.size()) {
    BaseIter::valid(false);
  }

  if(BaseIter::valid()) {
    if((unsigned int)cur_index_ >= BaseIter::mesh()->outgoing_hes_per_vertex_[_ref_h.idx()].size()) {
      BaseIter::valid(false);
    }
  }

  if(BaseIter::valid()) {
    BaseIter::cur_handle((
        BaseIter::mesh()->outgoing_hes_per_vertex_[_ref_h.idx()])[cur_index_]);
  }
}


VertexOHalfEdgeIter& VertexOHalfEdgeIter::operator--() {

    size_t n_outgoing_halfedges = BaseIter::mesh()->outgoing_hes_per_vertex_[BaseIter::ref_handle().idx()].size();

    if (cur_index_ == 0) {
        cur_index_ = n_outgoing_halfedges-1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else {
        --cur_index_;
    }

    BaseIter::cur_handle((BaseIter::mesh()->outgoing_hes_per_vertex_[BaseIter::ref_handle().idx()])[cur_index_]);

  return *this;
}


VertexOHalfEdgeIter& VertexOHalfEdgeIter::operator++() {

    size_t n_outgoing_halfedges = BaseIter::mesh()->outgoing_hes_per_vertex_[BaseIter::ref_handle().idx()].size();

    ++cur_index_;

    if (cur_index_ == n_outgoing_halfedges) {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }

    BaseIter::cur_handle((BaseIter::mesh()->outgoing_hes_per_vertex_[BaseIter::ref_handle().idx()])[cur_index_]);

  return *this;
}


//================================================================================================
// VertexVertexIter
//================================================================================================


VertexVertexIter::VertexVertexIter(const VertexHandle& _ref_h,
        const TopologyKernel* _mesh, int _max_laps) :
BaseIter(_mesh, _ref_h, _max_laps),
cur_index_(0) {

  if(!_mesh->has_vertex_bottom_up_incidences()) {
#ifndef NDEBUG
        std::cerr << "This iterator needs bottom-up incidences!" << std::endl;
#endif
        BaseIter::valid(false);
        return;
    }

  if((size_t)_ref_h.idx() >= BaseIter::mesh()->outgoing_hes_per_vertex_.size()) {
    BaseIter::valid(false);
  }

  if(BaseIter::valid()) {
    if((size_t)cur_index_ >= BaseIter::mesh()->outgoing_hes_per_vertex_[_ref_h.idx()].size()) {
      BaseIter::valid(false);
    }
  }

  if(BaseIter::valid()) {
    HalfEdgeHandle heh = BaseIter::mesh()->outgoing_hes_per_vertex_[_ref_h.idx()][cur_index_];
    BaseIter::cur_handle(BaseIter::mesh()->halfedge(heh).to_vertex());
  }
}


VertexVertexIter& VertexVertexIter::operator--() {

    size_t n_outgoing_halfedges = BaseIter::mesh()->outgoing_hes_per_vertex_[BaseIter::ref_handle().idx()].size();

    if (cur_index_ == 0) {
        cur_index_ = n_outgoing_halfedges-1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else {
        --cur_index_;
    }

    HalfEdgeHandle heh = BaseIter::mesh()->outgoing_hes_per_vertex_[BaseIter::ref_handle().idx()][cur_index_];
    BaseIter::cur_handle(BaseIter::mesh()->halfedge(heh).to_vertex());

  return *this;
}


VertexVertexIter& VertexVertexIter::operator++() {

    size_t n_outgoing_halfedges = BaseIter::mesh()->outgoing_hes_per_vertex_[BaseIter::ref_handle().idx()].size();

    ++cur_index_;

    if (cur_index_ == n_outgoing_halfedges) {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }


    HalfEdgeHandle heh = BaseIter::mesh()->outgoing_hes_per_vertex_[BaseIter::ref_handle().idx()][cur_index_];
    BaseIter::cur_handle(BaseIter::mesh()->halfedge(heh).to_vertex());

  return *this;
}


////================================================================================================
//// HalfEdgeHalfFaceIter
////================================================================================================


HalfEdgeHalfFaceIter::HalfEdgeHalfFaceIter(const HalfEdgeHandle& _ref_h,
        const TopologyKernel* _mesh, int _max_laps) :
BaseIter(_mesh, _ref_h, _max_laps),
cur_index_(0) {

	if(!_mesh->has_edge_bottom_up_incidences()) {
#ifndef NDEBUG
        std::cerr << "This iterator needs bottom-up incidences!" << std::endl;
#endif
        BaseIter::valid(false);
        return;
    }

	if((unsigned int)_ref_h.idx() >= BaseIter::mesh()->incident_hfs_per_he_.size()) {
		BaseIter::valid(false);
	}

	if(BaseIter::valid()) {
		if((unsigned int)cur_index_ >= BaseIter::mesh()->incident_hfs_per_he_[_ref_h.idx()].size()) {
			BaseIter::valid(false);
		}
	}

	if(BaseIter::valid()) {
		BaseIter::cur_handle((
				BaseIter::mesh()->incident_hfs_per_he_[_ref_h.idx()])[cur_index_]);
	}
}


HalfEdgeHalfFaceIter& HalfEdgeHalfFaceIter::operator--() {

    size_t n_outgoing_halffaces = BaseIter::mesh()->incident_hfs_per_he_[BaseIter::ref_handle().idx()].size();

    if (cur_index_ == 0) {
        cur_index_ = n_outgoing_halffaces-1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else {
        --cur_index_;
    }

    BaseIter::cur_handle((BaseIter::mesh()->incident_hfs_per_he_[BaseIter::ref_handle().idx()])[cur_index_]);

    return *this;
}


HalfEdgeHalfFaceIter& HalfEdgeHalfFaceIter::operator++() {


    size_t n_outgoing_halffaces = BaseIter::mesh()->incident_hfs_per_he_[BaseIter::ref_handle().idx()].size();

    ++cur_index_;

    if (cur_index_ == n_outgoing_halffaces) {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }

    BaseIter::cur_handle((BaseIter::mesh()->incident_hfs_per_he_[BaseIter::ref_handle().idx()])[cur_index_]);

    return *this;
}

////================================================================================================
//// VertexFaceIter
////================================================================================================

VertexFaceIter::VertexFaceIter(const VertexHandle& _ref_h,
        const TopologyKernel* _mesh, int _max_laps) :
BaseIter(_mesh, _ref_h, _max_laps) {

  if(!_mesh->has_full_bottom_up_incidences()) {
#ifndef NDEBUG
        std::cerr << "This iterator needs bottom-up incidences!" << std::endl;
#endif
        BaseIter::valid(false);
        return;
    }

    if((unsigned int)_ref_h.idx() >= BaseIter::mesh()->outgoing_hes_per_vertex_.size()) {
        BaseIter::valid(false);
        return;
    }

    // Build up face list
    const std::vector<HalfEdgeHandle>& incidentHalfedges = BaseIter::mesh()->outgoing_hes_per_vertex_[_ref_h.idx()];
    for(std::vector<HalfEdgeHandle>::const_iterator it = incidentHalfedges.begin(); it != incidentHalfedges.end(); ++it) {

        if(*it < 0 || (unsigned int)it->idx() >= BaseIter::mesh()->incident_hfs_per_he_.size()) continue;
            const std::vector<HalfFaceHandle>& incidentHalfFaces = BaseIter::mesh()->incident_hfs_per_he_[it->idx()];

        for (std::vector<HalfFaceHandle>::const_iterator hf_it = incidentHalfFaces.begin();
                hf_it != incidentHalfFaces.end(); ++hf_it) {
            faces_.push_back(BaseIter::mesh()->face_handle(*hf_it));
        }
    }
    // Remove all duplicate entries
    std::sort(faces_.begin(), faces_.end());
    faces_.resize(std::unique(faces_.begin(), faces_.end()) - faces_.begin());

    // Remove invalid handles
    if ((faces_.size() > 0) && !faces_.front().is_valid())
        faces_.erase(faces_.begin());

    cur_index_ = 0;
    BaseIter::valid(faces_.size() > 0);
    if(BaseIter::valid()) {
        BaseIter::cur_handle(faces_[cur_index_]);
    }
}

VertexFaceIter& VertexFaceIter::operator--() {

    if(cur_index_ == 0) {
        cur_index_ = faces_.size()-1;
        --lap_;
        if (lap_ < 0) {
            BaseIter::valid(false);
        }
    } else {
        --cur_index_;
    }

    BaseIter::cur_handle(faces_[cur_index_]);
    return *this;
}


VertexFaceIter& VertexFaceIter::operator++() {

    ++cur_index_;
    if(cur_index_ == faces_.size()) {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_) {
            BaseIter::valid(false);
        }
  }
    BaseIter::cur_handle(faces_[cur_index_]);
  return *this;
}

////================================================================================================
//// VertexCellIter
////================================================================================================

VertexCellIter::VertexCellIter(const VertexHandle& _ref_h,
        const TopologyKernel* _mesh, int _max_laps) :
BaseIter(_mesh, _ref_h, _max_laps) {

	if(!_mesh->has_full_bottom_up_incidences()) {
#ifndef NDEBUG
        std::cerr << "This iterator needs bottom-up incidences!" << std::endl;
#endif
        BaseIter::valid(false);
        return;
    }

    if((unsigned int)_ref_h.idx() >= BaseIter::mesh()->outgoing_hes_per_vertex_.size()) {
        BaseIter::valid(false);
        return;
    }

    // Build up cell list
    const std::vector<HalfEdgeHandle>& incidentHalfedges = BaseIter::mesh()->outgoing_hes_per_vertex_[_ref_h.idx()];
    for(std::vector<HalfEdgeHandle>::const_iterator it = incidentHalfedges.begin(); it != incidentHalfedges.end(); ++it) {

    	if(*it < 0 || (unsigned int)it->idx() >= BaseIter::mesh()->incident_hfs_per_he_.size()) continue;
        const std::vector<HalfFaceHandle>& incidentHalfFaces = BaseIter::mesh()->incident_hfs_per_he_[it->idx()];

    	for(std::vector<HalfFaceHandle>::const_iterator hf_it = incidentHalfFaces.begin();
                hf_it != incidentHalfFaces.end(); ++hf_it) {
    		if((unsigned int)hf_it->idx() < BaseIter::mesh()->incident_cell_per_hf_.size()) {
    			CellHandle c_idx = BaseIter::mesh()->incident_cell_per_hf_[hf_it->idx()];
    			if(c_idx != TopologyKernel::InvalidCellHandle)
                    cells_.push_back(c_idx);
    		}
    	}
    }
    // Remove all duplicate entries
    std::sort(cells_.begin(), cells_.end());
    cells_.resize(std::unique(cells_.begin(), cells_.end()) - cells_.begin());

    // Remove invalid handles
    if ((cells_.size() > 0) && !cells_.front().is_valid())
        cells_.erase(cells_.begin());

    cur_index_ = 0;
    BaseIter::valid(cells_.size()>0);
    if(BaseIter::valid()) {
        BaseIter::cur_handle(cells_[cur_index_]);
    }
}

VertexCellIter& VertexCellIter::operator--() {

    if(cur_index_ == 0) {
        cur_index_ = cells_.size()-1;
        --lap_;
        if (lap_ < 0) {
            BaseIter::valid(false);
        }
    } else {
        --cur_index_;
    }

    BaseIter::cur_handle(cells_[cur_index_]);
    return *this;
}


VertexCellIter& VertexCellIter::operator++() {

    ++cur_index_;
    if(cur_index_ == cells_.size()) {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_) {
            BaseIter::valid(false);
        }
	}
    BaseIter::cur_handle(cells_[cur_index_]);
	return *this;
}

////================================================================================================
//// HalfEdgeCellIter
////================================================================================================


HalfEdgeCellIter::HalfEdgeCellIter(const HalfEdgeHandle& _ref_h,
        const TopologyKernel* _mesh, int _max_laps) :
BaseIter(_mesh, _ref_h, _max_laps),
cur_index_(0) {

	if(!_mesh->has_edge_bottom_up_incidences() || !_mesh->has_face_bottom_up_incidences()) {
#ifndef NDEBUG
        std::cerr << "This iterator needs bottom-up incidences!" << std::endl;
#endif
        BaseIter::valid(false);
        return;
    }

    if((unsigned int)_ref_h.idx() >= BaseIter::mesh()->incident_hfs_per_he_.size()) {

        BaseIter::valid(false);
        return;
    }
    if((unsigned int)cur_index_ >= BaseIter::mesh()->incident_hfs_per_he_[_ref_h.idx()].size()) {

    	BaseIter::valid(false);
    	return;
    }
    if((unsigned int)((BaseIter::mesh()->incident_hfs_per_he_[_ref_h.idx()])[cur_index_]).idx() >=
    		BaseIter::mesh()->incident_cell_per_hf_.size()) {

    	BaseIter::valid(false);
    	return;
    }

    // collect cell handles
    const std::vector<HalfFaceHandle>& incidentHalffaces = BaseIter::mesh()->incident_hfs_per_he_[_ref_h.idx()];
    std::set<CellHandle> cells;
    for (unsigned int i = 0; i < incidentHalffaces.size(); ++i)
    {
        CellHandle ch = getCellHandle(i);
        if (ch.is_valid()) {
            if(cells.count(ch) == 0) {
                cells_.push_back(ch);
            }
            cells.insert(ch);
        }
    }

    BaseIter::valid(cells_.size() > 0);

    if(BaseIter::valid()) {
        BaseIter::cur_handle(cells_[cur_index_]);
    }
}


HalfEdgeCellIter& HalfEdgeCellIter::operator--() {

    if (cur_index_ == 0) {
        cur_index_ = cells_.size()-1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else {
        --cur_index_;
    }

    BaseIter::cur_handle(cells_[cur_index_]);
    return *this;
}


HalfEdgeCellIter& HalfEdgeCellIter::operator++() {

    ++cur_index_;

    if (cur_index_ == cells_.size()) {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }

    BaseIter::cur_handle(cells_[cur_index_]);
    return *this;
}

CellHandle HalfEdgeCellIter::getCellHandle(int _cur_index) const
{
    const std::vector<HalfFaceHandle>& halffacehandles = BaseIter::mesh()->incident_hfs_per_he_[BaseIter::ref_handle().idx()];
    HalfFaceHandle currentHalfface = halffacehandles[_cur_index];
    if(!currentHalfface.is_valid()) return CellHandle(-1);
    CellHandle cellhandle = BaseIter::mesh()->incident_cell_per_hf_[currentHalfface.idx()];
    return cellhandle;
}


////================================================================================================
//// CellVertexIter
////================================================================================================


CellVertexIter::CellVertexIter(const CellHandle& _ref_h,
        const TopologyKernel* _mesh, int _max_laps) :
BaseIter(_mesh, _ref_h, _max_laps) {

    OpenVolumeMeshCell c = BaseIter::mesh()->cell(_ref_h);
    std::vector<HalfFaceHandle>::const_iterator hf_iter = c.halffaces().begin();
    for(; hf_iter != c.halffaces().end(); ++hf_iter) {
        const OpenVolumeMeshFace& halfface = BaseIter::mesh()->halfface(*hf_iter);
        const std::vector<HalfEdgeHandle>& hes = halfface.halfedges();
        for(std::vector<HalfEdgeHandle>::const_iterator he_iter = hes.begin(); he_iter != hes.end(); ++he_iter) {
            incident_vertices_.push_back(BaseIter::mesh()->halfedge(*he_iter).to_vertex());
        }
    }

    // Remove all duplicate entries
    std::sort(incident_vertices_.begin(), incident_vertices_.end());
    incident_vertices_.resize(std::unique(incident_vertices_.begin(), incident_vertices_.end()) - incident_vertices_.begin());

    cur_index_ = 0;
    BaseIter::valid(incident_vertices_.size() > 0);

    if(BaseIter::valid()) {
        BaseIter::cur_handle(incident_vertices_[cur_index_]);
    }
}


CellVertexIter& CellVertexIter::operator--() {

    if (cur_index_ == 0) {
        cur_index_ = incident_vertices_.size()-1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else {
        --cur_index_;
    }

    BaseIter::cur_handle(incident_vertices_[cur_index_]);
    return *this;
}


CellVertexIter& CellVertexIter::operator++() {

    ++cur_index_;
    if (cur_index_ == incident_vertices_.size()){
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }

    BaseIter::cur_handle(incident_vertices_[cur_index_]);
	return *this;
}

////================================================================================================
//// CellCellIter
////================================================================================================


CellCellIter::CellCellIter(const CellHandle& _ref_h,
        const TopologyKernel* _mesh, int _max_laps) :
BaseIter(_mesh, _ref_h, _max_laps) {

    if(!_mesh->has_face_bottom_up_incidences()) {
#ifndef NDEBUG
        std::cerr << "This iterator needs bottom-up incidences!" << std::endl;
#endif
        BaseIter::valid(false);
        return;
    }

	std::vector<HalfFaceHandle>::const_iterator hf_iter = BaseIter::mesh()->cell(_ref_h).halffaces().begin();
    std::vector<HalfFaceHandle>::const_iterator hf_end  = BaseIter::mesh()->cell(_ref_h).halffaces().end();
    for(; hf_iter != hf_end; ++hf_iter) {

		HalfFaceHandle opp_hf = BaseIter::mesh()->opposite_halfface_handle(*hf_iter);
		CellHandle ch = BaseIter::mesh()->incident_cell_per_hf_[opp_hf.idx()];
		if(ch != TopologyKernel::InvalidCellHandle) {
            adjacent_cells_.push_back(ch);
		}
	}

    // Remove all duplicate entries
    std::sort(adjacent_cells_.begin(), adjacent_cells_.end());
    adjacent_cells_.resize(std::unique(adjacent_cells_.begin(), adjacent_cells_.end()) - adjacent_cells_.begin());

    cur_index_ = 0;
    BaseIter::valid(adjacent_cells_.size()>0);
	if(BaseIter::valid()) {
        BaseIter::cur_handle(adjacent_cells_[cur_index_]);
	}
}


CellCellIter& CellCellIter::operator--() {

    if (cur_index_ == 0) {
        cur_index_  = adjacent_cells_.size() - 1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else {
        --cur_index_;
    }

    BaseIter::cur_handle(adjacent_cells_[cur_index_]);
    return *this;
}


CellCellIter& CellCellIter::operator++() {

    ++cur_index_;
    if (cur_index_ == adjacent_cells_.size())
    {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }
    BaseIter::cur_handle(adjacent_cells_[cur_index_]);
    return *this;
}

////================================================================================================
//// HalfFaceVertexIter
////================================================================================================


HalfFaceVertexIter::HalfFaceVertexIter(const HalfFaceHandle& _ref_h,
        const TopologyKernel* _mesh, int _max_laps) :
BaseIter(_mesh, _ref_h, _max_laps) {

    if(!_ref_h.is_valid()) return;

    const OpenVolumeMeshFace& halfface = _mesh->halfface(_ref_h);
    const std::vector<HalfEdgeHandle>& hes = halfface.halfedges();
    for(std::vector<HalfEdgeHandle>::const_iterator he_it = hes.begin();
            he_it != hes.end(); ++he_it) {
        vertices_.push_back(_mesh->halfedge(*he_it).from_vertex());
    }

    cur_index_ = 0;

    BaseIter::valid(vertices_.size() > 0);
    if(BaseIter::valid()) {
        BaseIter::cur_handle(vertices_[cur_index_]);
    }
}


HalfFaceVertexIter& HalfFaceVertexIter::operator--() {

    if (cur_index_ == 0) {
        cur_index_  = vertices_.size() - 1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else {
        --cur_index_;
    }

    BaseIter::cur_handle(vertices_[cur_index_]);
    return *this;
}


HalfFaceVertexIter& HalfFaceVertexIter::operator++() {

    ++cur_index_;
    if (cur_index_ == vertices_.size())
    {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }
    BaseIter::cur_handle(vertices_[cur_index_]);
    return *this;
}

//================================================================================================
// BoundaryHalfFaceHalfFaceIter
//================================================================================================

BoundaryHalfFaceHalfFaceIter::BoundaryHalfFaceHalfFaceIter(const HalfFaceHandle& _ref_h,
        const TopologyKernel* _mesh, int _max_laps) :
BaseIter(_mesh, _ref_h, _max_laps) {

    if(!_mesh->has_face_bottom_up_incidences()) {
#ifndef NDEBUG
        std::cerr << "This iterator needs bottom-up incidences!" << std::endl;
#endif
        BaseIter::valid(false);
        return;
    }

    // Go over all incident halfedges
//    const std::vector<HalfEdgeHandle> halfedges = _mesh->halfface(_ref_h).halfedges();
    const OpenVolumeMeshFace& halfface = _mesh->halfface(_ref_h);
    const std::vector<HalfEdgeHandle>& halfedges = halfface.halfedges();
    for(std::vector<HalfEdgeHandle>::const_iterator he_it = halfedges.begin();
            he_it != halfedges.end(); ++he_it) {

        // Get outside halffaces
        OpenVolumeMesh::HalfEdgeHalfFaceIter hehf_it = _mesh->hehf_iter(_mesh->opposite_halfedge_handle(*he_it));
        for(; hehf_it.valid(); ++hehf_it) {

            if(_mesh->is_boundary(*hehf_it)) {
                neighbor_halffaces_.push_back(*hehf_it);
                common_edges_.push_back(_mesh->edge_handle(*he_it));
            }
        }
    }

    cur_index_ = 0;
    BaseIter::valid(neighbor_halffaces_.size() > 0);
    if(BaseIter::valid()) {
        BaseIter::cur_handle(neighbor_halffaces_[cur_index_]);
    }
}


BoundaryHalfFaceHalfFaceIter& BoundaryHalfFaceHalfFaceIter::operator--() {

    if (cur_index_ == 0)
    {
        cur_index_ = neighbor_halffaces_.size() - 1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else{
        --cur_index_;
    }

    BaseIter::cur_handle(neighbor_halffaces_[cur_index_]);
    return *this;
}


BoundaryHalfFaceHalfFaceIter& BoundaryHalfFaceHalfFaceIter::operator++() {

    ++cur_index_;
    if (cur_index_ == neighbor_halffaces_.size())
    {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }

    BaseIter::cur_handle(neighbor_halffaces_[cur_index_]);
    return *this;
}

////================================================================================================
//// VertexIter
////================================================================================================


VertexIter::VertexIter(const TopologyKernel* _mesh, const VertexHandle& _vh) :
BaseIter(_mesh, _vh),
cur_index_(_vh.idx()) {

    while ((unsigned int)cur_index_ < BaseIter::mesh()->n_vertices() && BaseIter::mesh()->is_deleted(VertexHandle(cur_index_)))
        ++cur_index_;
    if((unsigned int)cur_index_ >= BaseIter::mesh()->n_vertices()) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(VertexHandle(cur_index_));
}


VertexIter& VertexIter::operator--() {

    --cur_index_;
    while (cur_index_ >= 0 && BaseIter::mesh()->is_deleted(VertexHandle(cur_index_)))
        --cur_index_;
    if(cur_index_ < 0) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(VertexHandle(cur_index_));
    return *this;
}


VertexIter& VertexIter::operator++() {

	++cur_index_;
    while ((unsigned int)cur_index_ < BaseIter::mesh()->n_vertices() && BaseIter::mesh()->is_deleted(VertexHandle(cur_index_)))
        ++cur_index_;
    if((unsigned int)cur_index_ >= BaseIter::mesh()->n_vertices()) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(VertexHandle(cur_index_));
	return *this;
}

////================================================================================================
//// EdgeIter
////================================================================================================


EdgeIter::EdgeIter(const TopologyKernel* _mesh, const EdgeHandle& _eh) :
BaseIter(_mesh, _eh),
cur_index_(_eh.idx()) {

    while ((unsigned int)cur_index_ < BaseIter::mesh()->edges_.size() && BaseIter::mesh()->is_deleted(EdgeHandle(cur_index_)))
        ++cur_index_;
    if((unsigned int)cur_index_ >= BaseIter::mesh()->edges_.size()) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(EdgeHandle(cur_index_));
}


EdgeIter& EdgeIter::operator--() {

    --cur_index_;
    while (cur_index_ >= 0 && BaseIter::mesh()->is_deleted(EdgeHandle(cur_index_)))
        --cur_index_;
    if(cur_index_ < 0) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(EdgeHandle(cur_index_));
    return *this;
}


EdgeIter& EdgeIter::operator++() {

    ++cur_index_;
    while ((unsigned int)cur_index_ < BaseIter::mesh()->edges_.size() && BaseIter::mesh()->is_deleted(EdgeHandle(cur_index_)))
        ++cur_index_;
    if((unsigned int)cur_index_ >= BaseIter::mesh()->edges_.size()) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(EdgeHandle(cur_index_));
    return *this;
}

////================================================================================================
//// HalfEdgeIter
////================================================================================================


HalfEdgeIter::HalfEdgeIter(const TopologyKernel* _mesh, const HalfEdgeHandle& _heh) :
BaseIter(_mesh, _heh),
cur_index_(_heh.idx()) {

    while ((unsigned int)cur_index_ < BaseIter::mesh()->edges_.size() * 2 && BaseIter::mesh()->is_deleted(HalfEdgeHandle(cur_index_)))
        ++cur_index_;
    if((unsigned int)cur_index_ >= BaseIter::mesh()->edges_.size() * 2) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(HalfEdgeHandle(cur_index_));
}


HalfEdgeIter& HalfEdgeIter::operator--() {

    --cur_index_;
    while (cur_index_ >= 0 && BaseIter::mesh()->is_deleted(HalfEdgeHandle(cur_index_)))
        --cur_index_;
    if(cur_index_ < 0) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(HalfEdgeHandle(cur_index_));
    return *this;
}


HalfEdgeIter& HalfEdgeIter::operator++() {

    ++cur_index_;
    while ((unsigned int)cur_index_ < BaseIter::mesh()->edges_.size() * 2 && BaseIter::mesh()->is_deleted(HalfEdgeHandle(cur_index_)))
        ++cur_index_;
    if((unsigned int)cur_index_ >= BaseIter::mesh()->edges_.size() * 2) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(HalfEdgeHandle(cur_index_));
    return *this;
}

////================================================================================================
//// FaceIter
////================================================================================================


FaceIter::FaceIter(const TopologyKernel* _mesh, const FaceHandle& _fh) :
BaseIter(_mesh, _fh),
cur_index_(_fh.idx()) {

    while ((unsigned int)cur_index_ < BaseIter::mesh()->faces_.size() && BaseIter::mesh()->is_deleted(FaceHandle(cur_index_)))
        ++cur_index_;
    if((unsigned int)cur_index_ >= BaseIter::mesh()->faces_.size()) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(FaceHandle(cur_index_));
}


FaceIter& FaceIter::operator--() {

    --cur_index_;
    while (cur_index_ >= 0 && BaseIter::mesh()->is_deleted(FaceHandle(cur_index_)))
        --cur_index_;
    if(cur_index_ < 0) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(FaceHandle(cur_index_));
    return *this;
}


FaceIter& FaceIter::operator++() {

    ++cur_index_;
    while ((unsigned int)cur_index_ < BaseIter::mesh()->faces_.size() && BaseIter::mesh()->is_deleted(FaceHandle(cur_index_)))
        ++cur_index_;
    if((unsigned int)cur_index_ >= BaseIter::mesh()->faces_.size()) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(FaceHandle(cur_index_));
    return *this;
}

////================================================================================================
//// HalfFaceIter
////================================================================================================


HalfFaceIter::HalfFaceIter(const TopologyKernel* _mesh, const HalfFaceHandle& _hfh) :
BaseIter(_mesh, _hfh),
cur_index_(_hfh.idx()) {

    while ((unsigned int)cur_index_ < BaseIter::mesh()->faces_.size() * 2 && BaseIter::mesh()->is_deleted(HalfFaceHandle(cur_index_)))
        ++cur_index_;
    if((unsigned int)cur_index_ >= BaseIter::mesh()->faces_.size() * 2) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(HalfFaceHandle(cur_index_));
}


HalfFaceIter& HalfFaceIter::operator--() {

    --cur_index_;
    while (cur_index_ >= 0 && BaseIter::mesh()->is_deleted(HalfFaceHandle(cur_index_)))
        --cur_index_;
    if(cur_index_ < 0) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(HalfFaceHandle(cur_index_));
    return *this;
}


HalfFaceIter& HalfFaceIter::operator++() {

    ++cur_index_;
    while ((unsigned int)cur_index_ < BaseIter::mesh()->faces_.size() * 2 && BaseIter::mesh()->is_deleted(HalfFaceHandle(cur_index_)))
        ++cur_index_;
    if((unsigned int)cur_index_ >= BaseIter::mesh()->faces_.size() * 2) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(HalfFaceHandle(cur_index_));
    return *this;
}

////================================================================================================
//// CellIter
////================================================================================================


CellIter::CellIter(const TopologyKernel* _mesh, const CellHandle& _ch) :
BaseIter(_mesh, _ch),
cur_index_(_ch.idx()) {

    while ((unsigned int)cur_index_ < BaseIter::mesh()->cells_.size() && BaseIter::mesh()->is_deleted(CellHandle(cur_index_)))
        ++cur_index_;
    if((unsigned int)cur_index_ >= BaseIter::mesh()->cells_.size()) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(CellHandle(cur_index_));
}


CellIter& CellIter::operator--() {

    --cur_index_;
    while (cur_index_ >= 0 && BaseIter::mesh()->is_deleted(CellHandle(cur_index_)))
        --cur_index_;
    if(cur_index_ < 0) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(CellHandle(cur_index_));
    return *this;
}


CellIter& CellIter::operator++() {

    ++cur_index_;
    while ((unsigned int)cur_index_ < BaseIter::mesh()->cells_.size() && BaseIter::mesh()->is_deleted(CellHandle(cur_index_)))
        ++cur_index_;
    if((unsigned int)cur_index_ >= BaseIter::mesh()->cells_.size()) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(CellHandle(cur_index_));
    return *this;
}

namespace Internal {

////================================================================================================
//// VertexIHalfEdgeIterImpl
////================================================================================================

VertexIHalfEdgeIterImpl::VertexIHalfEdgeIterImpl(const VertexHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps) :
    BaseIter(_mesh, _ref_h, _max_laps),
    voh_iter_(_mesh->voh_iter(_ref_h, _max_laps))
{
    BaseIter::valid(voh_iter_.valid());
    if (BaseIter::valid()) {
        BaseIter::cur_handle(mesh()->opposite_halfedge_handle(*voh_iter_));
    }
}

VertexIHalfEdgeIterImpl& VertexIHalfEdgeIterImpl::operator--() {
    --voh_iter_;
    BaseIter::lap(voh_iter_.lap());
    BaseIter::valid(voh_iter_.valid());
    BaseIter::cur_handle(mesh()->opposite_halfedge_handle(*voh_iter_));
    return *this;
}

VertexIHalfEdgeIterImpl& VertexIHalfEdgeIterImpl::operator++() {
    ++voh_iter_;
    BaseIter::lap(voh_iter_.lap());
    BaseIter::valid(voh_iter_.valid());
    BaseIter::cur_handle(mesh()->opposite_halfedge_handle(*voh_iter_));
    return *this;
}

////================================================================================================
//// VertexEdgeIterImpl
////================================================================================================

VertexEdgeIterImpl::VertexEdgeIterImpl(const VertexHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps) :
    BaseIter(_mesh, _ref_h, _max_laps),
    voh_iter_(_mesh->voh_iter(_ref_h, _max_laps))
{
    BaseIter::valid(voh_iter_.valid());
    if (BaseIter::valid()) {
        BaseIter::cur_handle(mesh()->edge_handle(*voh_iter_));
    }
}

VertexEdgeIterImpl& VertexEdgeIterImpl::operator--() {
    --voh_iter_;
    BaseIter::lap(voh_iter_.lap());
    BaseIter::valid(voh_iter_.valid());
    BaseIter::cur_handle(mesh()->edge_handle(*voh_iter_));
    return *this;
}

VertexEdgeIterImpl& VertexEdgeIterImpl::operator++() {
    ++voh_iter_;
    BaseIter::lap(voh_iter_.lap());
    BaseIter::valid(voh_iter_.valid());
    BaseIter::cur_handle(mesh()->edge_handle(*voh_iter_));
    return *this;
}

////================================================================================================
//// VertexHalfFaceIterImpl
////================================================================================================

VertexHalfFaceIterImpl::VertexHalfFaceIterImpl(const VertexHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps) :
    BaseIter(_mesh, _ref_h, _max_laps),
    cur_index_(0)
{
    for (VertexEdgeIter ve_iter = _mesh->ve_iter(_ref_h); ve_iter.valid(); ++ve_iter) {
        for (EdgeHalfFaceIter ehf_iter = _mesh->ehf_iter(*ve_iter); ehf_iter.valid(); ++ehf_iter) {
            halffaces_.push_back(*ehf_iter);
        }
    }

    // Remove all duplicate entries
    std::sort(halffaces_.begin(), halffaces_.end());
    halffaces_.resize(std::unique(halffaces_.begin(), halffaces_.end()) - halffaces_.begin());

    BaseIter::valid(halffaces_.size() > 0);
    if (BaseIter::valid()) {
        BaseIter::cur_handle(halffaces_[cur_index_]);
    }
}

VertexHalfFaceIterImpl& VertexHalfFaceIterImpl::operator--() {
    if (cur_index_ == 0) {
        cur_index_ = halffaces_.size() - 1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else {
        --cur_index_;
    }
    BaseIter::cur_handle(halffaces_[cur_index_]);
    return *this;
}

VertexHalfFaceIterImpl& VertexHalfFaceIterImpl::operator++() {
    ++cur_index_;
    if (cur_index_ >= halffaces_.size()) {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }
    BaseIter::cur_handle(halffaces_[cur_index_]);
    return *this;
}


////================================================================================================
//// HalfEdgeFaceIterImpl
////================================================================================================

HalfEdgeFaceIterImpl::HalfEdgeFaceIterImpl(const HalfEdgeHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps) :
    BaseIter(_mesh, _ref_h, _max_laps),
    cur_index_(0)
{
    for (HalfEdgeHalfFaceIter hehf_iter = _mesh->hehf_iter(_ref_h); hehf_iter.valid(); ++hehf_iter) {
        faces_.push_back(_mesh->face_handle(*hehf_iter));
    }

    // Remove all duplicate entries
    std::sort(faces_.begin(), faces_.end());
    faces_.resize(std::unique(faces_.begin(), faces_.end()) - faces_.begin());

    BaseIter::valid(faces_.size() > 0);
    if (BaseIter::valid()) {
        BaseIter::cur_handle(faces_[cur_index_]);
    }
}

HalfEdgeFaceIterImpl& HalfEdgeFaceIterImpl::operator--() {
    if (cur_index_ == 0) {
        cur_index_ = faces_.size() - 1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else {
        --cur_index_;
    }
    BaseIter::cur_handle(faces_[cur_index_]);
    return *this;
}

HalfEdgeFaceIterImpl& HalfEdgeFaceIterImpl::operator++() {
    ++cur_index_;
    if (cur_index_ >= faces_.size()) {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }
    BaseIter::cur_handle(faces_[cur_index_]);
    return *this;
}


////================================================================================================
//// EdgeHalfFaceIterImpl
////================================================================================================

EdgeHalfFaceIterImpl::EdgeHalfFaceIterImpl(const EdgeHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps) :
    BaseIter(_mesh, _ref_h, _max_laps),
    cur_index_(0)
{
    const HalfEdgeHandle he = _mesh->halfedge_handle(_ref_h, 0);
    for (HalfEdgeHalfFaceIter hehf_iter = _mesh->hehf_iter(he); hehf_iter.valid(); ++hehf_iter) {
        halffaces_.push_back(*hehf_iter);
        halffaces_.push_back(_mesh->opposite_halfface_handle(*hehf_iter));
    }

    BaseIter::valid(halffaces_.size() > 0);
    if (BaseIter::valid()) {
        BaseIter::cur_handle(halffaces_[cur_index_]);
    }
}

EdgeHalfFaceIterImpl& EdgeHalfFaceIterImpl::operator--() {
    if (cur_index_ == 0) {
        cur_index_ = halffaces_.size() - 1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else {
        --cur_index_;
    }
    BaseIter::cur_handle(halffaces_[cur_index_]);
    return *this;
}

EdgeHalfFaceIterImpl& EdgeHalfFaceIterImpl::operator++() {
    ++cur_index_;
    if (cur_index_ >= halffaces_.size()) {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }
    BaseIter::cur_handle(halffaces_[cur_index_]);
    return *this;
}


////================================================================================================
//// EdgeFaceIterImpl
////================================================================================================

EdgeFaceIterImpl::EdgeFaceIterImpl(const EdgeHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps) :
    HalfEdgeFaceIterImpl(_mesh->halfedge_handle(_ref_h, 0), _mesh, _max_laps) {}


////================================================================================================
//// EdgeCellIterImpl
////================================================================================================

EdgeCellIterImpl::EdgeCellIterImpl(const EdgeHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps) :
    HalfEdgeCellIter(_mesh->halfedge_handle(_ref_h, 0), _mesh, _max_laps) {}


////================================================================================================
//// HalfFaceHalfEdgeIterImpl
////================================================================================================

HalfFaceHalfEdgeIterImpl::HalfFaceHalfEdgeIterImpl(const HalfFaceHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps) :
    BaseIter(_mesh, _ref_h, _max_laps),
    cur_index_(0)
{
    BaseIter::valid(_ref_h.is_valid() && _mesh->halfface(_ref_h).halfedges().size() > 0);
    if (BaseIter::valid()) {
        BaseIter::cur_handle(_mesh->halfface(_ref_h).halfedges()[cur_index_]);
    }
}

HalfFaceHalfEdgeIterImpl& HalfFaceHalfEdgeIterImpl::operator--() {
    const std::vector<HalfEdgeHandle> halfedges =
        mesh()->halfface(ref_handle()).halfedges();
    if (cur_index_ == 0) {
        cur_index_ = halfedges.size() - 1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else {
        --cur_index_;
    }
    BaseIter::cur_handle(halfedges[cur_index_]);
    return *this;
}

HalfFaceHalfEdgeIterImpl& HalfFaceHalfEdgeIterImpl::operator++() {
    const std::vector<HalfEdgeHandle> halfedges =
        mesh()->halfface(ref_handle()).halfedges();
    ++cur_index_;
    if (cur_index_ >= halfedges.size()) {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }
    BaseIter::cur_handle(halfedges[cur_index_]);
    return *this;
}


////================================================================================================
//// HalfFaceEdgeIterImpl
////================================================================================================

HalfFaceEdgeIterImpl::HalfFaceEdgeIterImpl(const HalfFaceHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps) :
    BaseIter(_mesh, _ref_h, _max_laps),
    cur_index_(0)
{
    BaseIter::valid(_ref_h.is_valid() && _mesh->halfface(_ref_h).halfedges().size() > 0);
    if (BaseIter::valid()) {
        HalfEdgeHandle he = _mesh->halfface(_ref_h).halfedges()[cur_index_];
        BaseIter::cur_handle(_mesh->edge_handle(he));
    }
}

HalfFaceEdgeIterImpl& HalfFaceEdgeIterImpl::operator--() {
    const std::vector<HalfEdgeHandle> halfedges =
        mesh()->halfface(ref_handle()).halfedges();
    if (cur_index_ == 0) {
        cur_index_ = halfedges.size() - 1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else {
        --cur_index_;
    }
    BaseIter::cur_handle(mesh()->edge_handle(halfedges[cur_index_]));
    return *this;
}

HalfFaceEdgeIterImpl& HalfFaceEdgeIterImpl::operator++() {
    const std::vector<HalfEdgeHandle> halfedges =
        mesh()->halfface(ref_handle()).halfedges();
    ++cur_index_;
    if (cur_index_ >= halfedges.size()) {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }
    BaseIter::cur_handle(mesh()->edge_handle(halfedges[cur_index_]));
    return *this;
}

////================================================================================================
//// FaceVertexIterImpl
////================================================================================================

FaceVertexIterImpl::FaceVertexIterImpl(const FaceHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps) :
    HalfFaceVertexIter(_mesh->halfface_handle(_ref_h, 0), _mesh, _max_laps) {}

////================================================================================================
//// FaceHalfEdgeIterImpl
////================================================================================================

FaceHalfEdgeIterImpl::FaceHalfEdgeIterImpl(const FaceHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps) :
    HalfFaceHalfEdgeIterImpl(_mesh->halfface_handle(_ref_h, 0), _mesh, _max_laps) {}


////================================================================================================
//// FaceEdgeIterImpl
////================================================================================================

FaceEdgeIterImpl::FaceEdgeIterImpl(const FaceHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps) :
    HalfFaceEdgeIterImpl(_mesh->halfface_handle(_ref_h, 0), _mesh, _max_laps) {}


////================================================================================================
//// CellHalfEdgeIterImpl
////================================================================================================

CellHalfEdgeIterImpl::CellHalfEdgeIterImpl(const CellHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps) :
    BaseIter(_mesh, _ref_h, _max_laps),
    cur_index_(0)
{
    for (CellHalfFaceIter chf_iter =_mesh->chf_iter(_ref_h); chf_iter.valid(); ++chf_iter) {
        for (HalfFaceHalfEdgeIter hfhe_iter =_mesh->hfhe_iter(*chf_iter); hfhe_iter.valid(); ++hfhe_iter) {
            halfedges_.push_back(*hfhe_iter);
        }
    }
    BaseIter::valid(halfedges_.size() > 0);
    if (BaseIter::valid()) {
        BaseIter::cur_handle(halfedges_[cur_index_]);
    }
}

CellHalfEdgeIterImpl& CellHalfEdgeIterImpl::operator--() {
    if (cur_index_ == 0) {
        cur_index_ = halfedges_.size() - 1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else {
        --cur_index_;
    }
    BaseIter::cur_handle(halfedges_[cur_index_]);
    return *this;
}

CellHalfEdgeIterImpl& CellHalfEdgeIterImpl::operator++() {
    ++cur_index_;
    if (cur_index_ >= halfedges_.size()) {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }
    BaseIter::cur_handle(halfedges_[cur_index_]);
    return *this;
}


////================================================================================================
//// CellEdgeIterImpl
////================================================================================================

CellEdgeIterImpl::CellEdgeIterImpl(const CellHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps) :
    BaseIter(_mesh, _ref_h, _max_laps),
    cur_index_(0)
{
    for (CellHalfEdgeIter che_iter = _mesh->che_iter(_ref_h); che_iter.valid(); ++che_iter) {
        edges_.push_back(_mesh->edge_handle(*che_iter));
    }

    // Remove all duplicate entries
    std::sort(edges_.begin(), edges_.end());
    edges_.resize(std::unique(edges_.begin(), edges_.end()) - edges_.begin());

    BaseIter::valid(edges_.size() > 0);
    if (BaseIter::valid()) {
        BaseIter::cur_handle(edges_[cur_index_]);
    }
}

CellEdgeIterImpl& CellEdgeIterImpl::operator--() {
    if (cur_index_ == 0) {
        cur_index_ = edges_.size() - 1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else {
        --cur_index_;
    }
    BaseIter::cur_handle(edges_[cur_index_]);
    return *this;
}

CellEdgeIterImpl& CellEdgeIterImpl::operator++() {
    ++cur_index_;
    if (cur_index_ >= edges_.size()) {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }
    BaseIter::cur_handle(edges_[cur_index_]);
    return *this;
}


////================================================================================================
//// CellHalfFaceIterImpl
////================================================================================================

CellHalfFaceIterImpl::CellHalfFaceIterImpl(const CellHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps) :
    BaseIter(_mesh, _ref_h, _max_laps),
    hf_iter_(BaseIter::mesh()->cell(_ref_h).halffaces().begin())
{
    BaseIter::valid(hf_iter_ != BaseIter::mesh()->cell(_ref_h).halffaces().end());
    if (BaseIter::valid()) {
        BaseIter::cur_handle(*hf_iter_);
    }
}

CellHalfFaceIterImpl& CellHalfFaceIterImpl::operator--() {
    const std::vector<HalfFaceHandle>& halffaces =
        BaseIter::mesh()->cell(ref_handle_).halffaces();
    if (hf_iter_ == halffaces.begin()) {
        hf_iter_ = halffaces.end();
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else {
        --hf_iter_;
    }
    BaseIter::cur_handle(*hf_iter_);
    return *this;
}

CellHalfFaceIterImpl& CellHalfFaceIterImpl::operator++() {
    ++hf_iter_;
    const std::vector<HalfFaceHandle>& halffaces =
        BaseIter::mesh()->cell(ref_handle_).halffaces();
    if (hf_iter_ == halffaces.end()) {
        hf_iter_ = halffaces.begin();
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }
    BaseIter::cur_handle(*hf_iter_);
    return *this;
}


////================================================================================================
//// CellFaceIterImpl
////================================================================================================

CellFaceIterImpl::CellFaceIterImpl(const CellHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps) :
    BaseIter(_mesh, _ref_h, _max_laps),
    hf_iter_(BaseIter::mesh()->cell(_ref_h).halffaces().begin())
{
    BaseIter::valid(hf_iter_ != BaseIter::mesh()->cell(_ref_h).halffaces().end());
    if (BaseIter::valid()) {
        BaseIter::cur_handle(BaseIter::mesh()->face_handle(*hf_iter_));
    }
}

CellFaceIterImpl& CellFaceIterImpl::operator--() {
    const std::vector<HalfFaceHandle>& halffaces =
        BaseIter::mesh()->cell(ref_handle_).halffaces();
    if (hf_iter_ == halffaces.begin()) {
        hf_iter_ = halffaces.end();
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else {
        --hf_iter_;
    }
    BaseIter::cur_handle(BaseIter::mesh()->face_handle(*hf_iter_));
    return *this;
}

CellFaceIterImpl& CellFaceIterImpl::operator++() {
    ++hf_iter_;
    const std::vector<HalfFaceHandle>& halffaces =
        BaseIter::mesh()->cell(ref_handle_).halffaces();
    if (hf_iter_ == halffaces.end()) {
        hf_iter_ = halffaces.begin();
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }
    BaseIter::cur_handle(BaseIter::mesh()->face_handle(*hf_iter_));
    return *this;
}

}


////================================================================================================
//// BoundaryItemIter
////================================================================================================

template <>
size_t BoundaryItemIter<VertexIter, VertexHandle>::n_items() const {
    return BaseIter::mesh()->n_vertices();
}

template <>
size_t BoundaryItemIter<HalfEdgeIter, HalfEdgeHandle>::n_items() const {
    return BaseIter::mesh()->n_halfedges();
}

template <>
size_t BoundaryItemIter<EdgeIter, EdgeHandle>::n_items() const {
    return BaseIter::mesh()->n_edges();
}

template <>
size_t BoundaryItemIter<HalfFaceIter, HalfFaceHandle>::n_items() const {
    return BaseIter::mesh()->n_halffaces();
}

template <>
size_t BoundaryItemIter<FaceIter, FaceHandle>::n_items() const {
    return BaseIter::mesh()->n_faces();
}

template <>
size_t BoundaryItemIter<CellIter, CellHandle>::n_items() const {
    return BaseIter::mesh()->n_cells();
}

template <>
bool BoundaryItemIter<VertexIter, VertexHandle>::has_incidences() const {
    return BaseIter::mesh()->has_full_bottom_up_incidences();
}

template <>
bool BoundaryItemIter<HalfEdgeIter, HalfEdgeHandle>::has_incidences() const {
    const TopologyKernel *mesh = BaseIter::mesh();
    return mesh->has_edge_bottom_up_incidences() && mesh->has_face_bottom_up_incidences();
}

template <>
bool BoundaryItemIter<EdgeIter, EdgeHandle>::has_incidences() const {
    const TopologyKernel *mesh = BaseIter::mesh();
    return mesh->has_edge_bottom_up_incidences() && mesh->has_face_bottom_up_incidences();
}

template <>
bool BoundaryItemIter<HalfFaceIter, HalfFaceHandle>::has_incidences() const {
    return BaseIter::mesh()->has_face_bottom_up_incidences();
}

template <>
bool BoundaryItemIter<FaceIter, FaceHandle>::has_incidences() const {
    return BaseIter::mesh()->has_face_bottom_up_incidences();
}

template <>
bool BoundaryItemIter<CellIter, CellHandle>::has_incidences() const {
    return true;
}

} // Namespace OpenVolumeMesh
