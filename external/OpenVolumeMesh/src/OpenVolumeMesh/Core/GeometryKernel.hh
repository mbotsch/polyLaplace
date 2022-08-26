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

#ifndef GEOMETRYKERNEL_HH_
#define GEOMETRYKERNEL_HH_

#include <cassert>
#include <iostream>
#include <type_traits>

#include "../Geometry/VectorT.hh"
#include "TopologyKernel.hh"

namespace OpenVolumeMesh {

template <class VecT, class TopologyKernelT = TopologyKernel>
class GeometryKernel : public TopologyKernelT {
public:

    typedef VecT PointT;
    typedef TopologyKernelT KernelT;

    /// Constructor
    GeometryKernel() = default;

    /// Destructor
    ~GeometryKernel() override = default;

    template<class OtherTopoKernel>
    void assign(const GeometryKernel<VecT, OtherTopoKernel> *other) {
        TopologyKernelT::assign(other);
        other->clone_vertices(vertices_);
    }

    /// Override of empty add_vertex function
    VertexHandle add_vertex() override { return add_vertex(VecT()); }

    /// Add a geometric point to the mesh
    VertexHandle add_vertex(const VecT& _p) {

        // Store vertex in list
        vertices_.push_back(_p);

        // Get handle of recently created vertex
        return KernelT::add_vertex();
    }

    /// Set the coordinates of point _vh
    void set_vertex(const VertexHandle& _vh, const VecT& _p) {

        assert(_vh.idx() < (int)vertices_.size());

        vertices_[_vh.uidx()] = _p;
    }

    /// Get point _vh's coordinates
    const VecT& vertex(const VertexHandle& _vh) const {
        return vertices_[_vh.uidx()];
    }

    VertexIter delete_vertex(const VertexHandle& _h) override {
        assert(_h.idx() < (int)TopologyKernelT::n_vertices());

        VertexIter nV = TopologyKernelT::delete_vertex(_h);

        if (TopologyKernelT::deferred_deletion_enabled())
        {

        }
        else
            vertices_.erase(vertices_.begin() + _h.idx());

        return nV;
    }

    void collect_garbage() override
    {
        if (!TopologyKernelT::needs_garbage_collection())
            return;

        if (TopologyKernelT::fast_deletion_enabled()) {
            TopologyKernelT::collect_garbage();
            vertices_.resize(TopologyKernel::n_vertices());
        } else {
            for (int i = (int)vertices_.size(); i > 0; --i)
                if (TopologyKernelT::is_deleted(VertexHandle(i-1)))
                {
                    vertices_.erase(vertices_.begin() + (i-1));
                }
            TopologyKernelT::collect_garbage();
        }

    }

    void swap_vertex_indices(VertexHandle _h1, VertexHandle _h2) override
    {
        assert(_h1.idx() >= 0 && _h1.idx() < (int)vertices_.size());
        assert(_h2.idx() >= 0 && _h2.idx() < (int)vertices_.size());

        if (_h1 == _h2)
            return;

        std::swap(vertices_[_h1.uidx()], vertices_[_h2.uidx()]);

        TopologyKernelT::swap_vertex_indices(_h1, _h2);
    }

protected:

    void delete_multiple_vertices(const std::vector<bool>& _tag) override{

        assert(_tag.size() == TopologyKernelT::n_vertices());

        std::vector<VecT> newVertices;

        typename std::vector<VecT>::const_iterator v_it = vertices_.begin();

        for(std::vector<bool>::const_iterator t_it = _tag.begin(),
                t_end = _tag.end(); t_it != t_end; ++t_it, ++v_it) {

            if(!(*t_it)) {
                // Not marked as deleted

                newVertices.push_back(*v_it);
            }
        }

        // Swap vertices
        vertices_.swap(newVertices);

        TopologyKernelT::delete_multiple_vertices(_tag);
    }

public:

    void clear(bool _clearProps = true) override {

        vertices_.clear();
        TopologyKernelT::clear(_clearProps);
    }

    typename PointT::value_type length(const HalfEdgeHandle& _heh) const {
        return vector(_heh).length();
    }

    typename PointT::value_type length(const EdgeHandle& _eh) const {
        return vector(_eh).length();
    }

    PointT vector(const HalfEdgeHandle& _heh) const {

        const typename TopologyKernelT::Edge& e = TopologyKernelT::halfedge(_heh);
        return (vertex(e.to_vertex()) - vertex(e.from_vertex()));
    }

    PointT vector(const EdgeHandle& _eh) const {

        const typename TopologyKernelT::Edge& e = TopologyKernelT::edge(_eh);
        return (vertex(e.to_vertex()) - vertex(e.from_vertex()));
    }

    PointT barycenter(const EdgeHandle& _eh) const {
        return PointT(0.5 * vertex(TopologyKernelT::edge(_eh).from_vertex()) +
                      0.5 * vertex(TopologyKernelT::edge(_eh).to_vertex()));
    }

    PointT barycenter(const FaceHandle& _fh) const {
        PointT p(typename PointT::value_type(0));
        typename PointT::value_type valence = 0;
        HalfFaceVertexIter hfv_it =
                TopologyKernelT::hfv_iter(TopologyKernelT::halfface_handle(_fh, 0));
        for(; hfv_it.valid(); ++hfv_it, valence += 1) {
            p += vertex(*hfv_it);
        }
        p /= valence;
        return p;
    }

    PointT barycenter(const CellHandle& _ch) const {
        PointT p(typename PointT::value_type(0));
        typename PointT::value_type valence = 0;
        CellVertexIter cv_it = TopologyKernelT::cv_iter(_ch);
        for(; cv_it.valid(); ++cv_it, valence += 1) {
            p += vertex(*cv_it);
        }
        p /= valence;
        return p;
    }

    /// Compute halfface normal assuming planarity (just uses first 2 edges)
    /// Note: NormalAttrib provides fast access to precomputed normals.
    PointT normal(const HalfFaceHandle& _hfh) const
    {
        if(TopologyKernelT::halfface(_hfh).halfedges().size() < 3) {
            std::cerr << "Warning: Degenerate face: "
                      << TopologyKernelT::face_handle(_hfh) << std::endl;
            return PointT {0.0, 0.0, 0.0};
        }

        const std::vector<HalfEdgeHandle>& halfedges = TopologyKernelT::halfface(_hfh).halfedges();
        std::vector<HalfEdgeHandle>::const_iterator he_it = halfedges.begin();

        const PointT &p1 = vertex(TopologyKernelT::halfedge(*he_it).from_vertex());
        const PointT &p2 = vertex(TopologyKernelT::halfedge(*he_it).to_vertex());
        ++he_it;
        const PointT &p3 = vertex(TopologyKernelT::halfedge(*he_it).to_vertex());

        const PointT n = (p2 - p1).cross(p3 - p2);
        return n.normalized();
    }

    void clone_vertices(std::vector<VecT>& _copy) const {
        _copy.clear();
        _copy.reserve(vertices_.size());
        std::copy(vertices_.begin(), vertices_.end(), std::back_inserter(_copy));
    }

    void swap_vertices(std::vector<VecT>& _copy) {
        if(_copy.size() != vertices_.size()) {
            std::cerr << "Vertex vectors differ in size! The size of the copy " <<
            		"is artificially set to the correct one. Some values may not be correctly initialized." << std::endl;
            _copy.resize(vertices_.size());
        }
        std::swap(vertices_, _copy);
    }

private:

    std::vector<VecT> vertices_;
};

} // Namespace OpenVolumeMesh

#endif /* GEOMETRYKERNEL_HH_ */
