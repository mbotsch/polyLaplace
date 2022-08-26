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

#pragma once

#ifndef NDEBUG
#include <iostream>
#endif
#include <string>
#include <vector>
#include <type_traits>

#include "../System/Compiler.hh"
#include "OpenVolumeMesh/Config/Export.hh"
#include "OpenVolumeMeshProperty.hh"
#include "PropertyHandles.hh"
#include "TypeName.hh"
#include "ForwardDeclarations.hh"

#if OVM_CXX_17
#include <optional>
#endif

namespace OpenVolumeMesh {

// Forward declarations
class BaseProperty;

class OVM_EXPORT ResourceManager {
public:
    ResourceManager() = default;
    ResourceManager(const ResourceManager &other);
    ResourceManager(ResourceManager &&other);
    ResourceManager& operator=(const ResourceManager &other);
    ResourceManager& operator=(ResourceManager &&other);
    virtual ~ResourceManager();

    template <class PropT, class HandleT> friend class PropertyPtr;

    /// Change size of stored vertex properties
    void resize_vprops(size_t _nv);

    /// Change size of stored edge properties
    void resize_eprops(size_t _ne);

    /// Change size of stored face properties
    void resize_fprops(size_t _nf);

    /// Change size of stored cell properties
    void resize_cprops(size_t _nc);

protected:

    void vertex_deleted(const VertexHandle& _h);

    void edge_deleted(const EdgeHandle& _h);

    void face_deleted(const FaceHandle& _h);

    void cell_deleted(const CellHandle& _h);

    void swap_cell_properties(CellHandle _h1, CellHandle _h2);

    void swap_face_properties(FaceHandle _h1, FaceHandle _h2);

    void swap_halfface_properties(HalfFaceHandle _h1, HalfFaceHandle _h2);

    void swap_edge_properties(EdgeHandle _h1, EdgeHandle _h2);

    void swap_halfedge_properties(HalfEdgeHandle _h1, HalfEdgeHandle _h2);

    void swap_vertex_properties(VertexHandle _h1, VertexHandle _h2);

    template <typename PropIterator, typename Handle>
    void swap_property_elements(PropIterator _begin, PropIterator _end, Handle _h1, Handle _h2)
    {
        PropIterator p_iter =  _begin;
        for (; p_iter != _end; ++p_iter)
            (*p_iter)->swap_elements(_h1.uidx(), _h2.uidx());
    }


public:

    void clear_vertex_props() { clearVec(vertex_props_); }

    void clear_edge_props() { clearVec(edge_props_); }

    void clear_halfedge_props() { clearVec(halfedge_props_); }

    void clear_face_props() { clearVec(face_props_); }

    void clear_halfface_props() { clearVec(halfface_props_); }

    void clear_cell_props() { clearVec(cell_props_); }

    void clear_mesh_props() { clearVec(mesh_props_); }

    /// Get number of vertices in mesh
    virtual size_t n_vertices() const = 0;
    /// Get number of edges in mesh
    virtual size_t n_edges() const = 0;
    /// Get number of halfedges in mesh
    virtual size_t n_halfedges() const = 0;
    /// Get number of faces in mesh
    virtual size_t n_faces() const = 0;
    /// Get number of halffaces in mesh
    virtual size_t n_halffaces() const = 0;
    /// Get number of cells in mesh
    virtual size_t n_cells() const = 0;


    /** Get or create property: if the property does not exist yet, create it.
     */
    template<typename T, typename EntityTag>
    PropertyTT<T, EntityTag> request_property(const std::string& _name = std::string(), const T _def = T());

#if OVM_CXX_17
    /** Create new property: if the property already exists, return no value.
     */
    template<typename T, typename EntityTag>
    std::optional<PropertyTT<T, EntityTag>> create_property(const std::string& _name = std::string(), const T _def = T());

    /** Get existing property: if the property does not exist, return no value.
     */
    template<typename T, typename EntityTag>
    std::optional<PropertyTT<T, EntityTag>> get_property(const std::string& _name = std::string());
#endif

    template<class T> VertexPropertyT<T> request_vertex_property(const std::string& _name = std::string(), const T _def = T());

    template<class T> EdgePropertyT<T> request_edge_property(const std::string& _name = std::string(), const T _def = T());

    template<class T> HalfEdgePropertyT<T> request_halfedge_property(const std::string& _name = std::string(), const T _def = T());

    template<class T> FacePropertyT<T> request_face_property(const std::string& _name = std::string(), const T _def = T());

    template<class T> HalfFacePropertyT<T> request_halfface_property(const std::string& _name = std::string(), const T _def = T());

    template<class T> CellPropertyT<T> request_cell_property(const std::string& _name = std::string(), const T _def = T());

    template<class T> MeshPropertyT<T> request_mesh_property(const std::string& _name = std::string(), const T _def = T());


private:

    void release_property(VertexPropHandle _handle);

    void release_property(EdgePropHandle _handle);

    void release_property(HalfEdgePropHandle _handle);

    void release_property(FacePropHandle _handle);

    void release_property(HalfFacePropHandle _handle);

    void release_property(CellPropHandle _handle);

    void release_property(MeshPropHandle _handle);

public:

    size_t n_vertex_props() const { return vertex_props_.size(); }

    size_t n_edge_props() const { return edge_props_.size(); }

    size_t n_halfedge_props() const { return halfedge_props_.size(); }

    size_t n_face_props() const { return face_props_.size(); }

    size_t n_halfface_props() const { return halfface_props_.size(); }

    size_t n_cell_props() const { return cell_props_.size(); }

    size_t n_mesh_props() const { return mesh_props_.size(); }

    template<typename T, class EntityTag>
    void set_persistent(PropertyTT<T, EntityTag>& _prop, bool _flag = true);

    typedef std::vector<BaseProperty*> Properties;

    Properties::const_iterator vertex_props_begin() const { return vertex_props_.begin(); }

    Properties::const_iterator vertex_props_end() const { return vertex_props_.end(); }

    Properties::const_iterator edge_props_begin() const { return edge_props_.begin(); }

    Properties::const_iterator edge_props_end() const { return edge_props_.end(); }

    Properties::const_iterator halfedge_props_begin() const { return halfedge_props_.begin(); }

    Properties::const_iterator halfedge_props_end() const { return halfedge_props_.end(); }

    Properties::const_iterator face_props_begin() const { return face_props_.begin(); }

    Properties::const_iterator face_props_end() const { return face_props_.end(); }

    Properties::const_iterator halfface_props_begin() const { return halfface_props_.begin(); }

    Properties::const_iterator halfface_props_end() const { return halfface_props_.end(); }

    Properties::const_iterator cell_props_begin() const { return cell_props_.begin(); }

    Properties::const_iterator cell_props_end() const { return cell_props_.end(); }

    Properties::const_iterator mesh_props_begin() const { return mesh_props_.begin(); }

    Properties::const_iterator mesh_props_end() const { return mesh_props_.end(); }

private:

    template <class FullPropT, class PropIterT>
    bool property_exists(const PropIterT& _begin, const PropIterT& _end, const std::string& _name) const
    {
        auto type_name = get_type_name<typename FullPropT::value_type>();

        if(_name.empty()) {
#ifndef NDEBUG
            std::cerr << "property_exists(): Checking for the existence of anonymous properties is" << std::endl;
            std::cerr << "ambiguous!" << std::endl;
#endif
            return false;
        }

        PropIterT it = _begin;
        for(; it != _end; ++it)
        {
            if((*it)->name() == _name
                && (*it)->internal_type_name() == type_name)
            {
                return true;
            }
        }
        return false;
    }

public:

    template <class PropT>
    bool vertex_property_exists(const std::string& _name) const {
        return property_exists<VertexPropertyT<PropT> >(vertex_props_begin(), vertex_props_end(), _name);
    }

    template <class PropT>
    bool edge_property_exists(const std::string& _name) const {
        return property_exists<EdgePropertyT<PropT> >(edge_props_begin(), edge_props_end(), _name);
    }

    template <class PropT>
    bool halfedge_property_exists(const std::string& _name) const {
        return property_exists<HalfEdgePropertyT<PropT> >(halfedge_props_begin(), halfedge_props_end(), _name);
    }

    template <class PropT>
    bool face_property_exists(const std::string& _name) const {
        return property_exists<FacePropertyT<PropT> >(face_props_begin(), face_props_end(), _name);
    }

    template <class PropT>
    bool halfface_property_exists(const std::string& _name) const {
        return property_exists<HalfFacePropertyT<PropT> >(halfface_props_begin(), halfface_props_end(), _name);
    }

    template <class PropT>
    bool cell_property_exists(const std::string& _name) const {
        return property_exists<CellPropertyT<PropT> >(cell_props_begin(), cell_props_end(), _name);
    }

    template <class PropT>
    bool mesh_property_exists(const std::string& _name) const {
        return property_exists<MeshPropertyT<PropT> >(mesh_props_begin(), mesh_props_end(), _name);
    }

protected:

    void delete_multiple_vertex_props(const std::vector<bool>& _tags);

    void delete_multiple_edge_props(const std::vector<bool>& _tags);

    void delete_multiple_face_props(const std::vector<bool>& _tags);

    void delete_multiple_cell_props(const std::vector<bool>& _tags);

private:

    template<class StdVecT>
    void resize_props(StdVecT& _vec, size_t _n);

    template<class StdVecT>
    void entity_deleted(StdVecT& _vec, const OpenVolumeMeshHandle& _h);

    template<class StdVecT>
    void remove_property(StdVecT& _vec, size_t _idx);

    template<typename T, typename EntityTag>
    PropertyTT<T, EntityTag> *internal_find_property(const std::string& _name);

    template<typename T, typename EntityTag>
    PropertyTT<T, EntityTag> internal_create_property(const std::string& _name, const T _def = T());

    template<class StdVecT>
    void clearVec(StdVecT& _vec);

    template<class StdVecT>
    void updatePropHandles(StdVecT& _vec);

    template<bool Move>
    void assignProperties(typename std::conditional<Move, Properties&, const Properties&>::type src,
                          Properties &dest);
    template<bool Move>
    void assignAllPropertiesFrom(typename std::conditional<Move, ResourceManager*, const ResourceManager*>::type src);

    Properties vertex_props_;

    Properties edge_props_;

    Properties halfedge_props_;

    Properties face_props_;

    Properties halfface_props_;

    Properties cell_props_;

    Properties mesh_props_;

    template<typename Entity>
    Properties &entity_props();

    template<typename Entity>
    size_t n();
};

}

#if defined(INCLUDE_TEMPLATES) && !defined(RESOURCEMANAGERT_CC)
#include "ResourceManagerT_impl.hh"
#endif

