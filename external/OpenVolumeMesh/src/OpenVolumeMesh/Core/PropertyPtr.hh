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

#include <string>

#include "PropertyHandles.hh"
#include "BaseProperty.hh"
#include "OpenVolumeMeshHandle.hh"
#include "../System/MemoryInclude.hh"
#include "../System/Deprecation.hh"

namespace OpenVolumeMesh {

class ResourceManager;

/**
 * \class PropertyPtr
 *
 * A smart-pointer-like class that counts the encapsulated
 * object's references and automatically deletes the memory
 * as soon as the object is not used anymore.
 */

template <class PropT, typename Entity>
class PropertyPtr : protected ptr::shared_ptr<PropT>, public BaseProperty {
public:

    friend class ResourceManager;

    typedef typename PropT::value_type                  value_type;
    typedef typename PropT::vector_type::const_iterator const_iterator;
    typedef typename PropT::vector_type::iterator       iterator;
    typedef typename PropT::reference                   reference;
    typedef typename PropT::const_reference             const_reference;

    using EntityHandleT = HandleT<Entity>;

    /// Constructor
    PropertyPtr() : BaseProperty(nullptr) {}
    /// Constructor
    PropertyPtr(PropT* _ptr, ResourceManager& _resMan, PropHandleT<Entity> _handle);

    /// Destructor
    ~PropertyPtr() override;

    PropertyPtr(const PropertyPtr<PropT, Entity>&) = default;
    PropertyPtr(PropertyPtr<PropT, Entity>&&) = default;
    PropertyPtr<PropT, Entity>& operator=(const PropertyPtr<PropT, Entity>&) = default;
    PropertyPtr<PropT, Entity>& operator=(PropertyPtr<PropT, Entity>&&) = default;

    using ptr::shared_ptr<PropT>::operator*;
    using ptr::shared_ptr<PropT>::operator->;
    using ptr::shared_ptr<PropT>::operator bool;

    const std::string& name() const override;

    void delete_element(size_t _idx) override;

    void swap_elements(size_t _idx0, size_t _idx1) override;

    void copy(size_t _src_idx, size_t _dst_idx) override;

    const_iterator begin() const { return ptr::shared_ptr<PropT>::get()->begin(); }
    iterator begin() { return ptr::shared_ptr<PropT>::get()->begin(); }
    size_t size() const override { return ptr::shared_ptr<PropT>::get()->size(); }

    const_iterator end() const { return ptr::shared_ptr<PropT>::get()->end(); }
    iterator end() { return ptr::shared_ptr<PropT>::get()->end(); }

#if OVM_ENABLE_DEPRECATED_APIS
    OVM_DEPRECATED("use handles to index properties")
    reference operator[](size_t _idx) { return (*ptr::shared_ptr<PropT>::get())[_idx]; }
    OVM_DEPRECATED("use handles to index properties")
    const_reference operator[](size_t _idx) const { return (*ptr::shared_ptr<PropT>::get())[_idx]; }
#endif

    reference operator[](const EntityHandleT& _h) { return (*ptr::shared_ptr<PropT>::get())[_h.uidx()]; }
    const_reference operator[](const EntityHandleT& _h) const { return (*ptr::shared_ptr<PropT>::get())[_h.uidx()]; }

    void serialize(std::ostream& _ostr) const override { ptr::shared_ptr<PropT>::get()->serialize(_ostr); }
    void deserialize(std::istream& _istr) override { ptr::shared_ptr<PropT>::get()->deserialize(_istr); }

    OpenVolumeMeshHandle handle() const override;

     bool persistent() const override { return ptr::shared_ptr<PropT>::get()->persistent(); }

     bool anonymous() const override { return ptr::shared_ptr<PropT>::get()->name().empty(); }

protected:
    const std::string &internal_type_name() const override { return ptr::shared_ptr<PropT>::get()->internal_type_name(); }

    void assign_values_from(const BaseProperty *other) override;
    void move_values_from(BaseProperty *other) override;

    void delete_multiple_entries(const std::vector<bool>& _tags) override;

    void resize(size_t _size) override;

    void set_handle(const OpenVolumeMeshHandle& _handle) override;
};

} // Namespace OpenVolumeMesh

#if defined(INCLUDE_TEMPLATES) && !defined(PROPERTYPTRT_CC)
#include "PropertyPtrT_impl.hh"
#endif

