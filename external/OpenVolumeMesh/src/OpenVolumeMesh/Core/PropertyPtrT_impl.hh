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

#define PROPERTYPTRT_CC

#include "PropertyPtr.hh"
#include "ResourceManager.hh"
#include "PropertyDefines.hh"

namespace OpenVolumeMesh {

template <class PropT, typename Entity>
PropertyPtr<PropT,Entity>::PropertyPtr(PropT* _ptr, ResourceManager& _resMan, PropHandleT<Entity> _handle) :
    ptr::shared_ptr<PropT>(_ptr), BaseProperty(&_resMan) {
    ptr::shared_ptr<PropT>::get()->set_handle(_handle);
}

template <class PropT, typename Entity>
PropertyPtr<PropT,Entity>::~PropertyPtr() {

    /*
     * If use count is 2 and prop is not set persistent,
     * remove it, since the resource manager is the
     * only one who stores the property.
     */
    if(resMan_ && !persistent() && ptr::shared_ptr<PropT>::use_count() == 2) {
        resMan_->release_property(PropHandleT<Entity>(handle().idx()));
    }
}

template <class PropT, typename Entity>
void PropertyPtr<PropT,Entity>::assign_values_from(const BaseProperty *other) {
    auto _other = static_cast<const PropertyPtr<PropT,Entity>*>(other);
    // FIXME: would be nice to perform a type check here
    ptr::shared_ptr<PropT>::get()->data_vector() = _other->get()->data_vector();
}

template <class PropT, typename Entity>
void PropertyPtr<PropT,Entity>::move_values_from(BaseProperty *other) {
    auto _other = static_cast<PropertyPtr<PropT,Entity>*>(other);
    // FIXME: would be nice to perform a type check here
    ptr::shared_ptr<PropT>::get()->data_vector() = std::move(_other->get()->data_vector());
}

template <class PropT, typename Entity>
void PropertyPtr<PropT,Entity>::resize(size_t _size) {
    ptr::shared_ptr<PropT>::get()->resize(_size);
}

template <class PropT, typename Entity>
const std::string& PropertyPtr<PropT,Entity>::name() const {
    return ptr::shared_ptr<PropT>::get()->name();
}

template <class PropT, typename Entity>
void PropertyPtr<PropT,Entity>::delete_element(size_t _idx) {
    ptr::shared_ptr<PropT>::get()->delete_element(_idx);
}

template <class PropT, typename Entity>
void PropertyPtr<PropT,Entity>::swap_elements(size_t _idx0, size_t _idx1) {
    ptr::shared_ptr<PropT>::get()->swap(_idx0, _idx1);
}

template <class PropT, typename Entity>
void PropertyPtr<PropT,Entity>::copy(size_t _src_idx, size_t _dst_idx) {
    ptr::shared_ptr<PropT>::get()->copy(_src_idx, _dst_idx);
}

template <class PropT, typename Entity>
void PropertyPtr<PropT,Entity>::set_handle(const OpenVolumeMeshHandle& _handle) {
    return ptr::shared_ptr<PropT>::get()->set_handle(_handle);
}

template <class PropT, typename Entity>
OpenVolumeMeshHandle PropertyPtr<PropT,Entity>::handle() const {
    return ptr::shared_ptr<PropT>::get()->handle();
}

template <class PropT, typename Entity>
void PropertyPtr<PropT,Entity>::delete_multiple_entries(const std::vector<bool>& _tags) {
    ptr::shared_ptr<PropT>::get()->delete_multiple_entries(_tags);
}

} // Namespace OpenVolumeMesh
