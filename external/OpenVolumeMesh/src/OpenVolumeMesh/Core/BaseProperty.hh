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

#ifndef BASEPROPERTY_HH_
#define BASEPROPERTY_HH_

#include <string>

#include "OpenVolumeMeshHandle.hh"
#include "OpenVolumeMesh/Config/Export.hh"

namespace OpenVolumeMesh {

class ResourceManager;

class OVM_EXPORT BaseProperty {
public:
    friend class ResourceManager;

    explicit BaseProperty(ResourceManager* _resMan) : resMan_(_resMan) {}

    BaseProperty(const BaseProperty& _other) = default;
    BaseProperty& operator=(const BaseProperty& _cpy) = delete;
    BaseProperty(BaseProperty&& _other) {
        resMan_ = _other.resMan_;
        _other.resMan_ = nullptr;
    }
    BaseProperty& operator=(BaseProperty&& _other) {
        resMan_ = _other.resMan_;
        _other.resMan_ = nullptr;
        return *this;
    }


    virtual ~BaseProperty();

    virtual const std::string& name() const = 0;

    virtual BaseProperty* clone(ResourceManager &_resMan, OpenVolumeMeshHandle _handle) const = 0;

    virtual void delete_element(size_t _idx) = 0;

    virtual void swap_elements(size_t _idx0, size_t _idx1) = 0;

    virtual void copy(size_t _src_idx, size_t _dst_idx) = 0;

    virtual void serialize(std::ostream& _ostr) const = 0;

    virtual void deserialize(std::istream& _istr) = 0;

    virtual OpenVolumeMeshHandle handle() const = 0;

    virtual bool persistent() const = 0;

    virtual bool anonymous() const = 0;

    virtual const std::string entityType() const = 0;

    virtual const std::string typeNameWrapper() const = 0;

    virtual size_t size() const = 0;

protected:

    virtual const std::string &internal_type_name() const = 0;

    /// Copy data from other property. `other` MUST point to an object with the same type as `this`!
    /// Currently no type check is performed.
    virtual void assign_values_from(const BaseProperty *other) = 0;

    /// Move data from other property. `other` MUST point to an object with the same type as `this`!
    /// Currently no type check is performed.
    virtual void move_values_from(BaseProperty *other) = 0;

    virtual void delete_multiple_entries(const std::vector<bool>& _tags) = 0;

    virtual void resize(size_t /*_size*/) = 0;

    virtual void set_handle(const OpenVolumeMeshHandle& /*_handle*/) = 0;

    void setResMan(ResourceManager *resMan) { resMan_ = resMan;}

    ResourceManager* resMan_;
};

} // Namespace OpenVolumeMesh

#endif /* BASEPROPERTY_HH_ */
