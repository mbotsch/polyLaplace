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

#define PROPERTYDEFINEST_CC

#include <istream>
#include <ostream>

#include "PropertyDefines.hh"

namespace OpenVolumeMesh {

template<typename T, typename Entity>
PropertyTT<T,Entity>::PropertyTT(const std::string& _name, const std::string& _internal_type_name, ResourceManager& _resMan, PropertyHandleT _handle, const T &_def) :
        PropertyPtr<OpenVolumeMeshPropertyT<T>, Entity>(new OpenVolumeMeshPropertyT<T>(_name, _internal_type_name, _def), _resMan, _handle) {

}

template<typename T, typename Entity>
PropertyTT<T,Entity>::PropertyTT(OpenVolumeMeshPropertyT<T> *_prop, ResourceManager &_resMan, PropertyHandleT _handle) :
        PropertyPtr<OpenVolumeMeshPropertyT<T>, Entity>(_prop, _resMan, _handle)
{
}

template<typename T, typename Entity>
BaseProperty *PropertyTT<T,Entity>::clone(ResourceManager &_resMan, OpenVolumeMeshHandle _handle) const
{
    auto prop_clone = ptr::shared_ptr<OpenVolumeMeshPropertyT<T>>::get()->clone();
    return new PropertyTT<T, Entity>(prop_clone, _resMan, PropertyHandleT(_handle.idx()));
}

template<typename T>
const std::string typeName() {
    throw std::runtime_error("Serialization is not supported for these data types!");
}

} // Namespace OpenVolumeMesh
