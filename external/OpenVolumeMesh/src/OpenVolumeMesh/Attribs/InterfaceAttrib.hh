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

#include "../Core/OpenVolumeMeshProperty.hh"
#include "../Core/OpenVolumeMeshHandle.hh"
#include "../Core/PropertyDefines.hh"
#include "OpenVolumeMesh/Config/Export.hh"
#include "../Core/TopologyKernel.hh"

namespace OpenVolumeMesh {

class TopologyKernel;

/**
 * @brief InterfaceAttrib stores if an entity is part of an interface,
 * e.g. an material boundary inside the volume.
 */
class OVM_EXPORT InterfaceAttrib {
    using boolref = std::vector<bool>::reference;
public:
    explicit InterfaceAttrib(TopologyKernel& _kernel);
    ~InterfaceAttrib() = default;

    bool operator[](const VertexHandle& _h) const {
        return v_interface_[_h];
    }

    boolref operator[](const VertexHandle& _h) {
        return v_interface_[_h];
    }

    bool operator[](const EdgeHandle& _h) const {
        return e_interface_[_h];
    }

    boolref operator[](const EdgeHandle& _h) {
        return e_interface_[_h];
    }

    bool operator[](const HalfEdgeHandle& _h) const {
        return e_interface_[kernel_.edge_handle(_h)];
    }

    boolref operator[](const HalfEdgeHandle& _h) {
        return e_interface_[kernel_.edge_handle(_h)];
    }

    bool operator[](const FaceHandle& _h) const {
        return f_interface_[_h];
    }

    boolref operator[](const FaceHandle& _h) {
        return f_interface_[_h];
    }

    bool operator[](const HalfFaceHandle& _h) const {
        return f_interface_[kernel_.face_handle(_h)];
    }

    boolref operator[](const HalfFaceHandle& _h) {
        return f_interface_[kernel_.face_handle(_h)];
    }

    TopologyKernel& kernel_;

    VertexPropertyT<bool> v_interface_;
    EdgePropertyT<bool> e_interface_;
    FacePropertyT<bool> f_interface_;
};


} // Namespace OpenVolumeMesh
