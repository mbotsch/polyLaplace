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

#include <istream>

#include "OpenVolumeMeshHandle.hh"

namespace OpenVolumeMesh {

bool operator==(const int& _lhs, const OpenVolumeMeshHandle& _rhs) {
    return _lhs == _rhs.idx();
}

bool operator==(const unsigned int& _lhs, const OpenVolumeMeshHandle& _rhs) {

    return _lhs == _rhs.uidx();
}

bool operator!=(const int& _lhs, const OpenVolumeMeshHandle& _rhs) {

    return !(_lhs == _rhs);
}

bool operator!=(const unsigned int& _lhs, const OpenVolumeMeshHandle& _rhs) {

    return !(_lhs == _rhs);
}

std::ostream& operator<<(std::ostream& _ostr, const OpenVolumeMeshHandle& _handle) {
    _ostr << _handle.idx();
    return _ostr;
}

std::istream& operator>>(std::istream& _istr, OpenVolumeMeshHandle& _handle) {
    int val = 0;
    _istr >> val;
    _handle.idx(val);
    return _istr;
}

} // Namespace OpenVolumeMesh
