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
 *   $Revision: 236 $                                                         *
 *   $Date: 2013-02-19 12:32:33 +0100 (Tue, 19 Feb 2013) $                    *
 *   $LastChangedBy: kremer $                                                *
 *                                                                           *
\*===========================================================================*/


#ifndef SERIALIZERS_HH
#define SERIALIZERS_HH

#include <iostream>
#include <map>
#include <vector>

#include <sstream>
#include <string>

#include <iterator>
#include "OpenVolumeMesh/Config/Export.hh"

namespace OpenVolumeMesh
{

template <typename ValueT>
std::ostream& serialize(std::ostream& os, const ValueT& rhs);

OVM_EXPORT
std::ostream& serialize(std::ostream& os, const std::string& rhs);

template <typename ValueT>
std::istream& deserialize(std::istream& is, ValueT& rhs);

OVM_EXPORT
std::istream& deserialize(std::istream& is, std::string& rhs);

template <typename KeyT, typename ValueT>
std::ostream& operator<<(std::ostream& os, const std::map< KeyT, ValueT >& rhs);

template <typename KeyT, typename ValueT>
std::istream& operator>>(std::istream& is, std::map< KeyT, ValueT >& rhs);

template <typename ValueT>
std::ostream& operator<<(std::ostream& os, const std::vector< ValueT >& rhs);

template <typename ValueT>
std::istream& operator>>(std::istream& is, std::vector< ValueT >& rhs);

OVM_EXPORT
std::istream& operator>>(std::istream& is, std::vector< bool >& rhs);

}

#if defined(INCLUDE_TEMPLATES) && !defined(SERIALIZERST_CC)
#include "SerializersT_impl.hh"
#endif

#endif // SERIALIZERS_HH
