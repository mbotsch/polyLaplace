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



#include "Serializers.hh"

namespace OpenVolumeMesh
{

std::ostream& serialize(std::ostream& _ostr, const std::string& _rhs)
{
    _ostr << _rhs.size() << ":";
    _ostr << _rhs;

    return _ostr;
}

std::istream& deserialize(std::istream& _istr, std::string& _rhs)
{
    size_t len;
    char delimiter;
    _istr >> len;  //deserialize size of string
    _istr >> delimiter;
    if (_istr && len) {
        std::vector<char> tmp(len);
        _istr.read(&tmp[0] , len); //deserialize characters of string
        _rhs.assign(&tmp[0], len);
    }

    return _istr;
}


std::istream& operator>>(std::istream& _istr, std::vector< bool >& _rhs)
{
    size_t size;
    _istr >> size;
    _rhs.resize(size);
    for (size_t i=0; i<size; i++)
    {
        bool b;
        _istr >> b;
        _rhs[i] = b;
    }

    return _istr;
}

}
