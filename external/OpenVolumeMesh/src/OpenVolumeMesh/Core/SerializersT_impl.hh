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

#define SERIALIZERST_CC

#include "Serializers.hh"

namespace OpenVolumeMesh
{

template<typename T>
T& decllval();

template <bool B> struct bool_type;
template <>       struct bool_type<true>  { char c[1]; };
template <>       struct bool_type<false> { char c[2]; };

typedef bool_type<true>  true_type;
typedef bool_type<false> false_type;

template <typename Stream, typename T>
class has_input_operator
{
private:
  template<class U> static true_type  test(char(*)[sizeof(decllval<Stream>() >> decllval<U>(), void(), 0)]);
  template<class U> static false_type test(...);

public:
  enum { bool_value = sizeof(true_type) == sizeof(test<T>(nullptr)) };
  typedef bool_type<bool_value> type;
  static type value;
};

template <typename Stream, typename T>
typename has_input_operator<Stream, T>::type has_input_operator<Stream, T>::value = typename has_input_operator<Stream, T>::type();


template <typename Stream, typename T>
class has_output_operator
{
private:
  template<class U> static true_type  test(char(*)[sizeof(decllval<Stream>() << decllval<U>(), void(), 0)]);
  template<class U> static false_type test(...);

public:
  enum { bool_value = sizeof(true_type) == sizeof(test<T>(nullptr)) };
  typedef bool_type<bool_value> type;
  static type value;
};

template <typename Stream, typename T>
typename has_output_operator<Stream, T>::type has_output_operator<Stream, T>::value = typename has_output_operator<Stream, T>::type();


template <typename ValueT>
std::ostream& serialize_helper(std::ostream& _ostr, ValueT& _rhs, true_type)
{
  _ostr << _rhs;
  return _ostr;
}

template <typename ValueT>
std::ostream& serialize_helper(std::ostream& _ostr, ValueT&, false_type)
{
  std::cout << "Warning: trying to serialize a type that does not have a serialize function" << std::endl;
  return _ostr;
}

template <typename ValueT>
std::ostream& serialize(std::ostream& _ostr, const ValueT& _rhs)
{
    return serialize_helper(_ostr, _rhs, has_output_operator<std::ostream, ValueT>::value);
}


template <typename ValueT>
std::istream& deserialize_helper(std::istream& _istr, ValueT& _rhs, true_type)
{
  _istr >> _rhs;
  return _istr;
}

template <typename ValueT>
std::istream& deserialize_helper(std::istream& _istr, ValueT&, false_type)
{
  std::cout << "Warning: trying to deserialize a type that does not have a deserialize function" << std::endl;
  return _istr;
}

template <typename ValueT>
std::istream& deserialize(std::istream& _istr, ValueT& _rhs)
{
  return deserialize_helper(_istr, _rhs, has_input_operator<std::istream, ValueT>::value);
}

template <typename KeyT, typename ValueT>
std::ostream& serialize(std::ostream& os, const std::map< KeyT, ValueT >& rhs)
{
    os << rhs.size() << std::endl;
    for (typename std::map< KeyT, ValueT >::const_iterator it = rhs.begin();
         it != rhs.end();
         ++it)
    {
        serialize(os,it->first) << std::endl;
        serialize(os, it->second) << std::endl;
    }

    return os;
}

template <typename KeyT, typename ValueT>
std::istream& deserialize(std::istream& is, std::map< KeyT, ValueT >& rhs)
{

    size_t size;
    is >> size;
    rhs.clear();
    for (size_t i=0; i<size; i++)
    {
        KeyT key;
        ValueT value;
        deserialize(is, key);
        deserialize(is, value);
        rhs[key] = value;
    }

    return is;
}

template <typename ValueT>
std::ostream& serialize(std::ostream& _ostr, const std::vector< ValueT >& _rhs)
{
    _ostr << _rhs.size() << std::endl;
    for (size_t i = 0; i < _rhs.size(); ++i)
        serialize(_ostr, _rhs[i]) << std::endl;
    return _ostr;
}

template <typename ValueT>
std::istream& deserialize(std::istream& _istr, std::vector< ValueT >& _rhs)
{
    size_t size;
    _istr >> size;
    _rhs.resize(size);
    for (size_t i=0; i<size; i++)
        deserialize(_istr,_rhs[i]);

    return _istr;
}


}
