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

//== INCLUDES =================================================================

#include <cassert>
#include <istream>
#include <ostream>
#include <numeric>
#include <string>
#include <vector>

#include "OpenVolumeMeshBaseProperty.hh"

#include "Serializers.hh"

namespace OpenVolumeMesh {

//== CLASS DEFINITION =========================================================

/** \class OpenVolumeMeshPropertyT
 *
 *  \brief Default property class for any type T.
 *
 *  The default property class for any type T.
 */

template<class T>
class OpenVolumeMeshPropertyT: public OpenVolumeMeshBaseProperty {
public:

    template <class PropT, class Entity> friend class PropertyPtr;

    typedef T                                         Value;
    typedef typename std::vector<T>                   vector_type;
    typedef T                                         value_type;
    typedef typename vector_type::reference           reference;
    typedef typename vector_type::const_reference     const_reference;

public:

	explicit OpenVolumeMeshPropertyT(
            const std::string& _name,
            const std::string& _internal_type_name,
            const T &_def = T())
        : OpenVolumeMeshBaseProperty(_name, _internal_type_name),
          def_(_def)
    {}


	OpenVolumeMeshPropertyT(const OpenVolumeMeshPropertyT& _rhs) = default;

public:
	// inherited from OpenVolumeMeshBaseProperty
	void reserve(size_t _n) override{
		data_.reserve(_n);
	}
	void resize(size_t _n) override {
                data_.resize(_n, def_);
	}
	size_t size() const override {
		return data_.size();
	}
	void clear() override {
		data_.clear();
		vector_type().swap(data_);
	}
	void push_back() override {
		data_.push_back(def_);
	}
	void swap(size_t _i0, size_t _i1) override {
        std::swap(data_[_i0], data_[_i1]);
    }

	virtual void copy(size_t _src_idx, size_t _dst_idx) {
		data_[_dst_idx] = data_[_src_idx];
	}
	void delete_element(size_t _idx) override {
		data_.erase(data_.begin() + static_cast<long>(_idx));
	}

public:

	size_t n_elements() const override {
		return data_.size();
	}
	size_t element_size() const override {
        return sizeof(T);
    }


#ifndef DOXY_IGNORE_THIS
	struct plus {
		size_t operator ()(size_t _b, const T& /*_v*/) {
			return _b + sizeof(T);
		}
	};
#endif

    size_t size_of() const override {
        if (element_size() != OpenVolumeMeshBaseProperty::UnknownSize)
            return this->OpenVolumeMeshBaseProperty::size_of(n_elements());
        return std::accumulate(data_.begin(), data_.end(), size_t(0), plus());
    }

	size_t size_of(size_t _n_elem) const override {
		return this->OpenVolumeMeshBaseProperty::size_of(_n_elem);
	}

	// Function to serialize a property
    void serialize(std::ostream& _ostr) const override {
        for(typename vector_type::const_iterator it = data_.begin();
                it != data_.end(); ++it) {
            OpenVolumeMesh::serialize(_ostr, *it) << std::endl;
        }
    }

    // Function to deserialize a property
    void deserialize(std::istream& _istr) override {
        for(unsigned int i = 0; i < n_elements(); ++i) {
            OpenVolumeMesh::deserialize(_istr, data_[i]);
        }
    }

public:
	// data access interface

	/// Get pointer to array (does not work for T==bool)
	const T* data() const {

		if (data_.empty())
			return 0;

		return &data_[0];
	}

	/// Get reference to property vector (be careful, improper usage, e.g. resizing, may crash)
	vector_type& data_vector() {

		return data_;
	}

	/// Access the i'th element. No range check is performed!
  reference operator[](size_t _idx) {
    assert(_idx < data_.size());
		return data_[_idx];
	}

	/// Const access to the i'th element. No range check is performed!
  const_reference operator[](size_t _idx) const {
    assert(_idx < data_.size());
		return data_[_idx];
	}

	/// Make a copy of self.
	OpenVolumeMeshPropertyT<T>* clone() const override {
		OpenVolumeMeshPropertyT<T>* p = new OpenVolumeMeshPropertyT<T>(*this);
		return p;
	}

	typename vector_type::const_iterator begin() const { return data_.begin(); }

	typename vector_type::iterator begin() { return data_.begin(); }

	typename vector_type::const_iterator end() const { return data_.end(); }

    typename vector_type::iterator end() { return data_.end(); }

protected:

    /// Delete multiple entries in list
    void delete_multiple_entries(const std::vector<bool>& _tags) override {

        assert(_tags.size() == data_.size());
        vector_type new_data;
        typename vector_type::iterator d_it = data_.begin();
        std::vector<bool>::const_iterator t_it = _tags.begin();
        std::vector<bool>::const_iterator t_end = _tags.end();
        for(; t_it != t_end; ++t_it, ++d_it) {
            if(!*t_it) {
                new_data.push_back(*d_it);
            }
        }
        data_.swap(new_data);
    }

private:

	vector_type data_;

	const T def_;
};


//-----------------------------------------------------------------------------
// Property specialization for bool type.
//-----------------------------------------------------------------------------

template<>
inline void OpenVolumeMeshPropertyT<bool>::swap(size_t _i0, size_t _i1)
{
    // std::vector<bool>::swap(reference x, reference y) exists, but
    // on libstdc++ with _GLIBCXX_DEBUG it doesn't compile
    // (2018-02-26, libstdc++ 8.2.0)

    auto tmp = data_[_i0];
    data_[_i0] = data_[_i1];
    data_[_i1] = tmp;
}

template<>
inline size_t OpenVolumeMeshPropertyT<bool>::size_of(size_t _n_elem) const
{
    return _n_elem / 8 + ((_n_elem % 8) != 0);
}

template<>
inline size_t OpenVolumeMeshPropertyT<bool>::size_of() const
{
    return size_of(n_elements());
}

template<>
inline size_t OpenVolumeMeshPropertyT<bool>::element_size() const
{
    return OpenVolumeMeshBaseProperty::UnknownSize;
}

template<>
inline void OpenVolumeMeshPropertyT<bool>::deserialize(std::istream& _istr)
{
    for(unsigned int i = 0; i < n_elements(); ++i) {
        value_type val;
        OpenVolumeMesh::deserialize(_istr, val);
        data_[i] = val;
    }
}


//-----------------------------------------------------------------------------
// Property specialization for std::string type.
//-----------------------------------------------------------------------------
template<>
inline size_t OpenVolumeMeshPropertyT<std::string>::size_of(size_t) const
{
    return OpenVolumeMeshBaseProperty::UnknownSize;
}

template<>
inline size_t OpenVolumeMeshPropertyT<std::string>::size_of() const
{
    return sizeof(data_);
}

template<>
inline size_t OpenVolumeMeshPropertyT<std::string>::element_size() const
{
    return OpenVolumeMeshBaseProperty::UnknownSize;
}


//-----------------------------------------------------------------------------

} // Namespace OpenVolumeMesh

