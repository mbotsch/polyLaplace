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
 *   $Revision$                                                         *
 *   $Date$                    *
 *   $LastChangedBy$                                                *
 *                                                                           *
\*===========================================================================*/

#ifndef ITERATORS_HH_
#define ITERATORS_HH_

#include <iterator>
#include <set>
#include <vector>

#ifndef NDEBUG
#include <iostream>
#endif

#include "OpenVolumeMeshHandle.hh"
#include "OpenVolumeMesh/Config/Export.hh"

namespace OpenVolumeMesh {

class TopologyKernel;

template <
class OH /* Output handle type */>
class BaseIterator {
public:

	// STL compliance
	typedef std::bidirectional_iterator_tag iterator_category;
	typedef int						        difference_type;
	typedef const OH				        value_type;
	typedef const OH*				        pointer;
	typedef const OH&				        reference;


    BaseIterator(const TopologyKernel* _mesh, const OH& _ch) :
        valid_(true), cur_handle_(_ch), mesh_(_mesh) {}

    explicit BaseIterator(const TopologyKernel* _mesh) :
        valid_(true), mesh_(_mesh) {}

    // STL compliance (needs to have default constructor)
	BaseIterator() : valid_(false), mesh_(nullptr) {}
    BaseIterator(const BaseIterator& _c) = default;
	virtual ~BaseIterator() = default;

    BaseIterator& operator=(const BaseIterator& _c) = default;

	bool operator== (const BaseIterator& _c) const {
        return (this->cur_handle_ == _c.cur_handle() &&
                this->valid_ == _c.valid() &&
				this->mesh_ == _c.mesh());
	}
	bool operator!= (const BaseIterator& _c) const {
		return !this->operator==(_c);
	}

	pointer operator->() const {
		return &cur_handle_;
	}

	reference operator*() const {
		return cur_handle_;
	}

	bool operator< (const BaseIterator& _c) const {
	    return cur_handle_.idx() < _c.cur_handle_.idx();
	}

	operator bool() const {
		return valid_;
	}

	void valid(bool _valid) {
		valid_ = _valid;
    }
    bool valid() const {
        return valid_;
    }
	void cur_handle(const OH& _h) {
		cur_handle_ = _h;
	}
	reference cur_handle() const {
		return cur_handle_;
    }
	const TopologyKernel* mesh() const {
		return mesh_;
	}

private:

    bool valid_;
    OH cur_handle_;
	const TopologyKernel* mesh_;
};


#if __cplusplus >= 201103L || _MSC_VER >= 1800 // an older MSVC version might be sufficient, didn't test

#include <type_traits>

template<class I>
using is_ovm_iterator = std::is_base_of<BaseIterator<typename std::remove_const<typename I::value_type>::type>, I>;

// provide begin() and end() for the iterator pairs provided in TopologyKernel,
// so we can use range-for, e.g. for(const auto &vh: mesh.vertices()) works.
template<class I>
typename std::enable_if<is_ovm_iterator<I>::value, I>::type
begin(const std::pair<I, I>& iterpair)
{
    return iterpair.first;
}

template<class I>
typename std::enable_if<is_ovm_iterator<I>::value, I>::type
end(const std::pair<I, I>& iterpair)
{
    return iterpair.second;
}

#endif // C++11

template <
class IH /*  Input handle type */,
class OH /* Output handle type */>
class BaseCirculator : public BaseIterator<OH> {
public:

    typedef BaseIterator<OH> BaseIter;

    BaseCirculator(const TopologyKernel* _mesh, const IH& _ih, const OH& _oh, int _max_laps = 1) :
        BaseIter(_mesh, _oh),
        lap_(0),
        max_laps_(_max_laps),
        ref_handle_(_ih)
    {}

    BaseCirculator(const TopologyKernel* _mesh, const IH& _ih, int _max_laps = 1) :
        BaseIter(_mesh, OH()),
        lap_(0),
        max_laps_(_max_laps),
        ref_handle_(_ih)
    {}

    // STL compliance (needs to have default constructor)
    BaseCirculator() :
        BaseIter(),
        lap_(0),
        max_laps_(1)
    {}
    BaseCirculator(const BaseCirculator& _c) = default;

    virtual ~BaseCirculator() = default;

    bool operator== (const BaseCirculator& _c) const {
        return (BaseIter::operator==(_c) &&
                this->lap() == _c.lap() &&
                this->ref_handle() == _c.ref_handle());
    }
    bool operator!= (const BaseCirculator& _c) const {
        return !this->operator==(_c);
    }

    bool operator< (const BaseCirculator& _c) const {
        if (lap_ == _c.lap_)
            return BaseIter::operator<(_c);
        else
            return lap_ < _c.lap_;
    }

    BaseCirculator& operator=(const BaseCirculator& _c) = default;

    const IH& ref_handle() const {
        return ref_handle_;
    }

    void lap(int _lap) {
        lap_ = _lap;
    }
    int lap() const {
        return lap_;
    }

    void max_laps(int _max_laps) {
        max_laps_ = _max_laps;
    }
    int max_laps() const {
        return max_laps_;
    }

protected:
    int lap_;
    int max_laps_;
    IH ref_handle_;

};

//===========================================================================

class OVM_EXPORT VertexOHalfEdgeIter :
    public BaseCirculator<
  VertexHandle,
  HalfEdgeHandle> {
public:
    typedef BaseCirculator<
      VertexHandle,
      HalfEdgeHandle> BaseIter;


  VertexOHalfEdgeIter(const VertexHandle& _vIdx,
            const TopologyKernel* _mesh, int _max_laps = 1);

  // Post increment/decrement operator
  VertexOHalfEdgeIter operator++(int) {
    VertexOHalfEdgeIter cpy = *this;
    ++(*this);
    return cpy;
  }
  VertexOHalfEdgeIter operator--(int) {
    VertexOHalfEdgeIter cpy = *this;
    --(*this);
    return cpy;
  }
  VertexOHalfEdgeIter operator+(int _n) {
    VertexOHalfEdgeIter cpy = *this;
    for(int i = 0; i < _n; ++i) {
      ++cpy;
    }
    return cpy;
  }
  VertexOHalfEdgeIter operator-(int _n) {
    VertexOHalfEdgeIter cpy = *this;
    for(int i = 0; i < _n; ++i) {
      --cpy;
    }
    return cpy;
  }
  VertexOHalfEdgeIter& operator+=(int _n) {
    for(int i = 0; i < _n; ++i) {
      ++(*this);
    }
    return *this;
  }
  VertexOHalfEdgeIter& operator-=(int _n) {
    for(int i = 0; i < _n; ++i) {
      --(*this);
    }
    return *this;
  }

  VertexOHalfEdgeIter& operator++();
  VertexOHalfEdgeIter& operator--();

private:

    size_t cur_index_;
};



//===========================================================================

class OVM_EXPORT VertexVertexIter :
    public BaseCirculator<
  VertexHandle,
  VertexHandle> {
public:
    typedef BaseCirculator<
      VertexHandle,
      VertexHandle> BaseIter;


  VertexVertexIter(const VertexHandle& _vIdx,
            const TopologyKernel* _mesh, int _max_laps = 1);

  // Post increment/decrement operator
  VertexVertexIter operator++(int) {
    VertexVertexIter cpy = *this;
    ++(*this);
    return cpy;
  }
  VertexVertexIter operator--(int) {
    VertexVertexIter cpy = *this;
    --(*this);
    return cpy;
  }
  VertexVertexIter operator+(int _n) {
    VertexVertexIter cpy = *this;
    for(int i = 0; i < _n; ++i) {
      ++cpy;
    }
    return cpy;
  }
  VertexVertexIter operator-(int _n) {
    VertexVertexIter cpy = *this;
    for(int i = 0; i < _n; ++i) {
      --cpy;
    }
    return cpy;
  }
  VertexVertexIter& operator+=(int _n) {
    for(int i = 0; i < _n; ++i) {
      ++(*this);
    }
    return *this;
  }
  VertexVertexIter& operator-=(int _n) {
    for(int i = 0; i < _n; ++i) {
      --(*this);
    }
    return *this;
  }

  VertexVertexIter& operator++();
  VertexVertexIter& operator--();

private:

    size_t cur_index_;
};

//===========================================================================

class OVM_EXPORT HalfEdgeHalfFaceIter : public BaseCirculator<
	HalfEdgeHandle,
	HalfFaceHandle> {
public:
    typedef BaseCirculator<
			HalfEdgeHandle,
			HalfFaceHandle> BaseIter;


    HalfEdgeHalfFaceIter(const HalfEdgeHandle& _heIdx, const TopologyKernel* _mesh, int _max_laps = 1);

	// Post increment/decrement operator
	HalfEdgeHalfFaceIter operator++(int) {
		HalfEdgeHalfFaceIter cpy = *this;
		++(*this);
		return cpy;
	}
	HalfEdgeHalfFaceIter operator--(int) {
		HalfEdgeHalfFaceIter cpy = *this;
		--(*this);
		return cpy;
	}
	HalfEdgeHalfFaceIter operator+(int _n) {
		HalfEdgeHalfFaceIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			++cpy;
		}
		return cpy;
	}
	HalfEdgeHalfFaceIter operator-(int _n) {
		HalfEdgeHalfFaceIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			--cpy;
		}
		return cpy;
	}
	HalfEdgeHalfFaceIter& operator+=(int _n) {
		for(int i = 0; i < _n; ++i) {
			++(*this);
		}
		return *this;
	}
	HalfEdgeHalfFaceIter& operator-=(int _n) {
		for(int i = 0; i < _n; ++i) {
			--(*this);
		}
		return *this;
	}

	HalfEdgeHalfFaceIter& operator++();
	HalfEdgeHalfFaceIter& operator--();

private:
    size_t cur_index_;
};

//===========================================================================

class OVM_EXPORT VertexFaceIter : public BaseCirculator<
  VertexHandle,
  FaceHandle> {
public:
    typedef BaseCirculator<
      VertexHandle,
      FaceHandle> BaseIter;

    VertexFaceIter(const VertexHandle& _vIdx, const TopologyKernel* _mesh, int _max_laps = 1);

  // Post increment/decrement operator
  VertexFaceIter operator++(int) {
    VertexFaceIter cpy = *this;
    ++(*this);
    return cpy;
  }
  VertexFaceIter operator--(int) {
    VertexFaceIter cpy = *this;
    --(*this);
    return cpy;
  }
  VertexFaceIter operator+(int _n) {
    VertexFaceIter cpy = *this;
    for(int i = 0; i < _n; ++i) {
      ++cpy;
    }
    return cpy;
  }
  VertexFaceIter operator-(int _n) {
    VertexFaceIter cpy = *this;
    for(int i = 0; i < _n; ++i) {
      --cpy;
    }
    return cpy;
  }
  VertexFaceIter& operator+=(int _n) {
    for(int i = 0; i < _n; ++i) {
      ++(*this);
    }
    return *this;
  }
  VertexFaceIter& operator-=(int _n) {
    for(int i = 0; i < _n; ++i) {
      --(*this);
    }
    return *this;
  }

  VertexFaceIter& operator++();
  VertexFaceIter& operator--();

private:
    std::vector<FaceHandle> faces_;
    size_t cur_index_;
};

//===========================================================================

class OVM_EXPORT VertexCellIter : public BaseCirculator<
	VertexHandle,
	CellHandle> {
public:
    typedef BaseCirculator<
			VertexHandle,
			CellHandle> BaseIter;

    VertexCellIter(const VertexHandle& _vIdx, const TopologyKernel* _mesh, int _max_laps = 1);

	// Post increment/decrement operator
	VertexCellIter operator++(int) {
		VertexCellIter cpy = *this;
		++(*this);
		return cpy;
	}
	VertexCellIter operator--(int) {
		VertexCellIter cpy = *this;
		--(*this);
		return cpy;
	}
	VertexCellIter operator+(int _n) {
		VertexCellIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			++cpy;
		}
		return cpy;
	}
	VertexCellIter operator-(int _n) {
		VertexCellIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			--cpy;
		}
		return cpy;
	}
	VertexCellIter& operator+=(int _n) {
		for(int i = 0; i < _n; ++i) {
			++(*this);
		}
		return *this;
	}
	VertexCellIter& operator-=(int _n) {
		for(int i = 0; i < _n; ++i) {
			--(*this);
		}
		return *this;
	}

	VertexCellIter& operator++();
	VertexCellIter& operator--();

private:
    std::vector<CellHandle> cells_;
    size_t cur_index_;
};

class OVM_EXPORT HalfEdgeCellIter : public BaseCirculator<
	HalfEdgeHandle,
	CellHandle> {
public:
    typedef BaseCirculator<
			HalfEdgeHandle,
			CellHandle> BaseIter;


    HalfEdgeCellIter(const HalfEdgeHandle& _heIdx, const TopologyKernel* _mesh, int _max_laps = 1);

	// Post increment/decrement operator
	HalfEdgeCellIter operator++(int) {
		HalfEdgeCellIter cpy = *this;
		++(*this);
		return cpy;
	}
	HalfEdgeCellIter operator--(int) {
		HalfEdgeCellIter cpy = *this;
		--(*this);
        return cpy;
    }
	HalfEdgeCellIter operator+(int _n) {
		HalfEdgeCellIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			++cpy;
		}
		return cpy;
	}
	HalfEdgeCellIter operator-(int _n) {
		HalfEdgeCellIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			--cpy;
		}
		return cpy;
	}
	HalfEdgeCellIter& operator+=(int _n) {
		for(int i = 0; i < _n; ++i) {
			++(*this);
		}
		return *this;
	}
	HalfEdgeCellIter& operator-=(int _n) {
		for(int i = 0; i < _n; ++i) {
			--(*this);
		}
		return *this;
	}

	HalfEdgeCellIter& operator++();
	HalfEdgeCellIter& operator--();

private:
    CellHandle getCellHandle(int _cur_index) const;

private:
    std::vector<CellHandle> cells_;
    size_t cur_index_;
};

//===========================================================================

class OVM_EXPORT CellVertexIter : public BaseCirculator<
	CellHandle,
	VertexHandle> {
public:
    typedef BaseCirculator<
			CellHandle,
			VertexHandle> BaseIter;

    CellVertexIter(const CellHandle& _cIdx, const TopologyKernel* _mesh, int _max_laps = 1);

	// Post increment/decrement operator
	CellVertexIter operator++(int) {
		CellVertexIter cpy = *this;
		++(*this);
		return cpy;
	}
	CellVertexIter operator--(int) {
		CellVertexIter cpy = *this;
		--(*this);
		return cpy;
	}
	CellVertexIter operator+(int _n) {
		CellVertexIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			++cpy;
		}
		return cpy;
	}
	CellVertexIter operator-(int _n) {
		CellVertexIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			--cpy;
		}
		return cpy;
	}
	CellVertexIter& operator+=(int _n) {
		for(int i = 0; i < _n; ++i) {
			++(*this);
		}
		return *this;
	}
	CellVertexIter& operator-=(int _n) {
		for(int i = 0; i < _n; ++i) {
			--(*this);
		}
		return *this;
	}

	CellVertexIter& operator++();
	CellVertexIter& operator--();

private:
	std::vector<VertexHandle> incident_vertices_;
    size_t cur_index_;
};

//===========================================================================

class OVM_EXPORT CellCellIter : public BaseCirculator<
	CellHandle,
	CellHandle> {
public:
    typedef BaseCirculator<
			CellHandle,
			CellHandle> BaseIter;

    CellCellIter(const CellHandle& _cIdx, const TopologyKernel* _mesh, int _max_laps = 1);

	// Post increment/decrement operator
	CellCellIter operator++(int) {
		CellCellIter cpy = *this;
		++(*this);
		return cpy;
	}
	CellCellIter operator--(int) {
		CellCellIter cpy = *this;
		--(*this);
		return cpy;
	}
	CellCellIter operator+(int _n) {
		CellCellIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			++cpy;
		}
		return cpy;
	}
	CellCellIter operator-(int _n) {
		CellCellIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			--cpy;
		}
		return cpy;
	}
	CellCellIter& operator+=(int _n) {
		for(int i = 0; i < _n; ++i) {
			++(*this);
		}
		return *this;
	}
	CellCellIter& operator-=(int _n) {
		for(int i = 0; i < _n; ++i) {
			--(*this);
		}
		return *this;
	}

	CellCellIter& operator++();
	CellCellIter& operator--();

private:
    std::vector<CellHandle> adjacent_cells_;
    size_t cur_index_;
};

//===========================================================================

class OVM_EXPORT HalfFaceVertexIter : public BaseCirculator<
    HalfFaceHandle,
    VertexHandle> {
public:
    typedef BaseCirculator<
            HalfFaceHandle,
            VertexHandle> BaseIter;

    HalfFaceVertexIter(const HalfFaceHandle& _hIdx, const TopologyKernel* _mesh, int _max_laps = 1);

    // Post increment/decrement operator
    HalfFaceVertexIter operator++(int) {
        HalfFaceVertexIter cpy = *this;
        ++(*this);
        return cpy;
    }
    HalfFaceVertexIter operator--(int) {
        HalfFaceVertexIter cpy = *this;
        --(*this);
        return cpy;
    }
    HalfFaceVertexIter operator+(int _n) {
        HalfFaceVertexIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            ++cpy;
        }
        return cpy;
    }
    HalfFaceVertexIter operator-(int _n) {
        HalfFaceVertexIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            --cpy;
        }
        return cpy;
    }
    HalfFaceVertexIter& operator+=(int _n) {
        for(int i = 0; i < _n; ++i) {
            ++(*this);
        }
        return *this;
    }
    HalfFaceVertexIter& operator-=(int _n) {
        for(int i = 0; i < _n; ++i) {
            --(*this);
        }
        return *this;
    }

    HalfFaceVertexIter& operator++();
    HalfFaceVertexIter& operator--();

private:
    std::vector<VertexHandle> vertices_;
    size_t cur_index_;
};

//===========================================================================

class OVM_EXPORT BoundaryHalfFaceHalfFaceIter : public BaseCirculator<HalfFaceHandle,
    HalfFaceHandle> {
private:
    typedef BaseCirculator<HalfFaceHandle,
            HalfFaceHandle> BaseIter;
public:
    BoundaryHalfFaceHalfFaceIter(const HalfFaceHandle& _ref_h,
            const TopologyKernel* _mesh, int _max_laps = 1);

    // Post increment/decrement operator
    BoundaryHalfFaceHalfFaceIter operator++(int) {
        BoundaryHalfFaceHalfFaceIter cpy = *this;
        ++(*this);
        return cpy;
    }
    BoundaryHalfFaceHalfFaceIter operator--(int) {
        BoundaryHalfFaceHalfFaceIter cpy = *this;
        --(*this);
        return cpy;
    }
    BoundaryHalfFaceHalfFaceIter operator+(int _n) {
        BoundaryHalfFaceHalfFaceIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            ++cpy;
        }
        return cpy;
    }
    BoundaryHalfFaceHalfFaceIter operator-(int _n) {
        BoundaryHalfFaceHalfFaceIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            --cpy;
        }
        return cpy;
    }
    BoundaryHalfFaceHalfFaceIter& operator+=(int _n) {
        for(int i = 0; i < _n; ++i) {
            ++(*this);
        }
        return *this;
    }
    BoundaryHalfFaceHalfFaceIter& operator-=(int _n) {
        for(int i = 0; i < _n; ++i) {
            --(*this);
        }
        return *this;
    }

    const EdgeHandle& common_edge() const { return common_edges_[cur_index_]; }

    BoundaryHalfFaceHalfFaceIter& operator++();
    BoundaryHalfFaceHalfFaceIter& operator--();

private:
    std::vector<HalfFaceHandle> neighbor_halffaces_;
    std::vector<EdgeHandle> common_edges_;
    size_t cur_index_;
};

//===========================================================================

class OVM_EXPORT VertexIter : public BaseIterator<VertexHandle> {
public:
    typedef BaseIterator<VertexHandle> BaseIter;


    VertexIter(const TopologyKernel* _mesh, const VertexHandle& _vh = VertexHandle(0));

	// Post increment/decrement operator
	VertexIter operator++(int) {
		VertexIter cpy = *this;
		++(*this);
		return cpy;
	}
	VertexIter operator--(int) {
		VertexIter cpy = *this;
		--(*this);
		return cpy;
	}
	VertexIter operator+(int _n) {
		VertexIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			++cpy;
		}
		return cpy;
	}
	VertexIter operator-(int _n) {
		VertexIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			--cpy;
		}
		return cpy;
	}
	VertexIter& operator+=(int _n) {
		for(int i = 0; i < _n; ++i) {
			++(*this);
		}
		return *this;
	}
	VertexIter& operator-=(int _n) {
		for(int i = 0; i < _n; ++i) {
			--(*this);
		}
		return *this;
	}

	VertexIter& operator++();
	VertexIter& operator--();

private:
    int cur_index_;
};

//===========================================================================

class OVM_EXPORT EdgeIter : public BaseIterator<EdgeHandle> {
public:
    typedef BaseIterator<EdgeHandle> BaseIter;


	EdgeIter(const TopologyKernel* _mesh, const EdgeHandle& _eh = EdgeHandle(0));

	// Post increment/decrement operator
	EdgeIter operator++(int) {
		EdgeIter cpy = *this;
		++(*this);
		return cpy;
	}
	EdgeIter operator--(int) {
		EdgeIter cpy = *this;
		--(*this);
		return cpy;
	}
	EdgeIter operator+(int _n) {
		EdgeIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			++cpy;
		}
		return cpy;
	}
	EdgeIter operator-(int _n) {
		EdgeIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			--cpy;
		}
		return cpy;
	}
	EdgeIter& operator+=(int _n) {
		for(int i = 0; i < _n; ++i) {
			++(*this);
		}
		return *this;
	}
	EdgeIter& operator-=(int _n) {
		for(int i = 0; i < _n; ++i) {
			--(*this);
		}
		return *this;
	}

	EdgeIter& operator++();
	EdgeIter& operator--();

private:
    int cur_index_;
};

//===========================================================================

class OVM_EXPORT HalfEdgeIter : public BaseIterator<HalfEdgeHandle> {
public:
    typedef BaseIterator<HalfEdgeHandle> BaseIter;


	HalfEdgeIter(const TopologyKernel* _mesh, const HalfEdgeHandle& _heh = HalfEdgeHandle(0));

	// Post increment/decrement operator
	HalfEdgeIter operator++(int) {
		HalfEdgeIter cpy = *this;
		++(*this);
		return cpy;
	}
	HalfEdgeIter operator--(int) {
		HalfEdgeIter cpy = *this;
		--(*this);
		return cpy;
	}
	HalfEdgeIter operator+(int _n) {
		HalfEdgeIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			++cpy;
		}
		return cpy;
	}
	HalfEdgeIter operator-(int _n) {
		HalfEdgeIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			--cpy;
		}
		return cpy;
	}
	HalfEdgeIter& operator+=(int _n) {
		for(int i = 0; i < _n; ++i) {
			++(*this);
		}
		return *this;
	}
	HalfEdgeIter& operator-=(int _n) {
		for(int i = 0; i < _n; ++i) {
			--(*this);
		}
		return *this;
	}

	HalfEdgeIter& operator++();
	HalfEdgeIter& operator--();

private:
    int cur_index_;
};

//===========================================================================

class OVM_EXPORT FaceIter : public BaseIterator<FaceHandle> {
public:
    typedef BaseIterator<FaceHandle> BaseIter;


	FaceIter(const TopologyKernel* _mesh, const FaceHandle& _fh = FaceHandle(0));

	// Post increment/decrement operator
	FaceIter operator++(int) {
		FaceIter cpy = *this;
		++(*this);
		return cpy;
	}
	FaceIter operator--(int) {
		FaceIter cpy = *this;
		--(*this);
		return cpy;
	}
	FaceIter operator+(int _n) {
		FaceIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			++cpy;
		}
		return cpy;
	}
	FaceIter operator-(int _n) {
		FaceIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			--cpy;
		}
		return cpy;
	}
	FaceIter& operator+=(int _n) {
		for(int i = 0; i < _n; ++i) {
			++(*this);
		}
		return *this;
	}
	FaceIter& operator-=(int _n) {
		for(int i = 0; i < _n; ++i) {
			--(*this);
		}
		return *this;
	}

	FaceIter& operator++();
	FaceIter& operator--();

private:
    int cur_index_;
};

//===========================================================================

class OVM_EXPORT HalfFaceIter : public BaseIterator<HalfFaceHandle> {
public:
    typedef BaseIterator<HalfFaceHandle> BaseIter;


	HalfFaceIter(const TopologyKernel* _mesh, const HalfFaceHandle& _hfh = HalfFaceHandle(0));

	// Post increment/decrement operator
	HalfFaceIter operator++(int) {
		HalfFaceIter cpy = *this;
		++(*this);
		return cpy;
	}
	HalfFaceIter operator--(int) {
		HalfFaceIter cpy = *this;
		--(*this);
		return cpy;
	}
	HalfFaceIter operator+(int _n) {
		HalfFaceIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			++cpy;
		}
		return cpy;
	}
	HalfFaceIter operator-(int _n) {
		HalfFaceIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			--cpy;
		}
		return cpy;
	}
	HalfFaceIter& operator+=(int _n) {
		for(int i = 0; i < _n; ++i) {
			++(*this);
		}
		return *this;
	}
	HalfFaceIter& operator-=(int _n) {
		for(int i = 0; i < _n; ++i) {
			--(*this);
		}
		return *this;
	}

	HalfFaceIter& operator++();
	HalfFaceIter& operator--();

private:
    int cur_index_;
};

//===========================================================================

class OVM_EXPORT CellIter : public BaseIterator<CellHandle> {
public:
    typedef BaseIterator<CellHandle> BaseIter;


	CellIter(const TopologyKernel* _mesh, const CellHandle& _ch = CellHandle(0));

	// Post increment/decrement operator
	CellIter operator++(int) {
		CellIter cpy = *this;
		++(*this);
		return cpy;
	}
	CellIter operator--(int) {
		CellIter cpy = *this;
		--(*this);
		return cpy;
	}
	CellIter operator+(int _n) {
		CellIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			++cpy;
		}
		return cpy;
	}
	CellIter operator-(int _n) {
		CellIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			--cpy;
		}
		return cpy;
	}
	CellIter& operator+=(int _n) {
		for(int i = 0; i < _n; ++i) {
			++(*this);
		}
		return *this;
	}
	CellIter& operator-=(int _n) {
		for(int i = 0; i < _n; ++i) {
			--(*this);
		}
		return *this;
	}

	CellIter& operator++();
	CellIter& operator--();

private:
    int cur_index_;
};

//===========================================================================

namespace Internal {

//===========================================================================

class OVM_EXPORT VertexIHalfEdgeIterImpl : public BaseCirculator<VertexHandle, HalfEdgeHandle> {
public:

    typedef BaseCirculator<VertexHandle, HalfEdgeHandle> BaseIter;
    typedef VertexHandle CenterEntityHandle;

    VertexIHalfEdgeIterImpl(const VertexHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps = 1);

    VertexIHalfEdgeIterImpl& operator++();
    VertexIHalfEdgeIterImpl& operator--();

private:
    VertexOHalfEdgeIter voh_iter_;
};

//===========================================================================

class OVM_EXPORT VertexEdgeIterImpl : public BaseCirculator<VertexHandle, EdgeHandle> {
public:

    typedef BaseCirculator<VertexHandle, EdgeHandle> BaseIter;
    typedef VertexHandle CenterEntityHandle;

    VertexEdgeIterImpl(const VertexHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps = 1);

    VertexEdgeIterImpl& operator++();
    VertexEdgeIterImpl& operator--();

private:
    VertexOHalfEdgeIter voh_iter_;
};

//===========================================================================

class OVM_EXPORT VertexHalfFaceIterImpl : public BaseCirculator<VertexHandle, HalfFaceHandle> {
public:

    typedef BaseCirculator<VertexHandle, HalfFaceHandle> BaseIter;
    typedef VertexHandle CenterEntityHandle;

    VertexHalfFaceIterImpl(const VertexHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps = 1);

    VertexHalfFaceIterImpl& operator++();
    VertexHalfFaceIterImpl& operator--();

private:
    std::vector<HalfFaceHandle> halffaces_;
    size_t cur_index_;
};

//===========================================================================

class OVM_EXPORT HalfEdgeFaceIterImpl : public BaseCirculator<HalfEdgeHandle, FaceHandle> {
public:

    typedef BaseCirculator<HalfEdgeHandle, FaceHandle> BaseIter;
    typedef HalfEdgeHandle CenterEntityHandle;

    HalfEdgeFaceIterImpl(const HalfEdgeHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps = 1);

    HalfEdgeFaceIterImpl& operator++();
    HalfEdgeFaceIterImpl& operator--();

private:
    std::vector<FaceHandle> faces_;
    size_t cur_index_;
};

//===========================================================================

class OVM_EXPORT EdgeHalfFaceIterImpl : public BaseCirculator<EdgeHandle, HalfFaceHandle> {
public:

    typedef BaseCirculator<EdgeHandle, HalfFaceHandle> BaseIter;
    typedef EdgeHandle CenterEntityHandle;

    EdgeHalfFaceIterImpl(const EdgeHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps = 1);

    EdgeHalfFaceIterImpl& operator++();
    EdgeHalfFaceIterImpl& operator--();

private:
    std::vector<HalfFaceHandle> halffaces_;
    size_t cur_index_;
};

//===========================================================================

class OVM_EXPORT EdgeFaceIterImpl : public HalfEdgeFaceIterImpl {
public:

    typedef EdgeHandle CenterEntityHandle;
    EdgeFaceIterImpl(const EdgeHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps = 1);

};

//===========================================================================

class OVM_EXPORT EdgeCellIterImpl : public HalfEdgeCellIter {
public:

    typedef EdgeHandle CenterEntityHandle;
    EdgeCellIterImpl(const EdgeHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps = 1);

};

//===========================================================================

class OVM_EXPORT HalfFaceHalfEdgeIterImpl : public BaseCirculator<HalfFaceHandle, HalfEdgeHandle> {
public:

    typedef BaseCirculator<HalfFaceHandle, HalfEdgeHandle> BaseIter;
    typedef HalfFaceHandle CenterEntityHandle;

    HalfFaceHalfEdgeIterImpl(const HalfFaceHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps = 1);

    HalfFaceHalfEdgeIterImpl& operator++();
    HalfFaceHalfEdgeIterImpl& operator--();

private:
    size_t cur_index_;
};

//===========================================================================

class OVM_EXPORT HalfFaceEdgeIterImpl : public BaseCirculator<HalfFaceHandle, EdgeHandle> {
public:

    typedef BaseCirculator<HalfFaceHandle, EdgeHandle> BaseIter;
    typedef HalfFaceHandle CenterEntityHandle;

    HalfFaceEdgeIterImpl(const HalfFaceHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps = 1);

    HalfFaceEdgeIterImpl& operator++();
    HalfFaceEdgeIterImpl& operator--();

private:
    size_t cur_index_;
};

//===========================================================================

class OVM_EXPORT FaceVertexIterImpl : public HalfFaceVertexIter {
public:

    typedef FaceHandle CenterEntityHandle;
    FaceVertexIterImpl(const FaceHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps = 1);

};

//===========================================================================

class OVM_EXPORT FaceHalfEdgeIterImpl : public HalfFaceHalfEdgeIterImpl {
public:

    typedef FaceHandle CenterEntityHandle;
    FaceHalfEdgeIterImpl(const FaceHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps = 1);

};

//===========================================================================

class OVM_EXPORT FaceEdgeIterImpl : public HalfFaceEdgeIterImpl {
public:

    typedef FaceHandle CenterEntityHandle;
    FaceEdgeIterImpl(const FaceHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps = 1);

};

//===========================================================================

class OVM_EXPORT CellHalfEdgeIterImpl : public BaseCirculator<CellHandle, HalfEdgeHandle> {
public:

    typedef BaseCirculator<CellHandle, HalfEdgeHandle> BaseIter;
    typedef CellHandle CenterEntityHandle;

    CellHalfEdgeIterImpl(const CellHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps = 1);

    CellHalfEdgeIterImpl& operator++();
    CellHalfEdgeIterImpl& operator--();

private:
    std::vector<HalfEdgeHandle> halfedges_;
    size_t cur_index_;
};

//===========================================================================

class OVM_EXPORT CellEdgeIterImpl : public BaseCirculator<CellHandle, EdgeHandle> {
public:

    typedef BaseCirculator<CellHandle, EdgeHandle> BaseIter;
    typedef CellHandle CenterEntityHandle;

    CellEdgeIterImpl(const CellHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps = 1);

    CellEdgeIterImpl& operator++();
    CellEdgeIterImpl& operator--();

private:
    std::vector<EdgeHandle> edges_;
    size_t cur_index_;
};

//===========================================================================

class OVM_EXPORT CellHalfFaceIterImpl : public BaseCirculator<CellHandle, HalfFaceHandle> {
public:

    typedef BaseCirculator<CellHandle, HalfFaceHandle> BaseIter;
    typedef CellHandle CenterEntityHandle;

    CellHalfFaceIterImpl(const CellHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps = 1);

    CellHalfFaceIterImpl& operator++();
    CellHalfFaceIterImpl& operator--();

private:
    std::vector<HalfFaceHandle>::const_iterator hf_iter_;
};

//===========================================================================

class OVM_EXPORT CellFaceIterImpl : public BaseCirculator<CellHandle, FaceHandle> {
public:

    typedef BaseCirculator<CellHandle, FaceHandle> BaseIter;
    typedef CellHandle CenterEntityHandle;

    CellFaceIterImpl(const CellHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps = 1);

    CellFaceIterImpl& operator++();
    CellFaceIterImpl& operator--();

private:
    std::vector<HalfFaceHandle>::const_iterator hf_iter_;
};

//===========================================================================

} // Namespace Internal

//===========================================================================

template <class CirculatorImpl>
class GenericCirculator : public CirculatorImpl {
public:

    GenericCirculator(const typename CirculatorImpl::CenterEntityHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps = 1) :
        CirculatorImpl(_ref_h, _mesh, _max_laps) {}

    GenericCirculator& operator++() {
        CirculatorImpl::operator++();
        return *this;
    }

    GenericCirculator& operator--() {
        CirculatorImpl::operator--();
        return *this;
    }

    // Post increment/decrement operator
    GenericCirculator operator++(int) {
        GenericCirculator cpy = *this;
        ++(*this);
        return cpy;
    }
    GenericCirculator operator--(int) {
        GenericCirculator cpy = *this;
        --(*this);
        return cpy;
    }
    GenericCirculator operator+(int _n) {
        GenericCirculator cpy = *this;
        for(int i = 0; i < _n; ++i) {
            ++cpy;
        }
        return cpy;
    }
    GenericCirculator operator-(int _n) {
        GenericCirculator cpy = *this;
        for(int i = 0; i < _n; ++i) {
            --cpy;
        }
        return cpy;
    }
    GenericCirculator& operator+=(int _n) {
        for(int i = 0; i < _n; ++i) {
            ++(*this);
        }
        return *this;
    }
    GenericCirculator& operator-=(int _n) {
        for(int i = 0; i < _n; ++i) {
            --(*this);
        }
        return *this;
    }

};

//===========================================================================


template class OVM_EXPORT GenericCirculator<Internal::VertexIHalfEdgeIterImpl>;
template class OVM_EXPORT GenericCirculator<Internal::VertexEdgeIterImpl>;
template class OVM_EXPORT GenericCirculator<Internal::VertexHalfFaceIterImpl>;

template class OVM_EXPORT GenericCirculator<Internal::HalfEdgeFaceIterImpl>;

template class OVM_EXPORT GenericCirculator<Internal::EdgeHalfFaceIterImpl>;
template class OVM_EXPORT GenericCirculator<Internal::EdgeFaceIterImpl>;
template class OVM_EXPORT GenericCirculator<Internal::EdgeCellIterImpl>;

template class OVM_EXPORT GenericCirculator<Internal::HalfFaceHalfEdgeIterImpl>;
template class OVM_EXPORT GenericCirculator<Internal::HalfFaceEdgeIterImpl>;

template class OVM_EXPORT GenericCirculator<Internal::FaceVertexIterImpl>;
template class OVM_EXPORT GenericCirculator<Internal::FaceHalfEdgeIterImpl>;
template class OVM_EXPORT GenericCirculator<Internal::FaceEdgeIterImpl>;

template class OVM_EXPORT GenericCirculator<Internal::CellHalfEdgeIterImpl>;
template class OVM_EXPORT GenericCirculator<Internal::CellEdgeIterImpl>;
template class OVM_EXPORT GenericCirculator<Internal::CellHalfFaceIterImpl>;
template class OVM_EXPORT GenericCirculator<Internal::CellFaceIterImpl>;





typedef GenericCirculator<Internal::VertexIHalfEdgeIterImpl> VertexIHalfEdgeIter;
typedef GenericCirculator<Internal::VertexEdgeIterImpl> VertexEdgeIter;
typedef GenericCirculator<Internal::VertexHalfFaceIterImpl> VertexHalfFaceIter;

typedef GenericCirculator<Internal::HalfEdgeFaceIterImpl> HalfEdgeFaceIter;

typedef GenericCirculator<Internal::EdgeHalfFaceIterImpl> EdgeHalfFaceIter;
typedef GenericCirculator<Internal::EdgeFaceIterImpl> EdgeFaceIter;
typedef GenericCirculator<Internal::EdgeCellIterImpl> EdgeCellIter;

typedef GenericCirculator<Internal::HalfFaceHalfEdgeIterImpl> HalfFaceHalfEdgeIter;
typedef GenericCirculator<Internal::HalfFaceEdgeIterImpl> HalfFaceEdgeIter;

typedef GenericCirculator<Internal::FaceVertexIterImpl> FaceVertexIter;
typedef GenericCirculator<Internal::FaceHalfEdgeIterImpl> FaceHalfEdgeIter;
typedef GenericCirculator<Internal::FaceEdgeIterImpl> FaceEdgeIter;

typedef GenericCirculator<Internal::CellHalfEdgeIterImpl> CellHalfEdgeIter;
typedef GenericCirculator<Internal::CellEdgeIterImpl> CellEdgeIter;
typedef GenericCirculator<Internal::CellHalfFaceIterImpl> CellHalfFaceIter;
typedef GenericCirculator<Internal::CellFaceIterImpl> CellFaceIter;

//===========================================================================

template <class Iter, class Handle>
class BoundaryItemIter : public BaseIterator<Handle> {
public:
    typedef BaseIterator<Handle> BaseIter;


    explicit BoundaryItemIter(const TopologyKernel* _mesh) :
    BaseIter(_mesh),
    it_(_mesh, Handle(0)),
    it_begin_(_mesh, Handle(0)),
    it_end_(_mesh, Handle((int)n_items())) {

        if(!has_incidences()) {
    #ifndef NDEBUG
            std::cerr << "This iterator needs bottom-up incidences!" << std::endl;
    #endif
            BaseIter::valid(false);
            return;
        }

        while(it_ != it_end_ && !BaseIter::mesh()->is_boundary(*it_)){
            ++it_;
        }
        BaseIter::valid(it_ != it_end_);
        if(BaseIter::valid()) {
            BaseIter::cur_handle(*it_);
        }
    }

    // Post increment/decrement operator
    BoundaryItemIter operator++(int) {
        BoundaryItemIter cpy = *this;
        ++(*this);
        return cpy;
    }
    BoundaryItemIter operator--(int) {
        BoundaryItemIter cpy = *this;
        --(*this);
        return cpy;
    }
    BoundaryItemIter operator+(int _n) {
        BoundaryItemIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            ++cpy;
        }
        return cpy;
    }
    BoundaryItemIter operator-(int _n) {
        BoundaryItemIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            --cpy;
        }
        return cpy;
    }
    BoundaryItemIter& operator+=(int _n) {
        for(int i = 0; i < _n; ++i) {
            ++(*this);
        }
        return *this;
    }
    BoundaryItemIter& operator-=(int _n) {
        for(int i = 0; i < _n; ++i) {
            --(*this);
        }
        return *this;
    }

    BoundaryItemIter& operator--() {
        --it_;
        while(it_ >= it_begin_ && !BaseIter::mesh()->is_boundary(*it_)){
            --it_;
        }
        if(it_ >= it_begin_) {
            BaseIter::cur_handle(*it_);
        } else {
            BaseIter::valid(false);
        }
        return *this;
    }

    BoundaryItemIter& operator++() {
        ++it_;
        while(it_ != it_end_ && !BaseIter::mesh()->is_boundary(*it_)){
            ++it_;
        }
        if(it_ != it_end_) {
            BaseIter::cur_handle(*it_);
        } else {
            BaseIter::valid(false);
        }
        return *this;
    }

private:
    size_t n_items() const;
    bool has_incidences() const;

private:
    Iter it_;
    const Iter it_begin_;
    const Iter it_end_;
};

//===========================================================================
typedef BoundaryItemIter<VertexIter, VertexHandle> BoundaryVertexIter;
typedef BoundaryItemIter<HalfEdgeIter, HalfEdgeHandle> BoundaryHalfEdgeIter;
typedef BoundaryItemIter<EdgeIter, EdgeHandle> BoundaryEdgeIter;
typedef BoundaryItemIter<HalfFaceIter, HalfFaceHandle> BoundaryHalfFaceIter;
typedef BoundaryItemIter<FaceIter, FaceHandle> BoundaryFaceIter;
typedef BoundaryItemIter<CellIter, CellHandle> BoundaryCellIter;

//===========================================================================

} // Namespace OpenVolumeMesh

#endif /* ITERATORS_HH_ */
