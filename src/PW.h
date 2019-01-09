/***************************************************************************
 *  Thunderegg, a library for solving Poisson's equation on adaptively 
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  Thunderegg Developers. See AUTHORS.md file at the
 *  top-level directory.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 ***************************************************************************/

#ifndef PW_H
#define PW_H
#include <petscsys.h>
template <typename X> class PW
{
	protected:
	X       obj       = nullptr;
	size_t *ref_count = nullptr;
	bool    owns      = true;
	void    decrement()
	{
		if (ref_count != nullptr) {
			--(*ref_count);
			if (*ref_count == 0) {
				if (obj != nullptr && owns) { PetscObjectDestroy((PetscObject *) &obj); }
				delete ref_count;
			}
		}
	}

	public:
	explicit PW(X obj=nullptr, bool owns=true)
	{
        this->obj=obj;
        this->owns=owns;
		ref_count  = new size_t;
		*ref_count = 1;
	}
	PW(const PW<X> &other)
	{
		ref_count = other.ref_count;
		obj       = other.obj;
		owns      = other.owns;
		++(*ref_count);
	}
	~PW() { decrement(); }
	PW &operator=(const PW<X> &other)
	{
		decrement();
		ref_count = other.ref_count;
		obj       = other.obj;
		++(*ref_count);
		return *this;
	}
	void reset(const X &obj, bool owns)
	{
		decrement();
		this->obj  = obj;
		this->owns = owns;
		ref_count  = new size_t;
		*ref_count = 1;
	}
	/*
	PW &operator=(const X &other)
	{
	    decrement();
	    ref_count    = new size_t;
	    obj          = other;
	    owns         = false;
	    (*ref_count) = 1;
	    return *this;
	}
	*/
	   operator X() const { return obj; }
	X *operator&() { return &obj; }
};
template <typename X> class PW_explicit : public PW<X>
{
	public:
	PW_explicit() = delete;
	PW_explicit(const PW_explicit<X> &other) : PW<X>(other) {}
	PW_explicit(const PW<X> &other) : PW<X>(other) {}
	operator X() = delete;
};
#endif
