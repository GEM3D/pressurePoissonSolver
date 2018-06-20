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
