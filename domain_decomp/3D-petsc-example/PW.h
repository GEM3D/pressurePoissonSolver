#ifndef PW_H
#define PW_H
#include <petscsys.h>
template <typename X> class PW
{
	protected:
	X       obj       = nullptr;
	size_t *ref_count = nullptr;
	void    decrement()
	{
		if (ref_count != nullptr) {
			--(*ref_count);
			if (*ref_count == 0) {
				if (obj != nullptr) {
					PetscObjectDestroy((PetscObject *) &obj);
				}
				delete ref_count;
			}
		}
	}

	public:
	PW()
	{
		ref_count  = new size_t;
		*ref_count = 1;
	}
	PW(const PW<X> &other)
	{
		ref_count = other.ref_count;
		obj       = other.obj;
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
	operator X() { return obj; }
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
