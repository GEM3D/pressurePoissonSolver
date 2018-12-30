#ifndef UTILS_H
#define UTILS_H
#include "SchurDomain.h"
#include <array>
#include <numeric>
#include <cmath>
namespace Utils
{
inline int index(const int &n, const int &xi, const int &yi, const int &zi)
{
	return xi + yi * n + zi * n * n;
}
inline int index(const int &n, const int &xi, const int &yi)
{
	return xi + yi * n;
}
template <size_t D>
class Slice
{
	private:
	double *start;
    std::array<int,D> strides;

	public:
	inline Slice() {}
	inline Slice(double *start, std::array<int,D> strides)
	{
		this->start   = start;
		this->strides = strides;
	}
	inline double &operator()(std::array<int,D> coord)
	{
		return start[std::inner_product(strides.begin(),strides.end(),coord.begin(),0)];
	}
};
template <size_t D>
inline Slice<D> getSlice(double *u_view, int n, Side<D+1> s)
{
	Slice<D> retval;
    std::array<int,D> strides;
    for(int i=0;i<s.toInt()/2;i++){
        strides[i]=std::pow(n,i);
    }
    for(int i = s.toInt()/2;i<(int)D;i++){
        strides[i]=std::pow(n,i+1);
    }
    if(s.isLowerOnAxis()){
        retval = Slice<D>(&u_view[0],strides);
    }else{
        retval = Slice<D>(&u_view[(n-1)*(int)std::pow(n,s.toInt()/2)],strides);
    }
	return retval;
}
template <size_t D> inline Slice<D-1> getSlice(SchurDomain<D> &d, double *u_view, Side<D> s)
{
	int start = d.local_index * std::pow(d.n,D);
	return getSlice<D-1>(&u_view[start], d.n, s);
}
} // namespace Utils
#endif
