#ifndef UTILS_H
#define UTILS_H
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
class Slice
{
	private:
	double *start;
	int     stridex;
	int     stridey;

	public:
	inline Slice() {}
	inline Slice(double *start, int stridex, int stridey)
	{
		this->start   = start;
		this->stridex = stridex;
		this->stridey = stridey;
	}
	inline double &operator()(int ix, int iy)
	{
		return start[stridex * ix + stridey * iy];
	}
};
inline Slice getSlice(double *u_view, int n, Side<3> s)
{
	Slice retval;
	switch (s.toInt()) {
		case Side<3>::west:
			retval = Slice(&u_view[0], n, n * n);
			break;
		case Side<3>::east:
			retval = Slice(&u_view[(n - 1)], n, n * n);
			break;
		case Side<3>::south:
			retval = Slice(&u_view[0], 1, n * n);
			break;
		case Side<3>::north:
			retval = Slice(&u_view[n * (n - 1)], 1, n * n);
			break;
		case Side<3>::bottom:
			retval = Slice(&u_view[0], 1, n);
			break;
		case Side<3>::top:
			retval = Slice(&u_view[n * n * (n - 1)], 1, n);
			break;
	}
	return retval;
}
template <size_t D> inline Slice getSlice(SchurDomain<D> &d, double *u_view, Side<3> s)
{
	int start = d.local_index * d.n * d.n * d.n;
	return getSlice(&u_view[start], d.n, s);
}
} // namespace Utils
#endif
