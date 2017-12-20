#ifndef MUELUWRAPPER_H
#define MUELUWRAPPER_H
#include <petscmat.h>
#include <string>
struct Vars;
class MueLuWrapper
{
	private:
	Vars *vars;
	int   num_rows;

	public:
	MueLuWrapper(Mat A, double tol, std::string config);
	~MueLuWrapper();

	void solve(Vec x, Vec b);
};
#endif
