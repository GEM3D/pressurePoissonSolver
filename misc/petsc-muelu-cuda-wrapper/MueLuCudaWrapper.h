#ifndef MUELUCUDAWRAPPER_H
#define MUELUCUDAWRAPPER_H
#include <petscmat.h>
#include <string>
struct Vars;
class MueLuCudaWrapper
{
	private:
	Vars *vars;
	int   num_rows;

	public:
	MueLuCudaWrapper(Mat A, double tol, std::string config);
	~MueLuCudaWrapper();

	void solve(Vec x, Vec b);
    static void initialize();
    static void finalize();
};
#endif
