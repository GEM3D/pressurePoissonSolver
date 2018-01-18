#ifndef PBMATRIX_H
#define PBMATRIX_H
#include "PW.h"
#include <functional>
#include <memory>
#include <petscmat.h>
#include <set>
#include <map>
#include <valarray>
struct PBlock {
	int                                    i;
	int                                    j;
	std::shared_ptr<std::valarray<double>> coeffs;
	std::function<int(int, int,int)> col_trans;
	std::function<int(int, int,int)> row_trans;
	friend bool operator<(const PBlock &l, const PBlock &r)
	{
		int lc = l.col_trans(4,1, 2);
		int lr = l.row_trans(4,1, 2);
		int rc = r.col_trans(4,1, 2);
		int rr = r.row_trans(4,1, 2);
		return std::tie(l.coeffs, l.i, l.j, lc, lr) < std::tie(r.coeffs, r.i, r.j, rc, rr);
	}
};
class PBMatrix
{
	private:
	int              n;
	std::set<PBlock> blocks;
    typedef   std::function<int(int,int,int)> trans_type;
    typedef  std::map<trans_type,std::valarray<int>,std::function<bool(trans_type,trans_type)>> map_type;
	map_type trans;
	int              local_size  = 0;
	int              global_size = 0;

	public:
	PBMatrix(int n, int local_size, int global_size);
	void insertBlock(int i, int j, std::shared_ptr<std::valarray<double>> coeffs,
	                 std::function<int(int, int,int)> col_trans,
	                 std::function<int(int, int,int)> row_trans);
	void apply(Vec x, Vec b);
	void       finalize();
	static int multiply(Mat A, Vec x, Vec b)
	{
		PBMatrix *pba = nullptr;
		MatShellGetContext(A, &pba);
		pba->apply(x, b);
		return 0;
	}
	PW_explicit<Mat> getMatrix()
	{
		PW<Mat> A;
		MatCreateShell(MPI_COMM_WORLD, local_size, local_size, global_size, global_size, this, &A);
		MatShellSetOperation(A, MATOP_MULT, (void (*)(void)) multiply);
		return A;
	}
};
#endif
