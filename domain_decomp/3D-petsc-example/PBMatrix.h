#ifndef PBMATRIX_H
#define PBMATRIX_H
#include "PW.h"
#include <functional>
#include <map>
#include <memory>
#include <petscmat.h>
#include <set>
#include <valarray>
struct PBlock {
	int                                    i;
	int                                    j;
	std::shared_ptr<std::valarray<double>> coeffs;
	std::function<int(int, int, int)> col_trans;
	std::function<int(int, int, int)> row_trans;
	friend bool operator<(const PBlock &l, const PBlock &r)
	{
		char lc = l.col_trans(2, 0, 0);
		char lr = l.row_trans(2, 0, 1);
		char rc = r.col_trans(2, 0, 0);
		char rr = r.row_trans(2, 0, 1);
		return std::tie(l.coeffs, l.i, l.j, lc, lr) < std::tie(r.coeffs, r.i, r.j, rc, rr);
	}
};
class PBMatrix
{
	private:
	int              n;
	std::set<PBlock> blocks;
	typedef std::function<int(int, int, int)> trans_type;
	typedef std::map<trans_type, std::valarray<int>, std::function<bool(trans_type, trans_type)>>
	         map_type;
	map_type trans;
	int      local_size  = 0;
	int      global_size = 0;

	public:
	PBMatrix(int n, int local_size, int global_size);
	void insertBlock(int i, int j, std::shared_ptr<std::valarray<double>> coeffs,
	                 std::function<int(int, int, int)>                    col_trans,
	                 std::function<int(int, int, int)>                    row_trans);
	void apply(Vec x, Vec b);
	void       finalize();
	static int multiply(Mat A, Vec x, Vec b)
	{
		PBMatrix *pba = nullptr;
		MatShellGetContext(A, &pba);
		pba->apply(x, b);
		return 0;
	}
	static int destroy(Mat A)
	{
		PBMatrix *pba = nullptr;
		MatShellGetContext(A, &pba);
		delete pba;
		return 0;
	}
	PW_explicit<Mat> getMatrix()
	{
		PW<Mat> A;
		MatCreateShell(MPI_COMM_WORLD, local_size, local_size, global_size, global_size, this, &A);
		MatShellSetOperation(A, MATOP_MULT, (void (*)(void)) multiply);
		MatShellSetOperation(A, MATOP_DESTROY, (void (*)(void)) destroy);
		return A;
	}
};
#endif
