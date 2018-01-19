#include "PBMatrix.h"
#include <iostream>
using namespace std;
PBMatrix::PBMatrix(int n, int local_size, int global_size)
{
	this->n           = n;
	this->local_size  = local_size;
	this->global_size = global_size;
}
void PBMatrix::insertBlock(int i, int j, std::shared_ptr<std::valarray<double>> coeffs,
                           std::function<int(int, int, int)>                    col_trans,
                           std::function<int(int, int, int)>                    row_trans)
{
	PBlock b = {i, j, coeffs, col_trans, row_trans};
	blocks.insert(b);
}
void PBMatrix::finalize()
{
	trans = map_type([](trans_type l, trans_type r) {
		char l1 = l(2, 0, 0);
		char l2 = l(2, 0, 1);
		char r1 = r(2, 0, 0);
		char r2 = r(2, 0, 1);
		return tie(l1, l2) < tie(r1, r2);
	});

	for (const PBlock &block : blocks) {
		if (trans.count(block.row_trans) == 0) {
			valarray<int> i_block(n * n);
			for (int yi = 0; yi < n; yi++) {
				for (int xi = 0; xi < n; xi++) {
					i_block[block.row_trans(n, xi, yi)] = yi * n + xi;
				}
			}
			trans[block.row_trans] = i_block;
		}
		if (trans.count(block.col_trans) == 0) {
			valarray<int> i_block(n * n);
			for (int yi = 0; yi < n; yi++) {
				for (int xi = 0; xi < n; xi++) {
					i_block[block.col_trans(n, xi, yi)] = yi * n + xi;
				}
			}
			trans[block.col_trans] = i_block;
		}
	}
}
void PBMatrix::apply(Vec x, Vec b)
{
	VecScale(b, 0);
	const double *x_view;
	double *      b_view;
	VecGetArrayRead(x, &x_view);
	VecGetArray(b, &b_view);
	for (const PBlock &block : blocks) {
		int i = block.i * n * n;
		int j = block.j * n * n;

		valarray<double> &coeffs = *block.coeffs;

		int *i_block = &trans[block.row_trans][0];
		int *j_block = &trans[block.col_trans][0];

		for (int li = 0; li < n * n; li++) {
			for (int lj = 0; lj < n * n; lj++) {
				b_view[i + i_block[li]] += coeffs[n * n * li + lj] * x_view[j + j_block[lj]];
			}
		}
	}
	VecRestoreArrayRead(x, &x_view);
	VecRestoreArray(b, &b_view);
}
