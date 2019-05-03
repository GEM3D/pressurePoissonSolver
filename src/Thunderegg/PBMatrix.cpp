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

#include "PBMatrix.h"
#include <iostream>
extern "C" {
// LU decomoposition of a general matrix
void dgetrf_(int *M, int *N, double *A, int *lda, int *IPIV, int *INFO);

// generate inverse of a matrix given its LU decomposition
void dgetri_(int *N, double *A, int *lda, int *IPIV, double *WORK, int *lwork, int *INFO);
}
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

	std::set<shared_ptr<valarray<double>>> bs;
	for (const PBlock &block : blocks) {
		bs.insert(block.coeffs);
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
	cerr << "NUM_BLOCKS: " << bs.size() << endl;
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
void inverse(double *A, int N)
{
	int *   IPIV  = new int[N + 1];
	int     LWORK = N * N;
	double *WORK  = new double[LWORK];
	int     INFO;

	dgetrf_(&N, &N, A, &N, IPIV, &INFO);
	dgetri_(&N, A, &N, IPIV, WORK, &LWORK, &INFO);

	delete[] IPIV;
	delete[] WORK;
}
struct PBlockMin {
	std::shared_ptr<std::valarray<double>> coeffs;
	std::function<int(int, int, int)> col_trans;
	std::function<int(int, int, int)> row_trans;
	PBlockMin(const PBlock &b)
	{
		coeffs    = b.coeffs;
		col_trans = b.col_trans;
		row_trans = b.row_trans;
	}
	friend bool operator<(const PBlockMin &l, const PBlockMin &r)
	{
		char lc = l.col_trans(2, 0, 0);
		char lr = l.row_trans(2, 0, 1);
		char rc = r.col_trans(2, 0, 0);
		char rr = r.row_trans(2, 0, 1);
		return std::tie(l.coeffs, lc, lr) < std::tie(r.coeffs, rc, rr);
	}
};
PBMatrix *PBMatrix::getDiagInv()
{
	PBMatrix *APB = new PBMatrix(n, local_size, global_size);
	map<int, set<PBlockMin>> diags;
	for (const PBlock &b : blocks) {
		if (b.i == b.j) {
			diags[b.i].insert(b);
		}
	}
	map<set<PBlockMin>, shared_ptr<valarray<double>>> invs;
	for (auto &p : diags) {
		int                           i   = p.first;
		shared_ptr<valarray<double>> &inv = invs[p.second];
		if (inv == nullptr) {
			inv.reset(new valarray<double>(n * n * n * n));
			valarray<double> &diag_coeffs = *inv;

			// sum up blocks
			for (const PBlockMin &b : p.second) {
				valarray<double> &coeffs = *b.coeffs;
				for (int row_yi = 0; row_yi < n; row_yi++) {
					for (int row_xi = 0; row_xi < n; row_xi++) {
						int i_dest = row_xi + row_yi * n;
						int i_orig = b.row_trans(n, row_xi, row_yi);
						for (int col_yi = 0; col_yi < n; col_yi++) {
							for (int col_xi = 0; col_xi < n; col_xi++) {
								int j_dest = col_xi + col_yi * n;
								int j_orig = b.col_trans(n, col_xi, col_yi);

								diag_coeffs[i_dest * n * n + j_dest]
								+= coeffs[i_orig * n * n + j_orig];
							}
						}
					}
				}
			}

			// invert
			inverse(&diag_coeffs[0], n * n);
		}

		auto trans = [](int n, int xi, int yi) { return xi + yi * n; };
		APB->insertBlock(i, i, inv, trans, trans);
	}
	APB->finalize();
	return APB;
}
BlockJacobiSmoother PBMatrix::getBlockJacobiSmoother()
{
	PBMatrix *APB = new PBMatrix(n, local_size, global_size);
	PBMatrix *R = new PBMatrix(n, local_size, global_size);
	map<int, set<PBlockMin>> diags;
	for (const PBlock &b : blocks) {
		if (b.i == b.j) {
			diags[b.i].insert(b);
		}else{
			R->blocks.insert(b);
		}
	}
	map<set<PBlockMin>, shared_ptr<valarray<double>>> invs;
	for (auto &p : diags) {
		int                           i   = p.first;
		shared_ptr<valarray<double>> &inv = invs[p.second];
		if (inv == nullptr) {
			inv.reset(new valarray<double>(n * n * n * n));
			valarray<double> &diag_coeffs = *inv;

			// sum up blocks
			for (const PBlockMin &b : p.second) {
				valarray<double> &coeffs = *b.coeffs;
				for (int row_yi = 0; row_yi < n; row_yi++) {
					for (int row_xi = 0; row_xi < n; row_xi++) {
						int i_dest = row_xi + row_yi * n;
						int i_orig = b.row_trans(n, row_xi, row_yi);
						for (int col_yi = 0; col_yi < n; col_yi++) {
							for (int col_xi = 0; col_xi < n; col_xi++) {
								int j_dest = col_xi + col_yi * n;
								int j_orig = b.col_trans(n, col_xi, col_yi);

								diag_coeffs[i_dest * n * n + j_dest]
								+= coeffs[i_orig * n * n + j_orig];
							}
						}
					}
				}
			}

			// invert
			inverse(&diag_coeffs[0], n * n);
		}

		auto trans = [](int n, int xi, int yi) { return xi + yi * n; };
		APB->insertBlock(i, i, inv, trans, trans);
	}
	APB->finalize();
	R->finalize();
	BlockJacobiSmoother jacobi;
	jacobi.D.reset(APB);
	jacobi.R.reset(R);
	return jacobi;
}
