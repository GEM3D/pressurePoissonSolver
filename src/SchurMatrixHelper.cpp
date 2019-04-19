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

#include "SchurMatrixHelper.h"
using namespace std;
enum class Rotation : char { x_cw, x_ccw, y_cw, y_ccw, z_cw, z_ccw };
struct Block {
	static const Side<3>          side_table[6][6];
	static const char             rots_table[6][6];
	static const vector<Rotation> main_rot_plan[6];
	static const vector<Rotation> aux_rot_plan_dirichlet[6];
	static const vector<Rotation> aux_rot_plan_neumann[16];
	static const char             rot_quad_lookup_left[4][4];
	static const char             rot_quad_lookup_right[4][4];
	static const char             quad_flip_lookup[4];
	char                          data = 0;
	IfaceType                     type;
	Side<3>                       main;
	Side<3>                       aux;
	int                           j;
	int                           i;
	bitset<6>                     neumann;
	Block(Side<3> main, int j, Side<3> aux, int i, bitset<6> neumann, IfaceType type)
	{
		this->main    = main;
		this->j       = j;
		this->aux     = aux;
		this->i       = i;
		this->neumann = neumann;
		this->type    = type;
		data          = 0;
		data |= main.isLowerOnAxis() << 7;
		data |= aux.isLowerOnAxis() << 3;
		rotate();
	}
	void applyRotation(const Rotation rot)
	{
		char r;
		// main rotation
		r = (data & ~(~0u << 2) << 4) >> 4;
		r = (r + rots_table[static_cast<int>(rot)][main.toInt()]) & 0b11;
		data &= ~(~(~0u << 2) << 4);
		data |= r << 4;
		// aux rotation
		r = data & ~(~0u << 2);
		r = (r + rots_table[static_cast<int>(rot)][aux.toInt()]) & 0b11;
		data &= ~0u << 2;
		data |= r;
		main                  = side_table[static_cast<int>(rot)][main.toInt()];
		aux                   = side_table[static_cast<int>(rot)][aux.toInt()];
		bitset<6> old_neumann = neumann;
		for (int i = 0; i < 6; i++) {
			neumann[side_table[(int) rot][i].toInt()] = old_neumann[i];
		}
	}
	void rotate()
	{
		for (Rotation rot : main_rot_plan[main.toInt()]) {
			applyRotation(rot);
		}
		if (neumann.to_ulong() == 0) {
			for (Rotation rot : aux_rot_plan_dirichlet[aux.toInt()]) {
				applyRotation(rot);
			}
		} else {
			for (Rotation rot : aux_rot_plan_neumann[neumann.to_ulong() >> 2]) {
				applyRotation(rot);
			}
		}
		// updated iface type
		auto rotateQuad = [&](int quad) {
			if (auxOrigLeft()) {
				quad = rot_quad_lookup_left[auxRot()][quad];
			} else {
				quad = rot_quad_lookup_right[auxRot()][quad];
			}
			if (auxFlipped()) { quad = quad_flip_lookup[quad]; }
			return quad;
		};
		switch (type.toInt()) {
			case IfaceType::fine_to_coarse:
			case IfaceType::fine_to_fine:
			case IfaceType::coarse_to_fine: {
				int quad = type.getOrthant();
				quad     = rotateQuad(quad);
				type.setOrthant(quad);
			} break;
			default:
				break;
		}
	}
	bool operator<(const Block &b) const
	{
		return std::tie(i, j, data) < std::tie(b.i, b.j, b.data);
	}
	bool operator==(const Block &b) const
	{
		return neumann.to_ulong() == b.neumann.to_ulong();
	}
	//
	bool mainLeft()
	{
		return main.isLowerOnAxis();
	}
	bool mainFlipped()
	{
		bool left      = main.isLowerOnAxis();
		bool orig_left = (data >> 7) & 0b1;
		return left != orig_left;
	}
	int mainRot()
	{
		return (data >> 4) & 0b11;
	}
	bool auxOrigLeft()
	{
		return (data >> 3) & 0b1;
	}
	bool auxLeft()
	{
		return aux.isLowerOnAxis();
	}
	bool auxFlipped()
	{
		bool left      = aux.isLowerOnAxis();
		bool orig_left = (data >> 3) & 0b1;
		return left != orig_left;
	}
	int auxRot()
	{
		return data & 0b11;
	}
};
struct BlockKey {
	IfaceType type;
	Side<3>   s;

	BlockKey() {}
	BlockKey(const Block &b)
	{
		type = b.type;
		s    = b.aux;
	}
	friend bool operator<(const BlockKey &l, const BlockKey &r)
	{
		return tie(l.s, l.type) < tie(r.s, r.type);
	}
};

const Side<3> Block::side_table[6][6]
= {{Side<3>::west, Side<3>::east, Side<3>::top, Side<3>::bottom, Side<3>::south, Side<3>::north},
   {Side<3>::west, Side<3>::east, Side<3>::bottom, Side<3>::top, Side<3>::north, Side<3>::south},
   {Side<3>::bottom, Side<3>::top, Side<3>::south, Side<3>::north, Side<3>::east, Side<3>::west},
   {Side<3>::top, Side<3>::bottom, Side<3>::south, Side<3>::north, Side<3>::west, Side<3>::east},
   {Side<3>::north, Side<3>::south, Side<3>::west, Side<3>::east, Side<3>::bottom, Side<3>::top},
   {Side<3>::south, Side<3>::north, Side<3>::east, Side<3>::west, Side<3>::bottom, Side<3>::top}};
const char Block::rots_table[6][6] = {{3, 1, 0, 0, 2, 2}, {1, 3, 2, 2, 0, 0}, {1, 3, 3, 1, 1, 3},
                                      {1, 3, 1, 3, 3, 1}, {0, 0, 0, 0, 3, 1}, {0, 0, 0, 0, 1, 3}};
const vector<Rotation> Block::main_rot_plan[6] = {{},
                                                  {Rotation::z_cw, Rotation::z_cw},
                                                  {Rotation::z_cw},
                                                  {Rotation::z_ccw},
                                                  {Rotation::y_ccw},
                                                  {Rotation::y_cw}};
const vector<Rotation> Block::aux_rot_plan_dirichlet[6]
= {{}, {}, {}, {Rotation::x_cw, Rotation::x_cw}, {Rotation::x_cw}, {Rotation::x_ccw}};
const vector<Rotation> Block::aux_rot_plan_neumann[16] = {{},
                                                          {},
                                                          {Rotation::x_cw, Rotation::x_cw},
                                                          {},
                                                          {Rotation::x_cw},
                                                          {Rotation::x_cw},
                                                          {Rotation::x_cw, Rotation::x_cw},
                                                          {Rotation::x_cw, Rotation::x_cw},
                                                          {Rotation::x_ccw},
                                                          {},
                                                          {Rotation::x_ccw},
                                                          {},
                                                          {Rotation::x_ccw},
                                                          {Rotation::x_cw},
                                                          {Rotation::x_ccw},
                                                          {}};
const char             Block::rot_quad_lookup_left[4][4]
= {{0, 1, 2, 3}, {1, 3, 0, 2}, {3, 2, 1, 0}, {2, 0, 3, 1}};
const char Block::rot_quad_lookup_right[4][4]
= {{0, 1, 2, 3}, {2, 0, 3, 1}, {3, 2, 1, 0}, {1, 3, 0, 2}};
const char Block::quad_flip_lookup[4] = {1, 0, 3, 2};

void SchurMatrixHelper::assembleMatrix(inserter insertBlock)
{
	set<Block> blocks;
	for (auto &p : sh->getIfaces()) {
		const IfaceSet<3> &ifs = p.second;
		int                i   = ifs.id_global;
		for (const Iface<3> &iface : ifs.ifaces) {
			Side<3> aux = iface.s;
			for (int s = 0; s < 6; s++) {
				int j = iface.global_id[s];
				if (j != -1) {
					Side<3> main = static_cast<Side<3>>(s);
					blocks.insert(Block(main, j, aux, i, iface.neumann, iface.type));
				}
			}
		}
	}
	Vec u, f, r, e, gamma, interp;
	VecCreateSeq(PETSC_COMM_SELF, n * n * n, &u);
	VecCreateSeq(PETSC_COMM_SELF, n * n * n, &f);
	VecCreateSeq(PETSC_COMM_SELF, n * n * n, &r);
	VecCreateSeq(PETSC_COMM_SELF, n * n * n, &e);
	VecCreateSeq(PETSC_COMM_SELF, n * n, &gamma);
	VecCreateSeq(PETSC_COMM_SELF, n * n, &interp);
	double *interp_view, *gamma_view;
	VecGetArray(interp, &interp_view);
	VecGetArray(gamma, &gamma_view);
	while (!blocks.empty()) {
		// the first in the set is the type of interface that we are going to solve for
		set<Block> todo;
		Block      curr_type = *blocks.begin();
		blocks.erase(blocks.begin());
		todo.insert(curr_type);
		set<Block> to_be_deleted;
		for (auto iter = blocks.begin(); iter != blocks.end(); iter++) {
			if (*iter == curr_type) {
				todo.insert(*iter);
				to_be_deleted.insert(*iter);
			}
		}
		for (Block i : to_be_deleted) {
			blocks.erase(i);
		}

		auto solver       = sh->getSolver();
		auto interpolator = sh->getInterpolator();
		// create domain representing curr_type
		SchurDomain<3> sd;
		sd.n = n;
		sd.domain.lengths.fill(1);
		sd.neumann                        = curr_type.neumann;
		sd.getIfaceInfoPtr(Side<3>::west) = new NormalIfaceInfo<3>();
		solver->addDomain(sd);
		std::deque<SchurDomain<3>> single_domain;
		single_domain.push_back(sd);

		map<BlockKey, shared_ptr<valarray<double>>> coeffs;
		// allocate blocks of coefficients
		for (const Block &b : todo) {
			shared_ptr<valarray<double>> ptr = coeffs[b];
			if (ptr.get() == nullptr) {
				coeffs[b] = shared_ptr<valarray<double>>(new valarray<double>(n * n * n * n));
			}
		}

		for (int j = 0; j < n * n; j++) {
			gamma_view[j] = 1;
			solver->domainSolve(single_domain, f, u, gamma);
			gamma_view[j] = 0;

			// fill the blocks
			for (auto &p : coeffs) {
				BlockKey  bk   = p.first;
				Side<3>   s    = bk.s;
				IfaceType type = bk.type;
				VecScale(interp, 0);
	            std::shared_ptr<Vector<3>> u_vec(new PetscVector<3>(u, n));
	            std::shared_ptr<Vector<2>> interp_vec(new PetscVector<2>(interp, n));
				interpolator->interpolate(sd, s, 0, type, u_vec, interp_vec);
				valarray<double> &block = *p.second;
				for (int i = 0; i < n * n; i++) {
					block[i * n * n + j] = -interp_view[i];
				}
				if (s == Side<3>::west) {
					switch (type.toInt()) {
						case IfaceType::normal:
							block[n * n * j + j] += 0.5;
							break;
						case IfaceType::coarse_to_coarse:
                        case IfaceType::fine_to_fine:
							block[n * n * j + j] += 1.0;
							break;
						default:
							break;
					}
				}
			}
		}

		// now insert these results into the matrix for each interface
		for (Block block : todo) {
			insertBlock(&block, coeffs[block]);
		}
	}
	VecRestoreArray(interp, &interp_view);
	VecRestoreArray(gamma, &gamma_view);
	VecDestroy(&u);
	VecDestroy(&f);
	VecDestroy(&r);
	VecDestroy(&e);
	VecDestroy(&gamma);
	VecDestroy(&interp);
}
PW_explicit<Mat> SchurMatrixHelper::formCRSMatrix()
{
	PW<Mat> A;
	MatCreate(MPI_COMM_WORLD, &A);
	int local_size  = sh->getIfaces().size() * n * n;
	int global_size = sh->getSchurVecGlobalSize();
	MatSetSizes(A, local_size, local_size, global_size, global_size);
	MatSetType(A, MATMPIAIJ);
	MatMPIAIJSetPreallocation(A, 10 * n * n, nullptr, 10 * n * n, nullptr);

	auto insertBlock = [&](Block *b, shared_ptr<valarray<double>> coeffs) {
		int global_i = b->i * n * n;
		int global_j = b->j * n * n;

		valarray<double> &            orig = *coeffs;
		const function<int(int, int)> transforms_left[4]
		= {[&](int xi, int yi) { return xi + yi * n; },
		   [&](int xi, int yi) { return n - yi - 1 + xi * n; },
		   [&](int xi, int yi) { return n - xi - 1 + (n - yi - 1) * n; },
		   [&](int xi, int yi) { return yi + (n - xi - 1) * n; }};

		const function<int(int, int)> transforms_right[4]
		= {[&](int xi, int yi) { return xi + yi * n; },
		   [&](int xi, int yi) { return yi + (n - xi - 1) * n; },
		   [&](int xi, int yi) { return n - xi - 1 + (n - yi - 1) * n; },
		   [&](int xi, int yi) { return n - yi - 1 + xi * n; }};

		const function<int(int, int)> transforms_left_inv[4]
		= {[&](int xi, int yi) {
			   xi = n - xi - 1;
			   return xi + yi * n;
		   },
		   [&](int xi, int yi) {
			   xi = n - xi - 1;
			   return n - yi - 1 + xi * n;
		   },
		   [&](int xi, int yi) {
			   xi = n - xi - 1;
			   return n - xi - 1 + (n - yi - 1) * n;
		   },
		   [&](int xi, int yi) {
			   xi = n - xi - 1;
			   return yi + (n - xi - 1) * n;
		   }};
		const function<int(int, int)> transforms_right_inv[4]
		= {[&](int xi, int yi) {
			   xi = n - xi - 1;
			   return xi + yi * n;
		   },
		   [&](int xi, int yi) {
			   xi = n - xi - 1;
			   return yi + (n - xi - 1) * n;
		   },
		   [&](int xi, int yi) {
			   xi = n - xi - 1;
			   return n - xi - 1 + (n - yi - 1) * n;
		   },
		   [&](int xi, int yi) {
			   xi = n - xi - 1;
			   return n - yi - 1 + xi * n;
		   }};

		valarray<double>        copy(n * n * n * n);
		function<int(int, int)> col_trans, row_trans;
		if (b->mainLeft()) {
			if (b->mainFlipped()) {
				col_trans = transforms_left_inv[b->mainRot()];
			} else {
				col_trans = transforms_left[b->mainRot()];
			}
		} else {
			if (b->mainFlipped()) {
				col_trans = transforms_right_inv[b->mainRot()];
			} else {
				col_trans = transforms_right[b->mainRot()];
			}
		}
		if (b->auxLeft()) {
			if (b->auxFlipped()) {
				row_trans = transforms_left_inv[b->auxRot()];
			} else {
				row_trans = transforms_left[b->auxRot()];
			}
		} else {
			if (b->auxFlipped()) {
				row_trans = transforms_right_inv[b->auxRot()];
			} else {
				row_trans = transforms_right[b->auxRot()];
			}
		}
		for (int row_yi = 0; row_yi < n; row_yi++) {
			for (int row_xi = 0; row_xi < n; row_xi++) {
				int i_dest = row_xi + row_yi * n;
				int i_orig = row_trans(row_xi, row_yi);
				for (int col_yi = 0; col_yi < n; col_yi++) {
					for (int col_xi = 0; col_xi < n; col_xi++) {
						int j_dest = col_xi + col_yi * n;
						int j_orig = col_trans(col_xi, col_yi);

						copy[i_dest * n * n + j_dest] = orig[i_orig * n * n + j_orig];
					}
				}
			}
		}
		vector<int> inds_i(n * n);
		iota(inds_i.begin(), inds_i.end(), global_i);
		vector<int> inds_j(n * n);
		iota(inds_j.begin(), inds_j.end(), global_j);

		MatSetValues(A, n * n, &inds_i[0], n * n, &inds_j[0], &copy[0], ADD_VALUES);
	};

	assembleMatrix(insertBlock);
	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
	return A;
}
PBMatrix *SchurMatrixHelper::formPBMatrix()
{
	int local_size  = sh->getIfaces().size() * n * n;
	int global_size = sh->getSchurVecGlobalSize();

	PBMatrix *APB         = new PBMatrix(n, local_size, global_size);
	auto      insertBlock = [&](Block *b, shared_ptr<valarray<double>> coeffs) {
        int global_i = b->i;
        int global_j = b->j;

        const function<int(int, int, int)> transforms_left[4]
        = {[](int n, int xi, int yi) { return xi + yi * n; },
           [](int n, int xi, int yi) { return n - yi - 1 + xi * n; },
           [](int n, int xi, int yi) { return n - xi - 1 + (n - yi - 1) * n; },
           [](int n, int xi, int yi) { return yi + (n - xi - 1) * n; }};

        const function<int(int, int, int)> transforms_right[4]
        = {[](int n, int xi, int yi) { return xi + yi * n; },
           [](int n, int xi, int yi) { return yi + (n - xi - 1) * n; },
           [](int n, int xi, int yi) { return n - xi - 1 + (n - yi - 1) * n; },
           [](int n, int xi, int yi) { return n - yi - 1 + xi * n; }};

        const function<int(int, int, int)> transforms_left_inv[4]
        = {[](int n, int xi, int yi) {
               xi = n - xi - 1;
               return xi + yi * n;
           },
           [](int n, int xi, int yi) {
               xi = n - xi - 1;
               return n - yi - 1 + xi * n;
           },
           [](int n, int xi, int yi) {
               xi = n - xi - 1;
               return n - xi - 1 + (n - yi - 1) * n;
           },
           [](int n, int xi, int yi) {
               xi = n - xi - 1;
               return yi + (n - xi - 1) * n;
           }};
        const function<int(int, int, int)> transforms_right_inv[4]
        = {[](int n, int xi, int yi) {
               xi = n - xi - 1;
               return xi + yi * n;
           },
           [](int n, int xi, int yi) {
               xi = n - xi - 1;
               return yi + (n - xi - 1) * n;
           },
           [](int n, int xi, int yi) {
               xi = n - xi - 1;
               return n - xi - 1 + (n - yi - 1) * n;
           },
           [](int n, int xi, int yi) {
               xi = n - xi - 1;
               return n - yi - 1 + xi * n;
           }};

        function<int(int, int, int)> col_trans, row_trans;
        if (b->mainLeft()) {
            if (b->mainFlipped()) {
                col_trans = transforms_left_inv[b->mainRot()];
            } else {
                col_trans = transforms_left[b->mainRot()];
            }
        } else {
            if (b->mainFlipped()) {
                col_trans = transforms_right_inv[b->mainRot()];
            } else {
                col_trans = transforms_right[b->mainRot()];
            }
        }
        if (b->auxLeft()) {
            if (b->auxFlipped()) {
                row_trans = transforms_left_inv[b->auxRot()];
            } else {
                row_trans = transforms_left[b->auxRot()];
            }
        } else {
            if (b->auxFlipped()) {
                row_trans = transforms_right_inv[b->auxRot()];
            } else {
                row_trans = transforms_right[b->auxRot()];
            }
        }
        APB->insertBlock(global_i, global_j, coeffs, col_trans, row_trans);
	};

	assembleMatrix(insertBlock);
	APB->finalize();

	return APB;
}
PW_explicit<Mat> SchurMatrixHelper::getPBMatrix()
{
	return formPBMatrix()->getMatrix();
}
PW_explicit<Mat> SchurMatrixHelper::getPBDiagInv()
{
	return formPBMatrix()->getDiagInv()->getMatrix();
}
void SchurMatrixHelper::getPBDiagInv(PC p)
{
	formPBMatrix()->getDiagInv()->getPrec(p);
}
