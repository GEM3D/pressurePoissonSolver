#ifndef DDMultiGrid_H
#define DDMultiGrid_H
#include "MyTypeDefs.h"
class Level
{
	private:
	Teuchos::RCP<RBMatrix>           A;
	Teuchos::RCP<Tpetra::Operator<>> op;
	Teuchos::RCP<matrix_type>        cA;
	Teuchos::RCP<BlockJacobiRelaxer> R;
	Teuchos::RCP<Amesos2::Solver<matrix_type, vector_type>> solver;

	Teuchos::RCP<vector_type> x;
	Teuchos::RCP<vector_type> b;
	Teuchos::RCP<vector_type> r;

	Teuchos::RCP<single_vector_type> s;

	bool peak = false;

	int    n     = -1;
	Level *child = nullptr;

	public:
	Level() = default;
	Level(DomainCollection *dc, int n, bool neumann = false)
	{
		this->n                    = n;
		Teuchos::RCP<map_type> map = dc->formMatrixMap(n);
		x                          = rcp(new vector_type(map, 1));
		b                          = rcp(new vector_type(map, 1));
		r                          = rcp(new vector_type(map, 1));
		if (n / 2 > 4) {
			child = new Level(dc, n / 2, neumann);
        }
        if(n / 2 > 4){
			dc->formRBMatrix(map, A, &s, n);
			if (neumann) {
				R                        = Teuchos::rcp(new BlockJacobiRelaxer(A, s));
				Teuchos::RCP<OpShift> os = Teuchos::rcp(new OpShift(A, s));
				this->op                 = os;
			} else {
				op = A;
				R  = Teuchos::rcp(new BlockJacobiRelaxer(A));
			}
		} else {
			dc->formCRSMatrix(map, cA, &s, n, true);
			std::cerr << "Hello" << std::endl;
			if (neumann) {
				cA->resumeFill();
				size_t                     size = cA->getNumEntriesInGlobalRow(0);
				Teuchos::ArrayView<int>    inds(new int[size], size);
				Teuchos::ArrayView<double> vals(new double[size], size);
				cA->getGlobalRowCopy(0, inds, vals, size);
				for (size_t i = 0; i < size; i++) {
					if (inds[i] == 0) {
						vals[i] = 1;
					} else {
						vals[i] = 0;
					}
				}
				cA->replaceGlobalValues(0, inds, vals);
				cA->fillComplete();
			}

			solver = Amesos2::create<matrix_type, vector_type>("KLU2", cA, x, b);

			Teuchos::RCP<Teuchos::ParameterList> amesoslist
			= Teuchos::rcp(new Teuchos::ParameterList);
			amesoslist->set("Trans", "TRANS", "Solve with transpose");
			amesoslist->set("Transpose", true);
			solver->setParameters(amesoslist);

			solver->symbolicFactorization().numericFactorization();
		}
	}
	~Level() { delete child; }
	void cycle(const vector_type &r_fine, vector_type &e_fine)
	{
		// coarsen vectors
		auto r_fine_view = r_fine.getLocalView<Kokkos::HostSpace>();
		auto b_view      = b->getLocalView<Kokkos::HostSpace>();
		b->putScalar(0);
		for (size_t i = 0; i < r_fine.getLocalLength(); i++) {
			b_view(i / 2, 0) += r_fine_view(i, 0);
		}
		b->scale(0.5);

		if (child != nullptr) {
			// relax
			x->putScalar(0);
			R->apply(*b, *x);

			// residual
			op->apply(*x, *r);
			r->update(1.0, *b, -1.0);

			// cycle
			child->cycle(*r, *x);

			// relax
			R->apply(*b, *x, Teuchos::NO_TRANS, 1.0, 1.0);
			if (false) {
				peak = true;
				// residual
				A->apply(*x, *r);
				r->update(1.0, *b, -1.0);

				// cycle
				child->cycle(*r, *x);

				// relax
				R->apply(*b, *x, Teuchos::NO_TRANS, 1.0, 1.0);
			}

		} else {
            if(solver.is_null()){
			// relax
			x->putScalar(0);
			R->apply(*b, *x);
				// relax
				R->apply(*b, *x, Teuchos::NO_TRANS, 1.0, 1.0);
            }else{
			solver->solve();
            }
		}
		// interpolate
		auto e_view = e_fine.getLocalView<Kokkos::HostSpace>();
		auto x_view = x->getLocalView<Kokkos::HostSpace>();
		for (size_t i = 0; i < r_fine.getLocalLength(); i++) {
			int iface_i = i % (n * 2);
			if (iface_i == 0) {
				e_view(i, 0) += 2 * x_view(i / 2, 0) + x_view(i / 2 + 1, 0);
			} else if (iface_i == 2 * n - 1) {
				e_view(i, 0) += 2 * x_view(i / 2, 0) + x_view(i / 2 - 1, 0);
			} else {
				if (iface_i % 2 == 1) {
					e_view(i, 0) += 3.0 / 4.0 * x_view(i / 2, 0) + 1.0 / 4.0 * x_view(i / 2 + 1, 0);
				} else {
					e_view(i, 0) += 1.0 / 4.0 * x_view(i / 2 - 1, 0) + 3.0 / 4.0 * x_view(i / 2, 0);
				}
			}
		}
	}
};
class DDMultiGrid : public Tpetra::Operator<scalar_type>
{
	private:
	DomainCollection *dc;
	Level *           bottom;
	int               max_coarse_size = 4;

	Teuchos::RCP<const map_type>     map;
	Teuchos::RCP<BlockJacobiRelaxer> R;
	Teuchos::RCP<vector_type>        r;
	Teuchos::RCP<Tpetra::Operator<>> A;

	public:
	DDMultiGrid(DomainCollection *dc, Teuchos::RCP<RBMatrix> A,
	            Teuchos::RCP<single_vector_type> s = Teuchos::null)
	{
		this->dc = dc;
		if (s.is_null()) {
			this->A = A;
			R       = Teuchos::rcp(new BlockJacobiRelaxer(A));
		} else {
			R                        = Teuchos::rcp(new BlockJacobiRelaxer(A, s));
			Teuchos::RCP<OpShift> os = Teuchos::rcp(new OpShift(A, s));
			this->A                  = os;
		}
		map    = dc->formMatrixMap(dc->n);
		bottom = new Level(dc, dc->n / 2, !s.is_null());
		r      = rcp(new vector_type(map, 1));
		R      = Teuchos::rcp(new BlockJacobiRelaxer(A));
	}
	~DDMultiGrid() { delete bottom; }
	void apply(const vector_type &b, vector_type &x, Teuchos::ETransp mode = Teuchos::NO_TRANS,
	           scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
	           scalar_type beta  = Teuchos::ScalarTraits<scalar_type>::zero()) const
	{
		// relax
		x.putScalar(0);
		R->apply(b, x);

		for (int i = 0; i < 1; i++) {
			// residual
			A->apply(x, *r);
			r->update(1.0, b, -1.0);

			// cycle
			bottom->cycle(*r, x);

			// relax
			R->apply(b, x, Teuchos::NO_TRANS, 1.0, 1.0);
		}
	}
	Teuchos::RCP<const map_type> getDomainMap() const { return map; }
	Teuchos::RCP<const map_type> getRangeMap() const { return map; }
};
#endif
