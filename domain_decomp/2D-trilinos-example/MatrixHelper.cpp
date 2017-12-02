#include "MatrixHelper.h"
#include "StencilHelper.h"
#include <vector>
using namespace std;
MatrixHelper::MatrixHelper(DomainCollection dc, Teuchos::RCP<const Teuchos::Comm<int>> comm)
{
	this->dc   = dc;
	this->comm = comm;
}
Teuchos::RCP<matrix_type> MatrixHelper::formCRSMatrix()
{
	Teuchos::RCP<matrix_type> A = Teuchos::rcp(new matrix_type(dc.getDomainRowMap(), 10));
	for (auto &p : dc.domains) {
		Domain &d     = p.second;
		int     n     = d.n;
		double  h_x   = d.x_length / n;
		double  h_y   = d.y_length / n;
		int     start = n * n * d.id;

		// center coeffs
		double coeff = -2.0 / (h_x * h_x) - 2.0 / (h_y * h_y);
		for (int y_i = 0; y_i < n; y_i++) {
			for (int x_i = 0; x_i < n; x_i++) {
				int row = start + x_i + n * y_i;
				A->insertGlobalValues(row, 1, &coeff, &row);
			}
		}
		// north coeffs
		coeff = 1.0 / (h_y * h_y);
		for (int y_i = 0; y_i < n - 1; y_i++) {
			for (int x_i = 0; x_i < n; x_i++) {
				int row = start + x_i + n * y_i;
				int col = start + x_i + n * (y_i + 1);
				A->insertGlobalValues(row, 1, &coeff, &col);
			}
		}
		// south coeffs
		for (int y_i = 1; y_i < n; y_i++) {
			for (int x_i = 0; x_i < n; x_i++) {
				int row = start + x_i + n * y_i;
				int col = start + x_i + n * (y_i - 1);
				A->insertGlobalValues(row, 1, &coeff, &col);
			}
		}
		coeff = 1.0 / (h_x * h_x);
		// east coeffs
		for (int y_i = 0; y_i < n; y_i++) {
			for (int x_i = 0; x_i < n - 1; x_i++) {
				int row = start + x_i + n * y_i;
				int col = start + x_i + 1 + n * y_i;
				A->insertGlobalValues(row, 1, &coeff, &col);
			}
		}
		// west coeffs
		for (int y_i = 0; y_i < n; y_i++) {
			for (int x_i = 1; x_i < n; x_i++) {
				int row = start + x_i + n * y_i;
				int col = start + x_i - 1 + n * y_i;
				A->insertGlobalValues(row, 1, &coeff, &col);
			}
		}
		// boundaries
		Side s = Side::north;
		do {
			StencilHelper *sh = getStencilHelper(d, s);
			for (int i = 0; i < n; i++) {
				int     row    = sh->row(i);
				int     size   = sh->size(i);
				double *coeffs = sh->coeffs(i);
				int *   cols   = sh->cols(i);
				A->insertGlobalValues(row, size, coeffs, cols);
			}
			delete sh;
			s++;
		} while (s != Side::north);
	}
	A->fillComplete();
	return A;
}
