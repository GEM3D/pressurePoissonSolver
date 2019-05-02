#include <Operators/Operator.h>
#include <Vector.h>
template <size_t D> class BiCGStab
{
	public:
	static int solve(std::shared_ptr<VectorGenerator<D>> vg, std::shared_ptr<const Operator<D>> A,
	                 std::shared_ptr<Vector<D>> x, std::shared_ptr<const Vector<D>> b,
	                 std::shared_ptr<const Operator<D>> Mr = nullptr)
	{
		std::shared_ptr<Vector<D>> resid = vg->getNewVector();
		std::shared_ptr<Vector<D>> ms;
		std::shared_ptr<Vector<D>> mp;
		if (Mr != nullptr) {
			ms = vg->getNewVector();
			mp = vg->getNewVector();
		}
		A->apply(x, resid);
		resid->scaleThenAdd(-1, b);
		double                     r0_norm = resid->twoNorm();
		std::shared_ptr<Vector<D>> rhat    = vg->getNewVector();
		rhat->copy(resid);
		std::shared_ptr<Vector<D>> p = vg->getNewVector();
		p->copy(resid);
		std::shared_ptr<Vector<D>> ap = vg->getNewVector();
		std::shared_ptr<Vector<D>> as = vg->getNewVector();

		std::shared_ptr<Vector<D>> s   = vg->getNewVector();
		double                     rho = rhat->dot(resid);

		double tolerance = 1e-12;
		int    num_its   = 0;
		int    max_it    = 1000;
		while (resid->twoNorm() / r0_norm > tolerance && num_its < max_it) {
			if (Mr != nullptr) {
				Mr->apply(p, mp);
				A->apply(mp, ap);
			} else {
				A->apply(p, ap);
			}
			double alpha = rho / rhat->dot(ap);
			s->copy(resid);
			s->addScaled(-alpha, ap);
			if (Mr != nullptr) {
				Mr->apply(s, ms);
				A->apply(ms, as);
			} else {
				A->apply(s, as);
			}
			double omega = as->dot(s) / as->dot(as);

			// update x and residual
			if (Mr != nullptr) {
				x->addScaled(alpha, mp, omega, ms);
			} else {
				x->addScaled(alpha, p, omega, s);
			}
			resid->addScaled(-alpha, ap, -omega, as);

			double rho_new = resid->dot(rhat);
			double beta    = rho_new * alpha / (rho * omega);
			p->addScaled(-omega, ap);
			p->scaleThenAdd(beta, resid);

			num_its++;
			rho = rho_new;
		}
		return num_its;
	}
};
