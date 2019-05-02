template <size_t D> class SchurDomainOp : public Operator<D>
{
	private:
	/**
	 * @brief PETSc Matrix object
	 */
	std::shared_ptr<SchurHelper<D>> helper;

	public:
	/**
	 * @brief Crate new WrapOp
	 *
	 * @param matrix the PETSc matrix
	 */
	SchurDomainOp(std::shared_ptr<SchurHelper<D>> helper)
	{
		this->helper = helper;
	}
	/**
	 * @brief Perform matrix/vector multiply.
	 *
	 * @param x the input vector.
	 * @param b the output vector.
	 */
	void apply(std::shared_ptr<const Vector<D>> x, std::shared_ptr<Vector<D>> b) const
	{
		helper->apply(x, b);
	}
};
