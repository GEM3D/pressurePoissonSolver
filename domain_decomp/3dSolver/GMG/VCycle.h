#ifndef GMGVCycle_H
#define GMGVCycle_H
namespace GMG
{
/**
 * @brief Implementation of a V-cycle
 */
class VCycle : public Cycle
{
	protected:
	/**
	 * @brief Implements V-cycle. Pre-smooth, visit coarser cycle and then post-smooth.
	 *
	 * @param level the current level that is being visited.
	 */
	void visit(const Level &level)
	{
		if (level.coarsest()) {
			smooth(level);
		} else {
			smooth(level);
			prepCoarser(level);
			visit(level.getCoarser());
			smooth(level);
		}
		if (!level.finest()) { prepFiner(level); }
	}

	public:
    /**
     * @brief Create new V-cycle
     *
     * @param finest_level a pointer to the finest level
     */
	VCycle(std::shared_ptr<Level> finest_level) : Cycle(finest_level) {}
};
} // namespace GMG
#endif
