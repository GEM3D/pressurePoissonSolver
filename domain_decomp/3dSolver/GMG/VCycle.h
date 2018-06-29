#ifndef GMGVCycle_H
#define GMGVCycle_H
#include "GMG/Cycle.h"
#include "tpl/json.hpp"
namespace GMG
{
/**
 * @brief Implementation of a V-cycle
 */
class VCycle : public Cycle
{
	private:
	int num_pre_sweeps    = 1;
	int num_post_sweeps   = 1;
	int num_coarse_sweeps = 1;

	protected:
	/**
	 * @brief Implements V-cycle. Pre-smooth, visit coarser level and then post-smooth.
	 *
	 * @param level the current level that is being visited.
	 */
	void visit(const Level &level)
	{
		if (level.coarsest()) {
			for (int i = 0; i < num_coarse_sweeps; i++) {
				smooth(level);
			}
		} else {
			for (int i = 0; i < num_pre_sweeps; i++) {
				smooth(level);
			}
			prepCoarser(level);
			visit(level.getCoarser());
			for (int i = 0; i < num_post_sweeps; i++) {
				smooth(level);
			}
		}
		if (!level.finest()) { prepFiner(level); }
	}

	public:
	/**
	 * @brief Create new V-cycle
	 *
	 * @param finest_level a pointer to the finest level
	 */
	VCycle(std::shared_ptr<Level> finest_level, nlohmann::json config_j) : Cycle(finest_level)
	{
		try {
			num_pre_sweeps = config_j.at("pre_sweeps");
		} catch (nlohmann::detail::out_of_range oor) {
		}
		try {
			num_post_sweeps = config_j.at("post_sweeps");
		} catch (nlohmann::detail::out_of_range oor) {
		}
		try {
			num_coarse_sweeps = config_j.at("coarse_sweeps");
		} catch (nlohmann::detail::out_of_range oor) {
		}
	}
};
} // namespace GMG
#endif
