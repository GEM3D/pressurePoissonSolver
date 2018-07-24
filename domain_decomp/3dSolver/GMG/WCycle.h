#ifndef GMGWCycle_H
#define GMGWCycle_H
#include "GMG/Cycle.h"
#include "tpl/json.hpp"
namespace GMG
{
/**
 * @brief Implementation of a W-cycle
 */
class WCycle : public Cycle
{
	private:
	int num_pre_sweeps    = 1;
	int num_post_sweeps   = 1;
	int num_coarse_sweeps = 1;
	int num_mid_sweeps = 1;

	protected:
	/**
	 * @brief Implements W-cycle. Pre-smooth, visit coarser level, smooth, visit coarse level, and then post-smooth.
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
			for (int i = 0; i < num_mid_sweeps; i++) {
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
	 * @brief Create new W-cycle
	 *
	 * @param finest_level a pointer to the finest level
	 */
	WCycle(std::shared_ptr<Level> finest_level, nlohmann::json config_j) : Cycle(finest_level)
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
		try {
			num_mid_sweeps = config_j.at("mid_sweeps");
		} catch (nlohmann::detail::out_of_range oor) {
		}
	}
};
} // namespace GMG
#endif
