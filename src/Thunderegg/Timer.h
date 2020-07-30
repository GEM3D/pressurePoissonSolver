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

#ifndef TIMER_H
#define TIMER_H
#include <chrono>
#include <map>
#include <string>
#include <vector>
namespace Tools
{
using namespace std;
using namespace std::chrono;
class Timer
{
	private:
	bool                                  verbose = false;
	vector<string>                        order;
	map<string, vector<double>>           times;
	map<string, steady_clock::time_point> starts;

	public:
	void start(string name)
	{
		if (verbose) { cout << "Starting " << name << endl; }
		MPI_Barrier(MPI_COMM_WORLD);
		starts[name] = steady_clock::now();
	}
	void stop(string name)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		duration<double> time = steady_clock::now() - starts[name];
		times[name].push_back(time.count());
		if (std::find(order.begin(), order.end(), name) == order.end()) { order.push_back(name); }
		if (verbose) { cout << "Stopped " << name << endl; }
	}
	friend ostream &operator<<(ostream &os, const Timer &timer)
	{
		os << endl;
		os << "TIMING RESULTS" << endl;
		os << "==============" << endl << endl;

		for (string name : timer.order) {
			vector<double> times = timer.times.at(name);

			os << name << endl;
			os << string(name.size(), '-') << endl;

			if (times.size() == 1) {
				os << "   time (sec): " << times[0] << endl << endl;
			} else {
				double average = 0;
				for (double t : times) {
					average += t;
				}
				average /= times.size();

				os << "average (sec): " << average << endl;
				if (times.size() < 10) {
					os << "  times (sec):";
					for (double t : times) {
						os << " " << t;
					}
					os << endl;
				}
				os << endl;
			}
		}
		return os;
	}
};
} // namespace Tools
#endif
