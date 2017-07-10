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
	vector<string> order;
	map<string, vector<double>>           times;
	map<string, steady_clock::time_point> starts;

	public:
	void start(string name)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		starts[name] = steady_clock::now();
	}
	void stop(string name)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		duration<double> time = steady_clock::now() - starts[name];
		times[name].push_back(time.count());
		if (std::find(order.begin(), order.end(), name) == order.end()) {
			order.push_back(name);
		}
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
				os << "  times (sec):";
				for (double t : times) {
					os << " " << t;
				}
				os << endl << endl;
			}
		}
		return os;
	}
};
}
#endif
