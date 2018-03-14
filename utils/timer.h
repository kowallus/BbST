#ifndef _TIMER_H
#define _TIMER_H

#include <chrono>

using namespace std;
using namespace chrono;

typedef struct {
	steady_clock::time_point start;
	steady_clock::time_point stop;
} stopWatch;

class ChronoStopWatch {

private:
	stopWatch timer;
public:
	ChronoStopWatch() {};
	void startTimer();
	void stopTimer();
	double getElapsedTime();
};

#endif /* TIMER_H_ */
