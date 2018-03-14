#include "timer.h"

void ChronoStopWatch::startTimer() {
	timer.start = steady_clock::now();
}

void ChronoStopWatch::stopTimer() {
	timer.stop = steady_clock::now();
}

double ChronoStopWatch::getElapsedTime() {

	nanoseconds time_span = duration_cast<nanoseconds>(timer.stop - timer.start);
	return (double)time_span.count() / 1000000000.0;
}
