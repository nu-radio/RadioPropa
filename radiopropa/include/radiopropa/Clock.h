#ifndef RADIOPROPA_CLOCK_H
#define RADIOPROPA_CLOCK_H

namespace radiopropa {

//class ClockImpl;

class Clock {
private:
	class Impl;
	Impl *impl;
public:
	Clock();
	virtual ~Clock();

	void reset();
	double getSecond();
	double getMillisecond();
	static Clock &getInstance();
};

} // namespace radiopropa

#endif // RADIOPROPA_CLOCK_H
