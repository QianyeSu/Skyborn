/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifndef CDO_TIMER_H
#define CDO_TIMER_H

#include <string>
#include <chrono>

namespace cdo
{

using namespace std::chrono;
using clock = steady_clock;

class timer
{
public:
  timer() { reset(); }

  void
  reset()
  {
    startPoint = clock::now();
  }

  double
  time_span() const
  {
    auto timeSpan = duration_cast<duration<double>>(clock::now() - startPoint);
    return timeSpan.count();
  }

  double
  elapsed(bool resetTimer = false)
  {
    auto timeSpan = time_span();
    if (resetTimer) reset();
    return timeSpan;
  }

private:
  clock::time_point startPoint;
};

// interval timer (stop watch)
static bool timerNeedInit{ true };
static double timerShift{ 0.0 };  // minimal internal time needed to do one measurement

// interval timer (stop watch)
class iTimer
{
private:
  double
  get_shift(void) const
  {
    constexpr int numTests = 100;
    double dt0 = 1.0;
    for (int i = 0; i < numTests; ++i)
      {
        auto now = clock::now();
        auto dt = get_time_val(now);
        dt0 = std::min(dt0, dt);
      }

    return dt0;
  }

  void
  timer_init(void)
  {
    timerShift = get_shift();
    timerNeedInit = false;
  }

  double
  get_time_val(const clock::time_point &_startPoint) const
  {
    auto dt = duration_cast<duration<double>>(clock::now() - _startPoint);
    return dt.count();
  }

  clock::time_point startPoint;
  bool isRunning{ false };

public:
  iTimer()
  {
    if (timerNeedInit) timer_init();
  }
  explicit iTimer(std::string const &_name) : name(_name)
  {
    if (timerNeedInit) timer_init();
  }

  void
  start()
  {
    if (isRunning) std::fprintf(stderr, "timer::start: timer::stop call missing\n");
    isRunning = true;
    startPoint = clock::now();
  }

  void
  stop()
  {
    if (!isRunning) std::fprintf(stderr, "timer::stop: timer::start call missing\n");

    auto dt = get_time_val(startPoint);
    dt -= timerShift;

    sum += dt;
    min = std::min(min, dt);
    max = std::max(max, dt);
    calls++;
    isRunning = false;
  }

  double
  elapsed()
  {
    if (isRunning) stop();
    return sum;
  }

  int calls{ 0 };
  int stat{ 0 };
  double sum{ 0.0 };
  double min{ 1.e30 };
  double max{ 0.0 };
  std::string name;
};

extern iTimer readTimer;
extern iTimer writeTimer;

inline double
get_wtime()
{
  return duration_cast<duration<double>>(clock::now().time_since_epoch()).count();
}

}  // namespace cdo

#endif /* CDO_TIMER_H */
