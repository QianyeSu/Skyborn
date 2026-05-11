#ifndef PROGRESS_H
#define PROGRESS_H

constexpr unsigned progressMinSize = 99999;

namespace progress
{

void set_context_function(const char *(*func)(void) );

}  // namespace progress

namespace cdo
{

extern bool ProgressInUse;

class Progress
{
private:
  bool isActiv{ false };
  bool contextActive{ false };
  int contextLen{ 0 };
  int value{ -1 };
  const char *context = "";

  void init();

public:
  Progress(int processID = 0)
  {
    if (processID == 0 && !ProgressInUse)
      {
        ProgressInUse = true;
        isActiv = true;
        init();
      }
  }
  ~Progress()
  {
    update(1.0);
    if (isActiv) ProgressInUse = false;
  }

  void update(double curval, double offset = 0.0, double refval = 1.0);
};

}  // namespace cdo

#endif
