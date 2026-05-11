#include <cstdio>
#include <algorithm>

#include "progress.h"
#include "cdo_options.h"

namespace progress
{

const char *(*getContext)(void) = nullptr;
/**
 * parameter p_context:
 *          will be displayed in status message and indicates sub process for
 *          which the progress is displayed.
 */
void
set_context_function(const char *(*func)(void) )
{
  getContext = func;
}

}  // namespace progress

namespace cdo
{

bool ProgressInUse = false;

void
Progress::init()
{
  if (progress::getContext != nullptr) context = progress::getContext();
}

void
Progress::update(double curval, double offset, double refval)
{
  if (!isActiv) return;
  if (!cdo::stdoutIsTerminal || Options::silentMode || Options::cdoVerbose) return;

  curval = std::clamp(curval, 0.0, 1.0);
  offset = std::clamp(offset, 0.0, 1.0);
  refval = std::clamp(refval, 0.0, 1.0);

  int newValue = (offset + refval * curval) * 100;

  if (value == -1)
  {
    contextLen = std::fprintf(stdout, "%s: %3d%%", context, 0);
    fflush(stdout);
    contextActive = true;
  }

  if (newValue != value)
  {
    value = newValue;
    std::fprintf(stdout, "\b\b\b\b%3d%%", value);
    fflush(stdout);
  }

  if (value == 100 && contextActive)
  {
    contextActive = false;
    while (contextLen--) std::fprintf(stdout, "\b \b");
    fflush(stdout);
  }
}

}  // namespace cdo
