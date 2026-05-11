#ifndef CDO_OUTPUT_H
#define CDO_OUTPUT_H

#include "mpmo.h"
#include "cdo_options.h"

// Debug Switches
extern int cdoDebug;
extern int cdoDebugExt;  //  Debug level for the KNMI extensions
                         // Subsystem Debug Switches
extern unsigned PROCESS;
extern unsigned PIPE;
extern unsigned PIPE_STREAM;
extern unsigned FILE_STREAM;
extern unsigned PTHREAD;
extern unsigned PROCESS_MANAGER;
extern unsigned PIPE;
extern unsigned CDO_NODE;
extern unsigned PARSER;
extern unsigned PROCESS_INT;
extern unsigned CDO_DEBUG;
extern unsigned FACTORY;
extern unsigned KVLIST;
extern unsigned MODULE_INFO;
extern unsigned ARGUMENTS;

extern std::string debug_option_string;

void print_debug_options();
void query_user_exit(std::string const &argument);

namespace cdo
{
void parse_debug_arguments(std::vector<std::string> const &tokens, unsigned &cdoDebugLevel, unsigned &cdiDebugLevel);
void print_debug_levels(unsigned cdoDebugLevel, unsigned cdiDebugLevel);
void set_debug(unsigned p_debug_level);
bool dbg();
extern void (*exitProgram)(std::string);
extern const char *(*getContext)(void);
void set_exit_function(void (*func)(std::string msg));
void set_context_function(const char *(*func)(void) );
}  // namespace cdo

void cdi_open_error(int cdiErrno, std::string const &format, const char *path);
std::string cdo_argv_to_string(std::vector<std::string> const &argv);

template <typename... Args>
void
cdo_abort(std::string const &format, Args const &...args)
{
  fflush(stdout);
  std::string errmsg = MpMO::PrintCerr(Red("\n%s (Abort): ") + format, cdo::getContext(), args...);
  if (MpMO::exitOnError) cdo::exitProgram(errmsg);
}

template <typename... Args>
void
cdo_error(std::string const &format, Args const &...args) noexcept
{
  fflush(stdout);
  MpMO::PrintCerr(Red("\n%s (Abort): ") + format, cdo::getContext(), args...);
  Options::cdoExitStatus = 1;
}

template <typename... Args>
void
cdo_warning(std::string const &format, Args const &...args) noexcept
{
  if (MpMO::warningsEnabled)
    {
      if (MpMO::pedantic)
        {
          MpMO::PrintCerr(Red("%s (Warning): ") + format, cdo::getContext(), args...);
          if (MpMO::exitOnError) cdo::exitProgram("cdo_warning (pedantic)");
        }
      else { MpMO::PrintCerr(Yellow("%s (Warning): ") + format, cdo::getContext(), args...); }
    }
}

template <typename... Args>
void
cdo_verbose(std::string const &format, Args const &...args) noexcept
{
  MpMO::PrintCerr(Green("%s: ") + format, cdo::getContext(), args...);
}

template <typename... Args>
void
cdo_print(std::string const &format, Args const &...args) noexcept
{
  if (!MpMO::silentMode) MpMO::Print(Green("%s: ") + format, cdo::getContext(), args...);
}

#ifdef WITH_CALLER_NAME
#define cdo_sys_error(...) MpMO::SysError_(__func__, __VA_ARGS__)
#else
#define cdo_sys_error(...) MpMO::SysError_("", __VA_ARGS__)
#endif

#endif
