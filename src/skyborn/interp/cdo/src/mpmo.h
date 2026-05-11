/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida
          Oliver Heidmann

*/
#ifndef MPMO_H
#define MPMO_H

#include <string>
#include <iostream>
#include <cstring>
#include <functional>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "mpmo_color.h"

namespace MpMO
{
/* ---- CAUTION ---- */
/* In this entire namespace the warnings -Wformat-nonliteral and -Wformat-security are disabled for all printf statements
 * this can cause security issues IF the users have access to the format strings used in this files.
 * in CDO this is not the case and as such it was decided to use this file as it is and to disable the warnings.
 *
 * */

extern unsigned DebugLevel;
extern bool silentMode;
void enable_silent_mode(bool enable);
extern bool warningsEnabled;
void enable_warnings(bool enable);
extern bool verbose;
void enable_verbose(bool enable);
extern bool exitOnError;
extern bool pedantic;
void enable_pedantic(bool enable);
extern int padding_width;

template <typename T>
T
Argument(T value) noexcept
{
  return value;
}

template <typename T>
T const *
Argument(std::basic_string<T> const &value) noexcept
{
  return value.c_str();
}

template <typename... Args>
void
Print(std::string const &format, Args const &...args) noexcept
{
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat-nonliteral"
#pragma GCC diagnostic ignored "-Wformat-security"
  if (!silentMode) printf((format + "\n").c_str(), Argument(args)...);
#pragma GCC diagnostic pop
}

template <typename... Args>
std::string
PrintCerr(std::string const &format, Args const &...args)
{
  int cx;
  cx = snprintf(nullptr, 0, format.c_str(), Argument(args)...);
  std::string msg = std::string(cx + 1, '\0');
  snprintf(&msg[0], cx + 1, format.c_str(), Argument(args)...);
  std::fprintf(stderr, "%s\n", msg.c_str());
  return std::string(msg);
}

std::string debug_scope_string(const char *p_file, const char *p_func, int p_line, const char *context);

void Debug_(const char *p_file, const char *p_func, int p_line, const char *context, int p_debugScope,
            std::function<void()> p_function);
void Debug_(const char *p_file, const char *p_func, int p_line, const char *context, std::function<void()> p_function);

template <typename... Args>
void
Debug_(const char *p_file, const char *p_func, int p_line, const char *context, int p_debugScope, std::string const &format,
       Args const &...args)
{
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat-nonliteral"
#pragma GCC diagnostic ignored "-Wformat-security"
  if (p_debugScope)
  {
    std::fprintf(stderr, (debug_scope_string(p_file, p_func, p_line, context) + format + "\n").c_str(), Argument(args)...);
  }
#pragma GCC diagnostic pop
}

template <typename... Args>
void
Debug_(const char *p_file, const char *p_func, int p_line, const char *context, std::string const &format, Args const &...args)
{
  Debug_(p_file, p_func, p_line, context, DebugLevel, format, args...);
}

template <typename... Args>
void
Error_(const char *caller, std::string const &format, Args const &...args) noexcept
{
  PrintCerr(Red("Error:") + "(%s)" + format, caller, Argument(args)...);
  if (exitOnError) std::exit(EXIT_FAILURE);
}

void Verbose_(std::function<void()> p_function) noexcept;

template <typename... Args>
void
Verbose_(std::string const &format, Args const &...args) noexcept
{
  if (verbose) PrintCerr(format, Argument(args)...);
}

template <typename... Args>
void
Warning_(const char *caller, std::string const &format, Args const &...args) noexcept
{
  (void) caller;  // quell warning if WITH_CALLER_NAME is not defined
  if (warningsEnabled)
  {
    if (pedantic)
    {
      PrintCerr(Red("Warning: ") + format, Argument(args)...);
      if (exitOnError) std::exit(EXIT_FAILURE);
    }
    else
    {
      PrintCerr(Yellow("Warning: ") + format, Argument(args)...);
    }
  }
}

template <typename... Args>
[[noreturn]] void
SysError_(const char *func, std::string const &format, Args const &...args) noexcept
{
  int saved_errno = errno;
  PrintCerr(Red("SysError: %s ") + format, func, Argument(args)...);
  if (saved_errno)
  {
    errno = saved_errno;
    std::perror("System error message");
  }
  std::exit(EXIT_FAILURE);
}

}  // namespace MpMO

#define Verbose(...) Verbose_(__VA_ARGS__)
#ifndef NO_DEBUG
#define Debug(...) MpMO::Debug_(__FILE__, __func__, __LINE__, cdo::getContext(), __VA_ARGS__)
#else
#define Debug(...)
#endif

#endif
