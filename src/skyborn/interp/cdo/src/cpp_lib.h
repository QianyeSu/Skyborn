/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef CPP_LIB_H
#define CPP_LIB_H

#include <version>
#include <ranges>
#include <span>

#if __cpp_lib_ranges >= 201911L
#define HAVE_LIB_RANGES 1
#else
#error C++20 Ranges library is Not Available!
#endif

#if __cpp_lib_span >= 202002L
#define HAVE_LIB_SPAN 1
#else
#error C++20 Span library is Not Available!
#endif

#if __cpp_lib_ranges_zip >= 202110L
#define HAVE_LIB_RANGES_ZIP 1
#endif

#if __cpp_lib_mdspan >= 202207L
#define HAVE_LIB_MDSPAN 1
#endif

#endif
