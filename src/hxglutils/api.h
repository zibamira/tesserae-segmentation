// AUTOMATICALLY GENERATED FILE.  DO NOT MODIFY.
#ifndef HXGLUTILS_API_EXPORT_MACRO_H
#define HXGLUTILS_API_EXPORT_MACRO_H

#ifdef HXGLUTILS_STATIC
#   define HXGLUTILS_API
#else
#   ifdef _WIN32
#       ifdef HXGLUTILS_EXPORTS
#           define HXGLUTILS_API __declspec(dllexport)
#       else
#           define HXGLUTILS_API __declspec(dllimport)
#       endif
#   else
#       define HXGLUTILS_API __attribute__ ((visibility("default")))
#   endif
#endif

#endif
