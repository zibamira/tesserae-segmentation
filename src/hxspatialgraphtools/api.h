// AUTOMATICALLY GENERATED FILE.  DO NOT MODIFY.
#ifndef HXSPATIALGRAPHTOOLS_API_EXPORT_MACRO_H
#define HXSPATIALGRAPHTOOLS_API_EXPORT_MACRO_H

#ifdef HXSPATIALGRAPHTOOLS_STATIC
#   define HXSPATIALGRAPHTOOLS_API
#else
#   ifdef _WIN32
#       ifdef HXSPATIALGRAPHTOOLS_EXPORTS
#           define HXSPATIALGRAPHTOOLS_API __declspec(dllexport)
#       else
#           define HXSPATIALGRAPHTOOLS_API __declspec(dllimport)
#       endif
#   else
#       define HXSPATIALGRAPHTOOLS_API __attribute__ ((visibility("default")))
#   endif
#endif

#endif
