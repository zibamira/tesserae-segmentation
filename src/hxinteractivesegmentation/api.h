// AUTOMATICALLY GENERATED FILE.  DO NOT MODIFY.
#ifndef HXINTERACTIVESEGMENTATION_API_EXPORT_MACRO_H
#define HXINTERACTIVESEGMENTATION_API_EXPORT_MACRO_H

#ifdef HXINTERACTIVESEGMENTATION_STATIC
#   define HXINTERACTIVESEGMENTATION_API
#else
#   ifdef _WIN32
#       ifdef HXINTERACTIVESEGMENTATION_EXPORTS
#           define HXINTERACTIVESEGMENTATION_API __declspec(dllexport)
#       else
#           define HXINTERACTIVESEGMENTATION_API __declspec(dllimport)
#       endif
#   else
#       define HXINTERACTIVESEGMENTATION_API __attribute__ ((visibility("default")))
#   endif
#endif

#endif
