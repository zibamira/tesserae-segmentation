// AUTOMATICALLY GENERATED FILE.  DO NOT MODIFY.
#ifndef HXLINEVIEWER_API_EXPORT_MACRO_H
#define HXLINEVIEWER_API_EXPORT_MACRO_H

#ifdef HXLINEVIEWER_STATIC
#   define HXLINEVIEWER_API
#else
#   ifdef _WIN32
#       ifdef HXLINEVIEWER_EXPORTS
#           define HXLINEVIEWER_API __declspec(dllexport)
#       else
#           define HXLINEVIEWER_API __declspec(dllimport)
#       endif
#   else
#       define HXLINEVIEWER_API __attribute__ ((visibility("default")))
#   endif
#endif

#endif
