// AUTOMATICALLY GENERATED FILE.  DO NOT MODIFY.
#ifndef HXCONTOURTREE_API_EXPORT_MACRO_H
#define HXCONTOURTREE_API_EXPORT_MACRO_H

#ifdef HXCONTOURTREE_STATIC
#   define HXCONTOURTREE_API
#else
#   ifdef _WIN32
#       ifdef HXCONTOURTREE_EXPORTS
#           define HXCONTOURTREE_API __declspec(dllexport)
#       else
#           define HXCONTOURTREE_API __declspec(dllimport)
#       endif
#   else
#       define HXCONTOURTREE_API __attribute__ ((visibility("default")))
#   endif
#endif

#endif
