// AUTOMATICALLY GENERATED FILE.  DO NOT MODIFY.  Place custom code in custominit.h.
#include "custominit.h"
void mcExitClass_HxModuleMaterial();
void mcInitClass_HxModuleMaterial();
void mcExitClass_HxLineRaycast();
void mcInitClass_HxLineRaycast();


extern "C"
#ifdef _WIN32
__declspec(dllexport)
#endif
void
amirapackage_hxlineviewer_init()
{
    static bool isInitialized = false;
    if (isInitialized)
      return;

    isInitialized = true;

    mcInitClass_HxModuleMaterial();
    mcInitClass_HxLineRaycast();

    custominit();
}


extern "C"
#ifdef _WIN32
__declspec(dllexport)
#endif
void
amirapackage_hxlineviewer_finish()
{
    static bool isFinished = false;
    if (isFinished)
      return;

    isFinished = true;


    mcExitClass_HxLineRaycast();
    mcExitClass_HxModuleMaterial();
}

#if defined(_WIN32)
#  include <windows.h>


BOOL WINAPI DllMain(
    __in  HINSTANCE hinstDLL,
    __in  DWORD fdwReason,
    __in  LPVOID lpvReserved
    )
{
    switch (fdwReason)
    {
    case DLL_PROCESS_ATTACH:
        amirapackage_hxlineviewer_init();
        break;
    case DLL_PROCESS_DETACH:
        amirapackage_hxlineviewer_finish();
        break;
    default:
        ;
    }
    return true;
}


#endif

#if defined(__GNUC__)
void
__attribute__((constructor)) soconstructor_hxlineviewer()
{
    amirapackage_hxlineviewer_init();
}

void
__attribute__((destructor)) sodestructor_hxlineviewer()
{
    amirapackage_hxlineviewer_finish();
}
#endif
