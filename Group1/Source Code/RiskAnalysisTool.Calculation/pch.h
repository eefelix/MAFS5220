// pch.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#ifdef _MSC_VER
#include <sdkddkver.h>
#endif

#if defined _WINDOWS || defined _WINDLL

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN             // Exclude rarely-used stuff from Windows headers
#endif

// Windows Header Files:
#include <Windows.h>

#endif

#ifdef __cplusplus
#include <cstdlib>
#include <cassert>

#include <memory>
#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>

#include <ql/qldefines.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/enable_shared_from_this.hpp>

#else

#include <stdlib.h>
#include <assert.h>

#endif

