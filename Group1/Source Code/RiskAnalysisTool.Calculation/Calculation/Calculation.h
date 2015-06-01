#pragma once

#if defined _RISKANALYSISTOOL_CALCULATION_DLL
#define _RISKANALYSISTOOL_CALCULATION_API __declspec(dllexport)
#elif defined _RISKANALYSISTOOL_CALCULATION_LIB
#define _RISKANALYSISTOOL_CALCULATION_API 
#else
#define _RISKANALYSISTOOL_CALCULATION_API __declspec(dllimport)
#endif

#ifdef __cplusplus

#ifndef EXTERN_C
#define EXTERN_C    extern "C"
#endif

#else /* __cplusplus */

#ifndef EXTERN_C
#define EXTERN_C    extern
#endif

#endif /* __cplusplus */