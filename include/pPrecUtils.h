/**
 * Writer: Xiaolin.liu
 * xiaolin.liu@mail.bnu.edu.cn
 *
 * This module contains basic functions for  calculation.
 * Functions list:
 * Kernel:
 * 20xx.xx.xx, LOC
 **/

#ifndef __INCLUDE_PPRECUTILS__
#define __INCLUDE_PPRECUTILS__
#include "pUtils.h"

INT get_PrecFlag();
void set_PrecFlag(INT flag);
#define PREC_FLAG (get_PrecFlag())
#define SET_PREC_FLAG(flag)                                                                                            \
    do                                                                                                                 \
    {                                                                                                                  \
        set_PrecFlag(ver);                                                                                             \
    } while (0)

void prec_SetSEOBPrecVariables(SEOBPrecVariables *vars, const REAL8 xVec[], const REAL8 pTVec[], const REAL8 s1Vec[],
                               const REAL8 s2Vec[]);

#endif
