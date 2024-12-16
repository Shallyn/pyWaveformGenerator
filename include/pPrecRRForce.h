/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_PPRECRRFORCE__
#define __INCLUDE_PPRECRRFORCE__
#include "pPrecUtils.h"

void prec_CalculateRRForce(SEOBPrecVariables *vars, 
    REAL8 eta, REAL8 prTDot,
    REAL8 fRRVec[], REAL8 fRRs1Vec[], REAL8 fRRs2Vec[]);

#endif

