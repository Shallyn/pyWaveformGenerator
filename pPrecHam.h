/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_PPRECHAM__
#define __INCLUDE_PPRECHAM__
#include "pHamiltonian.h"
#include "pPrecUtils.h"
int PrecHcapNumericalDerivative(double t,
                                const REAL8 values[],
                                REAL8 dvalues[],
                                void *funcParams);

#endif

