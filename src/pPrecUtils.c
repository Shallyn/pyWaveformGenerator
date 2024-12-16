/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "pPrecUtils.h"

INT g_PrecFlag = 1;
INT get_PrecFlag()
{
    return g_PrecFlag;
}
void set_PrecFlag(INT flag)
{
    g_PrecFlag = flag;
    return;
}

void prec_SetSEOBPrecVariables(SEOBPrecVariables *vars, 
    const REAL8 xVec[], const REAL8 pTVec[], 
    const REAL8 s1Vec[], const REAL8 s2Vec[])
{
    memcpy(vars->xVec, xVec, 3*sizeof(REAL8));
    memcpy(vars->pTVec, pTVec, 3*sizeof(REAL8));
    memcpy(vars->s1Vec, s1Vec, 3*sizeof(REAL8));
    memcpy(vars->s2Vec, s2Vec, 3*sizeof(REAL8));
    vars->r = sqrt(inner_product3d(xVec, xVec));
    vars->prT = inner_product3d(xVec, pTVec) / vars->r;
    vars->pT2 = inner_product3d(pTVec, pTVec);
    vars->xS1 = inner_product3d(xVec, s1Vec);
    vars->xS2 = inner_product3d(xVec, s2Vec);
    vars->pTS1 = inner_product3d(pTVec, s1Vec);
    vars->pTS2 = inner_product3d(pTVec, s2Vec);
    vars->S1S1 = inner_product3d(s1Vec, s1Vec);
    vars->S1S2 = inner_product3d(s1Vec, s2Vec);
    vars->S2S2 = inner_product3d(s2Vec, s2Vec);    
	vars->LVec[0] = xVec[1] * pTVec[2] - xVec[2] * pTVec[1];
	vars->LVec[1] = xVec[2] * pTVec[0] - xVec[0] * pTVec[2];
	vars->LVec[2] = xVec[0] * pTVec[1] - xVec[1] * pTVec[0];
    vars->xpS1 = inner_product3d(vars->LVec, s1Vec);
    vars->xpS2 = inner_product3d(vars->LVec, s2Vec);
    return;
}
