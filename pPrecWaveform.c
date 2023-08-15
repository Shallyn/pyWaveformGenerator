/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "myLog.h"
#include "pPrecWaveform.h"
#include "pHamiltonian.h"
#include <gsl/gsl_sf_gamma.h>

void prec_CalculateSEOBPrecWaveformVariables(SEOBPrecWaveformVariables *vars,
    REAL8 nchia, REAL8 nchis, REAL8 lchia, REAL8 lchis, REAL8 echia, REAL8 echis,
    REAL8 chi1chi1, REAL8 chi1chi2, REAL8 chi2chi2, REAL8 Jn, REAL8 Jl, REAL8 Je,
    REAL8 r, REAL8 prT, REAL8 prTDot)
{
    REAL8 sqr;
    REAL8 x, x2, x3, x4;
    REAL8 y, y2, y3, y4, y5, y6;
    REAL8 w;
    REAL8 z, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13, z14, z15, z16, z18, z20, z21, z22, z24, z30;
    vars->nchia = nchia;
    vars->nchis = nchis;
    vars->lchia = lchia;
    vars->lchis = lchis;
    vars->echia = echia;
    vars->echis = echis;
    vars->chi1chi1 = chi1chi1;
    vars->chi1chi2 = chi1chi2;
    vars->chi2chi2 = chi2chi2;
    vars->Jn = Jn;
    vars->Je = Je;
    vars->Jl = Jl;
    sqr = sqrt(r);
    vars->x = x = 1./sqr;
    vars->x2 = x2 = 1./r;
    vars->x3 = x3 = x2*x;
    vars->x4 = x4 = x3*x;
    vars->y = y = sqr*prT;
    vars->y2 = y2 = y*y;
    vars->y3 = y3 = y2*y;
    vars->y4 = y4 = y3*y;
    vars->y5 = y5 = y4*y;
    vars->y6 = y6 = y5*y;
    w = 1. + r*r*prTDot;
    if (w>0.) {z2 = cbrt(w); z = sqrt(z2);}
    if (w<0.){ z2 = z = 0.;}
    vars->z2 = z2;
    vars->z3 = z3 = z2*z;
    vars->z4 = z4 = z2*z2;
    vars->z5 = z5 = z4*z;
    vars->z6 = z6 = z4*z2;
    vars->z7 = z7 = z6*z;
    vars->z8 = z8 = z6*z2;
    vars->z9 = z9 = z8*z;
    vars->z10 = z10 = z8*z2;
    vars->z11 = z11 = z10*z;
    vars->z12 = z12 = z10*z2;
    vars->z13 = z13 = z12*z;
    vars->z14 = z14 = z12*z2;
    vars->z15 = z15 = z14*z;
    vars->z16 = z16 = z14*z2;
    vars->z18 = z18 = z16*z2;
    vars->z20 = z20 = z18*z2;
    vars->z21 = z21 = z20*z;
    vars->z22 = z22 = z20*z2;
    vars->z24 = z24 = z22*z2;
    vars->z30 = z30 = z24*z6;    
    return;
}

/* ---------------------------------------------------------------------------- */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/*                                    V1                                        */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/* ---------------------------------------------------------------------------- */
static COMPLEX16 prec_CalculateWaveformFactor_mode22_v1(SEOBPrecWaveformVariables *vars,
    REAL8 eta, REAL8 dm)
{
    // REAL8 x, x2, x3, x4, x5, x6, x7, x8;
    REAL8 eta2 = eta*eta;
    REAL8 nchia, nchis, lchia, lchis, echia, echis;
    REAL8 nchia2, nchis2, lchia2, lchis2, echia2, echis2;
    REAL8 chi1chi1, chi1chi2, chi2chi2;
    REAL8 Jn, Jl, Je;
    REAL8 Jn2, Jl2, Je2;
    REAL8 x, x2, x3, x4;
    REAL8 y, y2, y3, y4, y5, y6;
    REAL8 z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13, z14, z15, z16, z18, z20, z21, z22, z24, z30;
    nchia = vars->nchia;
    nchis = vars->nchis;
    lchia = vars->lchia;
    lchis = vars->lchis;
    echia = vars->echia;
    echis = vars->echis;
    nchia2 = nchia*nchia;
    nchis2 = nchis*nchis;
    lchia2 = lchia*lchia;
    lchis2 = lchis*lchis;
    echia2 = echia*echia;
    echis2 = echis*echis;
    chi1chi1 = vars->chi1chi1;
    chi1chi2 = vars->chi1chi2;
    chi2chi2 = vars->chi2chi2;
    Jn = vars->Jn;
    Jl = vars->Jl;
    Je = vars->Je;
    Jn2 = Jn*Jn;
    Jl2 = Jl*Jl;
    Je2 = Je*Je;
    x = vars->x;x2 = vars->x2;x3 = vars->x3;x4 = vars->x4;
    y = vars->y;y2 = vars->y2;y3 = vars->y3;y4 = vars->y4;y5 = vars->y5;y6 = vars->y6;
    z2 = vars->z2;z3 = vars->z3;z4 = vars->z4;z5 = vars->z5;z6 = vars->z6;
    z7 = vars->z7;z8 = vars->z8;z9 = vars->z9;z10 = vars->z10;z11 = vars->z11;
    z12 = vars->z12;z13 = vars->z13;z14 = vars->z14;z15 = vars->z15;z16 = vars->z16;
    z18 = vars->z18;z20 = vars->z20;z21 = vars->z21;z22 = vars->z22;z24 = vars->z24;
    z30 = vars->z30;
    COMPLEX16 fp0, fp2, fp3, fp4;
    fp0 = (1. - y2 - 2.*z2 + 2.*I*y*z3 + z6)/(2.*z2);
    fp2 = (x2*(14 + 172*z10 - 3*z12 - 8*z18 + y4*(-7 + z6) - 175*z6 + 42*lchia*z6 + 42*dm*lchis*z6 + 42*I*nchia*z6 + 42*I*dm*nchis*z6 - I*y*z3*(14 + 37*z12 + 13*z6) + y2*(-7 + 49*z12 + 18*z6) + I*y3*z3*(-7 + 19*z6) + eta*(-14 - 110*z10 + 30*z12 + 31*z18 + I*y*z3*(14 + 41*z12 - 101*z6) + 63*z6 + y2*(21 + 49*z12 + 37*z6) - y4*(7 + 38*z6) + I*y3*z3*(-7 + 97*z6))))/(84*z8);
    fp3 = (-(x3*(6*I*eta*Je*(Jl*z3*(lchis*y - nchis*z3) + Jn*(lchis + lchis*y2 - nchis*y*z3))*(-1 + y2 - 2*I*y*z3 - z6) - echis*(Jl2 + Jn2)*(eta*(-1 + y2 + 5*I*y*z3 + 9*z6 - 8*z8) - 2*(-1 + y2 + 4*I*y*z3 + 5*z6 - 4*z8)) + 2*dm*echia*(Jl2 + Jn2)*(-1 + y2 + 4*I*y*z3 + 5*z6 - 4*z8))))/(6*(Jl2 + Jn2)*z5);
    fp4 = (x4*(I*(3024)*Je*Jl*(echia*(-1 + 4*eta)*nchia + dm*echia*(-1 + 2*eta)*nchis + echis*(-1 + 2*eta)*(dm*nchia + nchis - 2*eta*nchis))*z12*(-1 + y2 - 2*I*y*z3 - z6) + Jl2*(I*(-3)*y5*z3*(35 + 245*z12 - 88*z6) + 6*y6*(-28 + 11*z12 + 65*z6) - 6*I*y3*z3*(70 - 902*z12 + 283*z18 + (109 - 84*chi2chi2*(-1 + dm) + 84*chi1chi1*(1 + dm) - 168*nchia2 - 336*dm*nchia*nchis - 168*nchis2)*z6) - 
        18*y4*(28 - 39*z12 + 20*z18 + (31 - 28*chi2chi2*(-1 + dm) + 28*chi1chi1*(1 + dm) - 56*nchia2 - 112*dm*nchia*nchis - 56*nchis2)*z6) - 12*(-56 + 344*z10 + (-289 - 63*chi2chi2*(-1 + dm) + 63*chi1chi1*(1 + dm) - 378*echia2 - 756*dm*echia*echis - 378*echis2 + 294*lchia - 378*lchia2 + 294*dm*lchis - 756*dm*lchia*lchis - 378*lchis2 + I*(630)*nchia + I*(1260)*lchia*nchia + I*(1260)*dm*lchis*nchia + 
        378*nchia2 + I*(630)*dm*nchis + I*(1260)*dm*lchia*nchis + I*(1260)*lchis*nchis + 756*dm*nchia*nchis + 378*nchis2)*z12 - 516*z16 + 6*(21*chi2chi2*(-1 + dm) - 21*chi1chi1*(1 + dm) + 2*(5 + 42*echia2 + 84*dm*echia*echis + 42*echis2 - 35*lchia + 21*lchia2 - 35*dm*lchis + 42*dm*lchia*lchis + 21*lchis2 - I*(98)*nchia + I*(21)*lchia*nchia + I*(21)*dm*lchis*nchia - I*(98)*dm*nchis + 
        I*(21)*dm*lchia*nchis + I*(21)*lchis*nchis))*z18 + 172*z22 + z24 + 4*z30 - 28*(-10 + 3*lchia + I*(3)*nchia + 9*nchia2 + 3*dm*(lchis + I*nchis + 6*nchia*nchis) + 9*nchis2)*z6) + I*y*z3*(-420 + (-7541 - 3024*chi2chi2*(-1 + dm) + 3024*chi1chi1*(1 + dm) + 15120*lchia - 3024*lchia2 + 15120*dm*lchis - 6048*dm*lchia*lchis - 3024*lchis2 + I*(9072)*nchia - I*(6048)*lchia*nchia - I*(6048)*dm*lchis*nchia - 
        3024*nchia2 + I*(9072)*dm*nchis - I*(6048)*dm*lchia*nchis - I*(6048)*lchis*nchis - 6048*dm*nchia*nchis - 3024*nchis2)*z12 - 11064*z18 + 381*z24 - 12*(29 + 252*nchia2 + 504*dm*nchia*nchis + 252*nchis2)*z6) - 6*y2*z6*(222 - 84*lchia - 84*dm*lchis - I*(84)*nchia + 672*nchia2 - I*(84)*dm*nchis + 1344*dm*nchia*nchis + 672*nchis2 + 1720*z10 - 2033*z12 + 303*z18 + 344*z4 - 750*z6 + 1092*lchia*z6 + 
        1092*dm*lchis*z6 + I*(1092)*nchia*z6 - I*(504)*lchia*nchia*z6 - I*(504)*dm*lchis*nchia*z6 - 840*nchia2*z6 + I*(1092)*dm*nchis*z6 - I*(504)*dm*lchia*nchis*z6 - I*(504)*lchis*nchis*z6 - 1680*dm*nchia*nchis*z6 - 840*nchis2*z6 - 84*chi2chi2*(-1 + dm)*(-1 + 5*z6) + 84*chi1chi1*(1 + dm)*(-1 + 5*z6)) + eta2*(I*(3)*y5*z3*(-35 + 43*z12 - 68*z6) + 6*y6*(-28 + 56*z12 - 13*z6) - 6*I*y3*z3*(-70 - 762*z12 + 151*z18 + 
        (-43 + 168*chi1chi1 + 336*chi1chi2 + 168*chi2chi2 - 672*nchis2)*z6) + I*y*z3*(-420 + (-6721 + 6048*chi1chi1 + 12096*chi1chi2 + 6048*chi2chi2 - 12096*lchis2 - I*(24192)*lchis*nchis - 12096*nchis2)*z12 - 3840*z18 + 309*z24 - 12*(-25 + 126*chi1chi1 + 252*chi1chi2 + 126*chi2chi2 - 504*nchis2)*z6) + 6*y4*(140 - 782*z12 + 297*z18 - 3*(25 + 56*chi1chi1 + 112*chi1chi2 + 56*chi2chi2 - 224*nchis2)*z6) +
        6*y2*(-224 + 220*z10 + (598 - 840*chi1chi1 - 1680*chi1chi2 - 840*chi2chi2 + I*(2016)*lchis*nchis + 3360*nchis2)*z12 + 1100*z16 - 605*z18 + 120*z24 - 6*(-71 + 14*chi1chi1 + 28*chi1chi2 + 14*chi2chi2 - 56*nchis2)*z6) + 6*(112 - 440*z10 - 2*(-733 + 378*chi1chi1 + 756*chi1chi2 + 378*chi2chi2 - I*(1008)*lchis*nchis - 1512*nchis2)*z12 - 220*z16 + (-409 + 504*chi1chi1 + 1008*chi1chi2 + 504*chi2chi2 - 
        2016*lchis2 - I*(2016)*lchis*nchis)*z18 + 220*z22 - 384*z24 + 103*z30 + 28*(-16 + 9*chi1chi1 + 18*chi1chi2 + 9*chi2chi2 - 36*nchis2)*z6)) + 2*eta*(3*y6*(-56 + 142*z12 - 11*z6) - I*(3)*y5*z3*(35 + 62*z12 + 53*z6) - I*y*z3*(-420 - (875 + 3024*chi2chi2*(-2 + dm) - 3024*chi1chi1*(2 + dm) + 1008*lchia + 6048*lchia2 + 1008*dm*lchis + 6048*dm*lchia*lchis + 6048*lchis2 + I*(4752)*nchia + 
        I*(12096)*lchia*nchia + I*(6048)*dm*lchis*nchia + 6048*nchia2 - I*(3744)*dm*nchis + I*(6048)*dm*lchia*nchis + I*(12096)*lchis*nchis + 6048*dm*nchia*nchis + 6048*nchis2)*z12 - 10989*z18 + 1050*z24 - 12*(23 + 126*chi1chi2 + 63*chi2chi2 - 63*chi2chi2*dm + 63*chi1chi1*(1 + dm) + 504*nchia2 - 252*dm*nchia*nchis - 252*nchis2)*z6) + 3*y4*(56 - 200*z12 + 501*z18 + 
        3*(-55 - 56*chi2chi2*(-2 + dm) + 56*chi1chi1*(2 + dm) - 224*nchia2 - 224*dm*nchia*nchis - 224*nchis2)*z6) + 3*y2*(224 - 124*z10 - 4*(716 + 210*chi2chi2*(-2 + dm) - 210*chi1chi1*(2 + dm) + 117*lchia + 93*dm*lchis + I*(3)*nchia + I*(504)*lchia*nchia + I*(252)*dm*lchis*nchia + 840*nchia2 + I*(39)*dm*nchis + I*(252)*dm*lchia*nchis + I*(504)*lchis*nchis + 840*dm*nchia*nchis + 840*nchis2)*z12 - 
        620*z16 - 453*z18 + 372*z24 + 4*(-37 + 126*chi1chi2 - 21*chi2chi2 + 21*chi1chi1*(-1 + dm) - 21*chi2chi2*dm + 21*lchia + 21*dm*lchis + I*(21)*nchia + 672*nchia2 + I*(21)*dm*nchis - 84*dm*nchia*nchis - 84*nchis2)*z6) - 3*(224 - 1128*z10 - 4*(-178 + 252*chi1chi2 + 252*chi2chi2 - 189*chi2chi2*dm + 63*chi1chi1*(4 + 3*dm) - 756*echia2 - 501*lchia - 756*lchia2 + 39*dm*lchis + 
        I*(69)*nchia + I*(2520)*lchia*nchia - I*(252)*dm*lchis*nchia + 756*nchia2 + I*(141)*dm*nchis - I*(252)*dm*lchia*nchis - I*(504)*lchis*nchis - 756*dm*nchia*nchis - 756*nchis2)*z12 + 316*z16 + 3*(271 - 168*chi2chi2*(-2 + dm) + 168*chi1chi1*(2 + dm) - 1344*echia2 - 396*lchia - 672*lchia2 + 4*dm*lchis - 672*dm*lchia*lchis - 672*lchis2 + I*(2108)*nchia - I*(672)*lchia*nchia - 
        I*(336)*dm*lchis*nchia + I*(412)*dm*nchis - I*(336)*dm*lchia*nchis - I*(672)*lchis*nchis)*z18 + 124*z22 - 180*z24 + 43*z30 + 84*(-11 + 6*chi1chi2 + 3*chi2chi2 - 3*chi2chi2*dm + 3*chi1chi1*(1 + dm) + 2*lchia + 2*dm*lchis + 2*I*nchia + 24*nchia2 + 2*I*dm*nchis - 12*dm*nchia*nchis - 12*nchis2)*z6) - 6*I*y3*(63 + 84*chi2chi2*(-2 + dm) - 84*chi1chi1*(2 + dm) + 
        336*nchia2 + 336*dm*nchia*nchis + 336*nchis2 + 178*z12 - 263*z6)*z9)) + Jn*(I*(3024)*Je*(-(echis*(-1 + 2*eta)*(dm*(-2*nchia*y + lchia*z3) + (-1 + 2*eta)*(2*nchis*y - lchis*z3))) + echia*(2*(-1 + 4*eta)*nchia*y + (1 - 4*eta)*lchia*z3 + dm*(-1 + 2*eta)*(2*nchis*y - lchis*z3)))*(-1 + y2 - 2*I*y*z3 - z6)*z9 + Jn*(I*(-3)*y5*z3*(35 + 245*z12 - 88*z6) + 
        6*y6*(-28 + 11*z12 + 65*z6) - 6*I*y3*z3*(70 - 902*z12 + 283*z18 + (109 - 84*chi2chi2*(-1 + dm) + 84*chi1chi1*(1 + dm) - 168*nchia2 - 336*dm*nchia*nchis - 168*nchis2)*z6) - 18*y4*(28 - 39*z12 + 20*z18 + (31 - 28*chi2chi2*(-1 + dm) + 28*chi1chi1*(1 + dm) - 56*nchia2 - 112*dm*nchia*nchis - 56*nchis2)*z6) - 12*(-56 + 344*z10 + (-289 - 63*chi2chi2*(-1 + dm) + 
        63*chi1chi1*(1 + dm) - 378*echia2 - 756*dm*echia*echis - 378*echis2 + 294*lchia - 378*lchia2 + 294*dm*lchis - 756*dm*lchia*lchis - 378*lchis2 + I*(630)*nchia + I*(1260)*lchia*nchia + I*(1260)*dm*lchis*nchia + 378*nchia2 + I*(630)*dm*nchis + I*(1260)*dm*lchia*nchis + I*(1260)*lchis*nchis + 756*dm*nchia*nchis + 378*nchis2)*z12 - 516*z16 + 6*(21*chi2chi2*(-1 + dm) - 
        21*chi1chi1*(1 + dm) + 2*(5 + 42*echia2 + 84*dm*echia*echis + 42*echis2 - 35*lchia + 21*lchia2 - 35*dm*lchis + 42*dm*lchia*lchis + 21*lchis2 - I*(98)*nchia + I*(21)*lchia*nchia + I*(21)*dm*lchis*nchia - I*(98)*dm*nchis + I*(21)*dm*lchia*nchis + I*(21)*lchis*nchis))*z18 + 172*z22 + z24 + 4*z30 - 28*(-10 + 3*lchia + I*(3)*nchia + 9*nchia2 + 3*dm*(lchis + I*nchis + 6*nchia*nchis) + 9*nchis2)*z6) + 
        I*y*z3*(-420 + (-7541 - 3024*chi2chi2*(-1 + dm) + 3024*chi1chi1*(1 + dm) + 15120*lchia - 3024*lchia2 + 15120*dm*lchis - 6048*dm*lchia*lchis - 3024*lchis2 + I*(9072)*nchia - I*(6048)*lchia*nchia - I*(6048)*dm*lchis*nchia - 3024*nchia2 + I*(9072)*dm*nchis - I*(6048)*dm*lchia*nchis - I*(6048)*lchis*nchis - 6048*dm*nchia*nchis - 3024*nchis2)*z12 - 11064*z18 + 381*z24 - 
        12*(29 + 252*nchia2 + 504*dm*nchia*nchis + 252*nchis2)*z6) - 6*y2*z6*(222 - 84*lchia - 84*dm*lchis - I*(84)*nchia + 672*nchia2 - I*(84)*dm*nchis + 1344*dm*nchia*nchis + 672*nchis2 + 1720*z10 - 2033*z12 + 303*z18 + 344*z4 - 750*z6 + 1092*lchia*z6 + 1092*dm*lchis*z6 + I*(1092)*nchia*z6 - I*(504)*lchia*nchia*z6 - I*(504)*dm*lchis*nchia*z6 - 840*nchia2*z6 + I*(1092)*dm*nchis*z6 - 
        I*(504)*dm*lchia*nchis*z6 - I*(504)*lchis*nchis*z6 - 1680*dm*nchia*nchis*z6 - 840*nchis2*z6 - 84*chi2chi2*(-1 + dm)*(-1 + 5*z6) + 84*chi1chi1*(1 + dm)*(-1 + 5*z6)) + eta2*(I*(3)*y5*z3*(-35 + 43*z12 - 68*z6) + 6*y6*(-28 + 56*z12 - 13*z6) - 6*I*y3*z3*(-70 - 762*z12 + 151*z18 + (-43 + 168*chi1chi1 + 336*chi1chi2 + 168*chi2chi2 - 672*nchis2)*z6) + 
        I*y*z3*(-420 + (-6721 + 6048*chi1chi1 + 12096*chi1chi2 + 6048*chi2chi2 - 12096*lchis2 - I*(24192)*lchis*nchis - 12096*nchis2)*z12 - 3840*z18 + 309*z24 - 12*(-25 + 126*chi1chi1 + 252*chi1chi2 + 126*chi2chi2 - 504*nchis2)*z6) + 6*y4*(140 - 782*z12 + 297*z18 - 3*(25 + 56*chi1chi1 + 112*chi1chi2 + 56*chi2chi2 - 224*nchis2)*z6) + 
        6*y2*(-224 + 220*z10 + (598 - 840*chi1chi1 - 1680*chi1chi2 - 840*chi2chi2 + I*(2016)*lchis*nchis + 3360*nchis2)*z12 + 1100*z16 - 605*z18 + 120*z24 - 6*(-71 + 14*chi1chi1 + 28*chi1chi2 + 14*chi2chi2 - 56*nchis2)*z6) + 6*(112 - 440*z10 - 2*(-733 + 378*chi1chi1 + 756*chi1chi2 + 378*chi2chi2 - I*(1008)*lchis*nchis - 1512*nchis2)*z12 - 220*z16 + (-409 + 504*chi1chi1 + 1008*chi1chi2 + 504*chi2chi2 - 2016*lchis2 - 
        I*(2016)*lchis*nchis)*z18 + 220*z22 - 384*z24 + 103*z30 + 28*(-16 + 9*chi1chi1 + 18*chi1chi2 + 9*chi2chi2 - 36*nchis2)*z6)) + 2*eta*(3*y6*(-56 + 142*z12 - 11*z6) - I*(3)*y5*z3*(35 + 62*z12 + 53*z6) - I*y*z3*(-420 - (875 + 3024*chi2chi2*(-2 + dm) - 3024*chi1chi1*(2 + dm) + 1008*lchia + 6048*lchia2 + 1008*dm*lchis + 6048*dm*lchia*lchis + 6048*lchis2 + I*(4752)*nchia + I*(12096)*lchia*nchia + 
        I*(6048)*dm*lchis*nchia + 6048*nchia2 - I*(3744)*dm*nchis + I*(6048)*dm*lchia*nchis + I*(12096)*lchis*nchis + 6048*dm*nchia*nchis + 6048*nchis2)*z12 - 10989*z18 + 1050*z24 - 12*(23 + 126*chi1chi2 + 63*chi2chi2 - 63*chi2chi2*dm + 63*chi1chi1*(1 + dm) + 504*nchia2 - 252*dm*nchia*nchis - 252*nchis2)*z6) + 3*y4*(56 - 200*z12 + 501*z18 + 
        3*(-55 - 56*chi2chi2*(-2 + dm) + 56*chi1chi1*(2 + dm) - 224*nchia2 - 224*dm*nchia*nchis - 224*nchis2)*z6) + 3*y2*(224 - 124*z10 - 4*(716 + 210*chi2chi2*(-2 + dm) - 210*chi1chi1*(2 + dm) + 117*lchia + 93*dm*lchis + I*(3)*nchia + I*(504)*lchia*nchia + I*(252)*dm*lchis*nchia + 840*nchia2 + I*(39)*dm*nchis + I*(252)*dm*lchia*nchis + I*(504)*lchis*nchis + 840*dm*nchia*nchis + 840*nchis2)*z12 - 
        620*z16 - 453*z18 + 372*z24 + 4*(-37 + 126*chi1chi2 - 21*chi2chi2 + 21*chi1chi1*(-1 + dm) - 21*chi2chi2*dm + 21*lchia + 21*dm*lchis + I*(21)*nchia + 672*nchia2 + I*(21)*dm*nchis - 84*dm*nchia*nchis - 84*nchis2)*z6) - 3*(224 - 1128*z10 - 4*(-178 + 252*chi1chi2 + 252*chi2chi2 - 189*chi2chi2*dm + 63*chi1chi1*(4 + 3*dm) - 756*echia2 - 501*lchia - 756*lchia2 + 39*dm*lchis + 
        I*(69)*nchia + I*(2520)*lchia*nchia - I*(252)*dm*lchis*nchia + 756*nchia2 + I*(141)*dm*nchis - I*(252)*dm*lchia*nchis - I*(504)*lchis*nchis - 756*dm*nchia*nchis - 756*nchis2)*z12 + 316*z16 + 3*(271 - 168*chi2chi2*(-2 + dm) + 168*chi1chi1*(2 + dm) - 1344*echia2 - 396*lchia - 672*lchia2 + 4*dm*lchis - 672*dm*lchia*lchis - 672*lchis2 + I*(2108)*nchia - I*(672)*lchia*nchia - 
        I*(336)*dm*lchis*nchia + I*(412)*dm*nchis - I*(336)*dm*lchia*nchis - I*(672)*lchis*nchis)*z18 + 124*z22 - 180*z24 + 43*z30 + 84*(-11 + 6*chi1chi2 + 3*chi2chi2 - 3*chi2chi2*dm + 3*chi1chi1*(1 + dm) + 2*lchia + 2*dm*lchis + 2*I*nchia + 24*nchia2 + 2*I*dm*nchis - 12*dm*nchia*nchis - 12*nchis2)*z6) - 6*I*y3*(63 + 84*chi2chi2*(-2 + dm) - 84*chi1chi1*(2 + dm) + 
        336*nchia2 + 336*dm*nchia*nchis + 336*nchis2 + 178*z12 - 263*z6)*z9)))))/(6048*(Jl2 + Jn2)*z14);
    return fp0 + fp2 + fp3 + fp4;
}

static COMPLEX16 prec_CalculateWaveformFactor_mode21(SEOBPrecWaveformVariables *vars,
    REAL8 eta, REAL8 dm)
{
    REAL8 eta2 = eta*eta;
    REAL8 nchia, nchis, lchia, lchis, echia, echis;
    REAL8 nchia2, nchis2, lchia2, lchis2, echia2, echis2;
    REAL8 chi1chi1, chi1chi2, chi2chi2;
    REAL8 Jn, Jl, Je;
    REAL8 Jn2, Jl2, Je2;
    REAL8 x, x2, x3, x4;
    REAL8 y, y2, y3, y4, y5, y6;
    REAL8 z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13, z14, z15, z16, z18, z20, z21, z22, z24, z30;
    nchia = vars->nchia;
    nchis = vars->nchis;
    lchia = vars->lchia;
    lchis = vars->lchis;
    echia = vars->echia;
    echis = vars->echis;
    nchia2 = nchia*nchia;
    nchis2 = nchis*nchis;
    lchia2 = lchia*lchia;
    lchis2 = lchis*lchis;
    echia2 = echia*echia;
    echis2 = echis*echis;
    chi1chi1 = vars->chi1chi1;
    chi1chi2 = vars->chi1chi2;
    chi2chi2 = vars->chi2chi2;
    Jn = vars->Jn;
    Jl = vars->Jl;
    Je = vars->Je;
    Jn2 = Jn*Jn;
    Jl2 = Jl*Jl;
    Je2 = Je*Je;
    x = vars->x;x2 = vars->x2;x3 = vars->x3;x4 = vars->x4;
    y = vars->y;y2 = vars->y2;y3 = vars->y3;y4 = vars->y4;y5 = vars->y5;y6 = vars->y6;
    z2 = vars->z2;z3 = vars->z3;z4 = vars->z4;z5 = vars->z5;z6 = vars->z6;
    z7 = vars->z7;z8 = vars->z8;z9 = vars->z9;z10 = vars->z10;z11 = vars->z11;
    z12 = vars->z12;z13 = vars->z13;z14 = vars->z14;z15 = vars->z15;z16 = vars->z16;
    z18 = vars->z18;z20 = vars->z20;z21 = vars->z21;z22 = vars->z22;z24 = vars->z24;
    z30 = vars->z30;

    COMPLEX16 fp1, fp2, fp3, fp4_1, fp4_2, fp4_3;
    fp1 = dm*(-1. + 1./z4);
    fp2 = (3*(echia + dm*echis)*x*(-1 + z8))/(2*z7);
    fp3 = (x2*(14*dm*(2 + eta*(-2 + y2) + y2) + 2*dm*(55 - 12*eta)*z12 - 
        I*(42)*(3*dm*lchia + (3 - 8*eta)*lchis)*y*z3 + (168*eta*lchis + 2*dm*(-69 + 26*eta - I*(126)*nchia) - 
        I*(42)*(6 + eta)*nchis + dm*(-64 + 23*eta)*y2)*z6 - I*(3)*dm*(-83 + 12*eta)*y*z9))/(42*z10);
    fp4_1 = (dm*echis*x3*(-147*(2 + eta*(-2 + y2) + y2) + 15*(-35 + 33*eta)*z12 - 21*(-3 - eta + 5*(1 + eta)*y2)*z14 - 
        2*(427 + 79*eta)*z16 - 21*(1 + eta)*z20 + (-(eta*(673 + 45*y2)) + 7*(239 + I*(216)*nchia + 51*y2))*z6 - 
        21*(2 + eta*(-2 + y2) + y2)*z8 - I*(36)*(49 + 2*eta)*y*z9))/(168*z13);
    fp4_2 = (x3*(I*(1512)*echis*nchis*z6 - echia*(147*(2 + eta*(-2 + y2) + y2) + 
        15*(35 - 121*eta)*z12 + 21*(-3 - eta + 5*(1 + eta)*y2)*z14 + 2*(427 + 
        383*eta)*z16 + 21*(1 + eta)*z20 + (-1673 + 1385*eta - I*(1512)*nchia 
        + I*(6048)*eta*nchia - 3*(119 + 57*eta)*y2)*z6 + 21*(2 + eta*(-2 + 
        y2) + y2)*z8 - I*(12)*(-147 + 104*eta)*y*z9)))/(168*z13);
    fp4_3 = (I*(1)*dm*x3*(9*echia*(Jl2 + Jn2)*nchis + eta*Je*(Jn*lchis*(1 + y2) + (Jl*lchis - Jn*nchis)*y*z3 - Jl*nchis*z6)))/((Jl2 + Jn2)*z7);
    return fp1 + fp2 + fp3 + fp4_1 + fp4_2 + fp4_3;
}

static COMPLEX16 prec_CalculateWaveformFactor_mode33(SEOBPrecWaveformVariables *vars,
    REAL8 eta, REAL8 dm)
{
    REAL8 eta2 = eta*eta;
    REAL8 nchia, nchis, lchia, lchis, echia, echis;
    REAL8 nchia2, nchis2, lchia2, lchis2, echia2, echis2;
    REAL8 chi1chi1, chi1chi2, chi2chi2;
    REAL8 Jn, Jl, Je;
    REAL8 Jn2, Jl2, Je2;
    REAL8 x, x2, x3, x4;
    REAL8 y, y2, y3, y4, y5, y6;
    REAL8 z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13, z14, z15, z16, z18, z20, z21, z22, z24, z30;
    nchia = vars->nchia;
    nchis = vars->nchis;
    lchia = vars->lchia;
    lchis = vars->lchis;
    echia = vars->echia;
    echis = vars->echis;
    nchia2 = nchia*nchia;
    nchis2 = nchis*nchis;
    lchia2 = lchia*lchia;
    lchis2 = lchis*lchis;
    echia2 = echia*echia;
    echis2 = echis*echis;
    chi1chi1 = vars->chi1chi1;
    chi1chi2 = vars->chi1chi2;
    chi2chi2 = vars->chi2chi2;
    Jn = vars->Jn;
    Jl = vars->Jl;
    Je = vars->Je;
    Jn2 = Jn*Jn;
    Jl2 = Jl*Jl;
    Je2 = Je*Je;
    x = vars->x;x2 = vars->x2;x3 = vars->x3;x4 = vars->x4;
    y = vars->y;y2 = vars->y2;y3 = vars->y3;y4 = vars->y4;y5 = vars->y5;y6 = vars->y6;
    z2 = vars->z2;z3 = vars->z3;z4 = vars->z4;z5 = vars->z5;z6 = vars->z6;
    z7 = vars->z7;z8 = vars->z8;z9 = vars->z9;z10 = vars->z10;z11 = vars->z11;
    z12 = vars->z12;z13 = vars->z13;z14 = vars->z14;z15 = vars->z15;z16 = vars->z16;
    z18 = vars->z18;z20 = vars->z20;z21 = vars->z21;z22 = vars->z22;z24 = vars->z24;
    z30 = vars->z30;

    COMPLEX16 fp1, fp3, fp4;
    fp1 = dm*(2*(I*(2)*y - I*(1)*y3 - z3 - 3*y2*z3 + I*(3)*y*z6 + z9))/(9*z3);
    fp3 = (x2*(I*(36)*dm*y - I*(36)*dm*eta*y + I*(36)*dm*eta*y3 - I*(9)*dm*y5 - I*(9)*dm*eta*y5 + 
        567*dm*z11 - 324*dm*eta*z11 + I*(63)*dm*y*z12 - I*(207)*dm*eta*y*z12 + I*(96)*dm*y3*z12 + 
        I*(168)*dm*eta*y3*z12 - 45*dm*z15 + 126*dm*eta*z15 + 102*dm*y2*z15 + 30*dm*eta*y2*z15 - I*(45)*dm*y*z18 + 
        I*(63)*dm*eta*y*z18 - 6*dm*z21 + 30*dm*eta*z21 - I*(254)*dm*y*z6 + I*(100)*dm*eta*y*z6 - I*(72)*eta*lchis*y*z6 + 
        72*eta*nchis*y*z6 - I*(6)*dm*y3*z6 + I*(48)*dm*eta*y3*z6 - I*(3)*dm*y5*z6 - I*(39)*dm*eta*y5*z6 - 516*dm*z9 + 
        168*dm*eta*z9 - 288*eta*lchis*z9 - I*(288)*eta*nchis*z9 - 36*dm*y2*z9 + 234*dm*eta*y2*z9 - 36*dm*y4*z9 - 144*dm*eta*y4*z9))/(162*z9);
    fp4 = (x3*(echia*(Jl2 + Jn2)*(I*(2)*(-1 + 4*eta)*y*(-2 + y2) + (-2 + 5*eta*(2 - 5*y2) + 6*y2)*z3 + I*(1)*(-24 + 101*eta)*y*z6 + 2*(1 - 5*eta)*z9) + 
        dm*(3*eta*Je*(Jn*lchis*(1 + y2) + (Jl*lchis - Jn*nchis)*y*z3 - Jl*nchis*z6)*(2*y*(-2 + y2) + I*(1)*(7 - 6*y2)*z3 - 6*y*z6 + I*(2)*z9) + 
        echis*(Jl2 + Jn2)*(I*(1)*(-2 + eta)*y*(-2 + y2) + (-2 + 3*eta + 6*y2 - 4*eta*y2)*z3 + I*(1)*(-24 + 17*eta)*y*z6 + (2 - 3*eta)*z9))))/(9*(Jl2 + Jn2)*z6);
    return fp1 + fp3 + fp4;
}

static COMPLEX16 prec_CalculateWaveformFactor_mode44(SEOBPrecWaveformVariables *vars,
    REAL8 eta, REAL8 dm)
{
    REAL8 eta2 = eta*eta;
    REAL8 nchia, nchis, lchia, lchis, echia, echis;
    REAL8 nchia2, nchis2, lchia2, lchis2, echia2, echis2;
    REAL8 chi1chi1, chi1chi2, chi2chi2;
    REAL8 Jn, Jl, Je;
    REAL8 Jn2, Jl2, Je2;
    REAL8 x, x2, x3, x4;
    REAL8 y, y2, y3, y4, y5, y6;
    REAL8 z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13, z14, z15, z16, z18, z20, z21, z22, z24, z30;
    nchia = vars->nchia;
    nchis = vars->nchis;
    lchia = vars->lchia;
    lchis = vars->lchis;
    echia = vars->echia;
    echis = vars->echis;
    nchia2 = nchia*nchia;
    nchis2 = nchis*nchis;
    lchia2 = lchia*lchia;
    lchis2 = lchis*lchis;
    echia2 = echia*echia;
    echis2 = echis*echis;
    Jn = vars->Jn;
    Jl = vars->Jl;
    Je = vars->Je;
    Jn2 = Jn*Jn;
    Jl2 = Jl*Jl;
    Je2 = Je*Je;
    x = vars->x;x2 = vars->x2;x3 = vars->x3;x4 = vars->x4;
    y = vars->y;y2 = vars->y2;y3 = vars->y3;y4 = vars->y4;y5 = vars->y5;y6 = vars->y6;
    z2 = vars->z2;z3 = vars->z3;z4 = vars->z4;z5 = vars->z5;z6 = vars->z6;
    z7 = vars->z7;z8 = vars->z8;z9 = vars->z9;z10 = vars->z10;z11 = vars->z11;
    z12 = vars->z12;z13 = vars->z13;z14 = vars->z14;z15 = vars->z15;z16 = vars->z16;
    z18 = vars->z18;z20 = vars->z20;z21 = vars->z21;z22 = vars->z22;z24 = vars->z24;
    z30 = vars->z30;

    COMPLEX16 fp2, fp4;
    fp2 = (7 + 6*y4 + 6*z12 - I*(24)*y3*z3 - 64*z4 + 51*z6 - 18*y2*(1 + 2*z6) + I*(6)*y*z3*(9 + 4*z6))/(64*z4);
    fp4 = (x2*(220*(-1 + 3*eta)*(2 + eta*(-2 + y2) + y2)*(7 + 6*(-3 + y2)*y2) + (-57943 + 17478*y2 + 
        60*(eta*(4351 - 3042*eta - 759*lchia - I*(759)*nchia + 759*dm*(lchis + I*(1)*nchis)) + 
        3*eta*(-898 + 1503*eta)*y2 + (356 - eta*(262 + 2103*eta))*y4))*z12 + 
        I*(6)*y*(-1647 + 10930*eta - 22755*eta2 + 4840*(1 + eta)*(-1 + 3*eta)*y2)*z15 - 
        3*(-6961 + 60*(713 - 636*eta)*eta + 6520*y2 + 380*eta*(-58 + 3*eta)*y2)*z18 + I*(120)*(49 + eta*(-278 + 267*eta))*y*z21 + 
        60*(8 + eta*(-106 + 183*eta))*z24 - I*(330)*(-1 + 3*eta)*y*(-3 + 2*y)*(3 + 2*y)*(2 + eta*(-2 + y2) + y2)*z3 - 
        20*(-1983 + eta*(6934 - 3315*eta + 396*lchia + I*(396)*nchia - 396*dm*(lchis + I*(1)*nchis)) + 
        1414*y2 + 2*eta*(-2347 + 1335*eta - 99*lchia - I*(99)*nchia + 99*dm*(lchis + I*(1)*nchis))*y2 + 
        3*(18 + eta*(-310 + 453*eta))*y4 + 3*(14 + (62 - 249*eta)*eta)*y6)*z6 + I*(6)*y*(15289 + 1163*y2 + 
        5*(eta*(-9448 + 2727*eta - 660*lchia - I*(660)*nchia + 660*dm*(lchis + I*(1)*nchis)) + eta*(-3154 + 5277*eta)*y2 + 
        4*(61 + (58 - 597*eta)*eta)*y4))*z9))/(42240*(-1 + 3*eta)*z10);
    return fp2 + fp4;
}

static COMPLEX16 prec_CalculateWaveformFactor_mode55(SEOBPrecWaveformVariables *vars,
    REAL8 eta, REAL8 dm)
{
    REAL8 eta2 = eta*eta;
    REAL8 nchia, nchis, lchia, lchis, echia, echis;
    REAL8 nchia2, nchis2, lchia2, lchis2, echia2, echis2;
    REAL8 chi1chi1, chi1chi2, chi2chi2;
    REAL8 Jn, Jl, Je;
    REAL8 Jn2, Jl2, Je2;
    REAL8 x, x2, x3, x4;
    REAL8 y, y2, y3, y4, y5, y6;
    REAL8 z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13, z14, z15, z16, z18, z20, z21, z22, z24, z30;
    nchia = vars->nchia;
    nchis = vars->nchis;
    lchia = vars->lchia;
    lchis = vars->lchis;
    echia = vars->echia;
    echis = vars->echis;
    nchia2 = nchia*nchia;
    nchis2 = nchis*nchis;
    lchia2 = lchia*lchia;
    lchis2 = lchis*lchis;
    echia2 = echia*echia;
    echis2 = echis*echis;
    Jn = vars->Jn;
    Jl = vars->Jl;
    Je = vars->Je;
    Jn2 = Jn*Jn;
    Jl2 = Jl*Jl;
    Je2 = Je*Je;
    x = vars->x;x2 = vars->x2;x3 = vars->x3;x4 = vars->x4;
    y = vars->y;y2 = vars->y2;y3 = vars->y3;y4 = vars->y4;y5 = vars->y5;y6 = vars->y6;
    z2 = vars->z2;z3 = vars->z3;z4 = vars->z4;z5 = vars->z5;z6 = vars->z6;
    z7 = vars->z7;z8 = vars->z8;z9 = vars->z9;z10 = vars->z10;z11 = vars->z11;
    z12 = vars->z12;z13 = vars->z13;z14 = vars->z14;z15 = vars->z15;z16 = vars->z16;
    z18 = vars->z18;z20 = vars->z20;z21 = vars->z21;z22 = vars->z22;z24 = vars->z24;
    z30 = vars->z30;
    COMPLEX16 fac;
    fac = (I*(2)*y*(41 - 48*y2 + 12*y4) + I*(120)*y*z12 + 24*z15 + 4*(43 - 99*y2 + 30*y4)*z3 - 625*z5 - I*(48)*y*(-13 + 5*y2)*z6 + 3*(143 - 80*y2)*z9)/(625*z5);
    return fac;
}


/**
 * This function calculates hlm mode factorized-resummed waveform
 * for given dynamical variables.
 * Eq. 17 and the entire Appendix of PRD 86, 024011 (2012) + changes
 * described in ￼the section "Factorized waveforms" of https://dcc.ligo.org/T1400476
 */
INT
prec_EOBGetPrecEccSpinFactorizedWaveform_v1(
            COMPLEX16 * hlm,	/**< OUTPUT, hlm waveforms */
            REAL8Vector * values,	/**< dyanmical variables: \f$(r,\phi,p_r,p_\phi)\f$ */
            REAL8Vector * cartvalues,	/**< dyanmical variables */
            const REAL8 v,	/**< velocity */
            const REAL8 Hreal,	/**< real Hamiltonian */
            const INT l,	/**< l mode index */
            const INT m,	/**< m mode index */
            SEOBPrecWaveformVariables *vars,
            SpinEOBParams * params	/**< Spin EOB parameters */
)
{
    int		debugPK = 0;
    /* Status of function calls */
    INT status;
    INT i;

    REAL8 eta;
    REAL8 r , pp, Omega, v2, Omegav2, vh, vh3, k, hathatk, eulerlogxabs;
    //pr
    REAL8 rcrossp_x, rcrossp_y, rcrossp_z;
    REAL8 Slm, deltalm;
    COMPLEX16 auxflm = 0.0;
    REAL8 hathatksq4, hathatk4pi, Tlmprefac, Tlmprodfac;
    COMPLEX16 Tlm, rholmPwrl,rholm;
    COMPLEX16 hNewton;
    COMPLEX16 facEcc;
    gsl_sf_result	z2;

    /* Non-Keplerian velocity */
    REAL8		vPhi    , vPhi2;

    /* Pre-computed coefficients */

    FacWaveformCoeffs *hCoeffs = params->hCoeffs;

    if (abs(m) > (INT) l) 
    {
        return CEV_FAILURE;
    }
    if (m == 0) 
    {
        return CEV_FAILURE;
    }
    eta = params->eta;

    /* Check our eta was sensible */
    if (eta > 0.25) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Eta seems to be > 0.25 - this isn't allowed!");
        return CEV_FAILURE;
    }
    /*
    * else if ( eta == 0.25 && m % 2 ) { // If m is odd and dM = 0, hLM
    * will be zero memset( hlm, 0, sizeof( COMPLEX16 ) ); return
    * XLAL_SUCCESS; }
    */

    r = values->data[0];
    pp = values->data[3];

    rcrossp_x = cartvalues->data[1] * cartvalues->data[5] - cartvalues->data[2] * cartvalues->data[4];
    rcrossp_y = cartvalues->data[2] * cartvalues->data[3] - cartvalues->data[0] * cartvalues->data[5];
    rcrossp_z = cartvalues->data[0] * cartvalues->data[4] - cartvalues->data[1] * cartvalues->data[3];

    //pp = sqrt(rcrossp_x * rcrossp_x + rcrossp_y * rcrossp_y + rcrossp_z * rcrossp_z);

    v2 = v * v;
    Omega = v2 * v;
    Omegav2 = Omega*v2;
    vh3 = Hreal * Omega;
    vh = cbrt(vh3);
    eulerlogxabs = CST_GAMMA + log(2.0 * (REAL8) m * v);

    FacWaveformCoeffs oldCoeffs = *(params->hCoeffs); //RC: the function XLALSimIMRSpinPrecEOBNonKeplerCoeff is calculating again the coefficients hCoeffs, for the omega. These are different because
    // for the dynamics we are using different coefficients for the 21 mode. I store the old coefficients and the I recover them ofter vPhi is computed.

    /* Calculate the non-Keplerian velocity */
    if (params->alignedSpins) 
    {

        vPhi = XLALSimIMRSpinAlignedEOBNonKeplerCoeff(values->data, params);

        if (IS_REAL8_FAIL_NAN(vPhi)) 
        {
            return CEV_FAILURE;
        }
        vPhi = r * cbrt(vPhi);

        // if (debugPK)
        //     XLAL_PRINT_INFO("In XLALSimIMRSpinEOBFluxCalculateNewtonianMultipole, getting rW = %.12e\n", vPhi);
        vPhi *= Omega;
        vPhi2 = vPhi * vPhi;
    } 
    else 
    {
        vPhi = XLALSimIMRSpinPrecEOBNonKeplerCoeff(cartvalues->data, params);
        if (IS_REAL8_FAIL_NAN(vPhi)) 
        {
            return CEV_FAILURE;
        }
        vPhi = r * cbrt(vPhi);

        // if (debugPK)
        //     XLAL_PRINT_INFO("In XLALSimIMRSpinEOBFluxCalculateNewtonianMultipole, getting rW = %.12e\n",
        //             vPhi);
        vPhi *= Omega;
        vPhi2 = vPhi * vPhi;
    }
    *hCoeffs = oldCoeffs; //RC: Here I recover the old coefficients
    /*
        * Calculate the newtonian multipole, 1st term in Eq. 17, given by
        * Eq. A1
        */
    // if (debugPK) {
    //     XLAL_PRINT_INFO("\nValues inside XLALSimIMRSpinEOBFluxGetSpinFactorizedWaveform:\n");
    //     for (i = 0; i < 11; i++)
    //         XLAL_PRINT_INFO("values[%d] = %.12e\n", i, cartvalues->data[i]);

    //     XLAL_PRINT_INFO("Calculating hNewton, with v = %.12e, vPhi = %.12e, r = %.12e, Phi = %.12e, l = %d, m = %d\n",
    //             v, vPhi, r, values->data[1], (UINT4) l, (UINT4) m);
    // }
    status = XLALSimIMRSpinEOBCalculateNewtonianMultipole(&hNewton,
            vPhi2, r, cartvalues->data[12] + cartvalues->data[13], (UINT) l, m, params);
                                
    if (status!=CEV_SUCCESS) 
    {
        return CEV_FAILURE;
    }
    /* Calculate the source term, 2nd term in Eq. 17, given by Eq. A5, Hreal is given by Eq.5 and Heff is in Eq.2 */
    if (((l + m) % 2) == 0) 
    {
        Slm = (Hreal * Hreal - 1.) / (2. * eta) + 1.;
    } 
    else 
    {
        Slm = v * pp;
        //Slm = v * sqrt(rcrossp_x * rcrossp_x + rcrossp_y * rcrossp_y + rcrossp_z * rcrossp_z);
    }
    // if (debugPK)
    //     XLAL_PRINT_INFO("In XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform: Hreal = %e, Slm = %e, eta = %e\n", Hreal, Slm, eta);

    /*
    * Calculate the Tail term, 3rd term in Eq. 17,
    * given by Eq. A6, and Eq. (42) of
    * http://arxiv.org/pdf/1212.4357.pdf (or PRD 87 084035 (2013))
    */
    k = m * Omega;
    hathatk = Hreal * k;

    gsl_sf_result	lnr1, arg1;
    status = gsl_sf_lngamma_complex_e(l + 1.0, -2.0 * hathatk, &lnr1, &arg1);
    if (status != GSL_SUCCESS) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Error in GSL function");
        return CEV_FAILURE;
    }
    status = gsl_sf_fact_e(l, &z2);
    if (status != GSL_SUCCESS) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Error in GSL function");
        return CEV_FAILURE;
    }
    Tlm = cexp((lnr1.val + CST_PI * hathatk) + I * (
            arg1.val + 2.0 * hathatk * log(4.0 * k / sqrt(CST_E))));
    Tlm /= z2.val;


    // if (debugPK){
    //     hathatksq4 = 4. * hathatk * hathatk;
    //     hathatk4pi = 4. * LAL_PI * hathatk;
    //     /* Calculating the prefactor of Tlm, outside the multiple product */
    //     Tlmprefac = sqrt(hathatk4pi / (1. - exp(-hathatk4pi))) / z2.val;

    //     /* Calculating the multiple product factor */
    //     for (Tlmprodfac = 1., i = 1; i <= l; i++)
    //         Tlmprodfac *= (hathatksq4 + (REAL8) i * i);

    //     REAL8		Tlmold;
    //     Tlmold = Tlmprefac * sqrt(Tlmprodfac);
    //     XLAL_PRINT_INFO("Tlm = %e + i%e, |Tlm| = %.16e (should be %.16e)\n", creal(Tlm), cimag(Tlm), cabs(Tlm), Tlmold);
    // }

    /* Calculate the residue phase and amplitude terms */
    /*
        * deltalm is the 4th term in Eq. 17, delta 22 given by Eq. A15,
        * others
        */
    /* rholm is the 5th term in Eq. 17, given by Eqs. A8 - A14 */
    /*
        * auxflm is a special part of the 5th term in Eq. 17, given by Eq.
        * A15
        */
    /*
        * Actual values of the coefficients are defined in the another function
        * see file LALSimIMRSpinEOBFactorizedWaveformCoefficientsPrec.c
        */
    facEcc = 0.0;
    REAL8 dm;
    if (1. - 4.*eta < 0.0) dm = 0.0;
    else dm = sqrt(1. - 4.*eta);
    switch (l) {
    case 2:
        switch (abs(m)) {
        case 2:
            deltalm = vh3 * (hCoeffs->delta22vh3 + vh3 * (hCoeffs->delta22vh6
                    + vh * vh * (hCoeffs->delta22vh9 * vh)))
                + Omega*(hCoeffs->delta22v5 * v2 + Omega*(hCoeffs->delta22v6  + hCoeffs->delta22v8 *v2));
            rholm = 1. + v2 * (hCoeffs->rho22v2 + v * (hCoeffs->rho22v3
                                + v * (hCoeffs->rho22v4
                    + v * (hCoeffs->rho22v5 + v * (hCoeffs->rho22v6
                                + hCoeffs->rho22v6l * eulerlogxabs + v * (hCoeffs->rho22v7
                                                        + v * (hCoeffs->rho22v8 + hCoeffs->rho22v8l * eulerlogxabs
                                                            + (hCoeffs->rho22v10 + hCoeffs->rho22v10l * eulerlogxabs) * v2)))))));
            facEcc = prec_CalculateWaveformFactor_mode22_v1(vars, eta, dm);
            //FIXME
                // if (debugPK){
                //     XLAL_PRINT_INFO("PK:: rho22v2 = %.12e, rho22v3 = %.12e, rho22v4 = %.12e,\n rho22v5 = %.16e, rho22v6 = %.16e, rho22v6LOG = %.16e, \n rho22v7 = %.12e, rho22v8 = %.16e, rho22v8LOG = %.16e, \n rho22v10 = %.16e, rho22v10LOG = %.16e\n, rho22v6 = %.12e, rho22v8 = %.12e, rho22v10 = %.12e\n",
                //         hCoeffs->rho22v2, hCoeffs->rho22v3, hCoeffs->rho22v4,
                //         hCoeffs->rho22v5, hCoeffs->rho22v6, hCoeffs->rho22v6l,
                //         hCoeffs->rho22v7, hCoeffs->rho22v8, hCoeffs->rho22v8l,
                //         hCoeffs->rho22v10, hCoeffs->rho22v10l,
                //         hCoeffs->rho22v6 + hCoeffs->rho22v6l * eulerlogxabs,
                //         hCoeffs->rho22v8 + hCoeffs->rho22v8l * eulerlogxabs,
                //         hCoeffs->rho22v10 + hCoeffs->rho22v10l * eulerlogxabs);}
            break;
        case 1:
            {
                deltalm = vh3 * (hCoeffs->delta21vh3 + vh3 * (hCoeffs->delta21vh6
                                            + vh * (hCoeffs->delta21vh7 + (hCoeffs->delta21vh9) * vh * vh)))
                    + Omegav2*(hCoeffs->delta21v5  + hCoeffs->delta21v7 * v2);
                rholm = 1. + v * (hCoeffs->rho21v1
                            + v * (hCoeffs->rho21v2 + v * (hCoeffs->rho21v3 + v * (hCoeffs->rho21v4
                                                        + v * (hCoeffs->rho21v5 + v * (hCoeffs->rho21v6 + hCoeffs->rho21v6l * eulerlogxabs
                                                                    + v * (hCoeffs->rho21v7 + hCoeffs->rho21v7l * eulerlogxabs
                                                                            + v * (hCoeffs->rho21v8 + hCoeffs->rho21v8l * eulerlogxabs
                                                                                + (hCoeffs->rho21v10 + hCoeffs->rho21v10l * eulerlogxabs) * v2))))))));
                auxflm = v * (hCoeffs->f21v1 + v2 * (hCoeffs->f21v3 + v * hCoeffs->f21v4 + v2 * (hCoeffs->f21v5 + v * hCoeffs->f21v6 + v2 * hCoeffs->f21v7c)));
                facEcc = prec_CalculateWaveformFactor_mode21(vars, eta, dm);
            }
            break;
        default:
            return CEV_FAILURE;
            break;
        }
        break;
    case 3:
        switch (m) {
        case 3:
            deltalm = vh3 * (hCoeffs->delta33vh3 + vh3 * (hCoeffs->delta33vh6 + hCoeffs->delta33vh9 * vh3))
                + hCoeffs->delta33v5 * v * v2 * v2 + hCoeffs->delta33v7 * v2 * v2 * v2 * v;	
            //R.C: delta33v7 is set to 0, whoever is adding it here as a coefficient is evil, TODO: double check that is zero and then remove it
            //RC: This terms are in Eq.A6 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
            rholm = 1. + v2 * (hCoeffs->rho33v2 + v * (hCoeffs->rho33v3 + v * (hCoeffs->rho33v4 + v * (hCoeffs->rho33v5 + v * (hCoeffs->rho33v6 +
                    hCoeffs->rho33v6l * eulerlogxabs + v * (hCoeffs->rho33v7 + v * (hCoeffs->rho33v8 + hCoeffs->rho33v8l * eulerlogxabs +
                    v2*(hCoeffs->rho33v10 + hCoeffs->rho33v10l*eulerlogxabs))))))));
            //RC: This terms are in Eq.A10 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
            auxflm = v * (v2 * (hCoeffs->f33v3 + v * (hCoeffs->f33v4 + v * (hCoeffs->f33v5  + v * hCoeffs->f33v6)))) + _Complex_I * vh3 * vh3 * hCoeffs->f33vh6;
            facEcc = prec_CalculateWaveformFactor_mode33(vars, eta, dm);
            break;
        case 2:
            deltalm = vh3 * (hCoeffs->delta32vh3 + vh * (hCoeffs->delta32vh4 + vh * vh * (hCoeffs->delta32vh6
                            + hCoeffs->delta32vh9 * vh3)));
            rholm = 1. + v * (hCoeffs->rho32v
                        + v * (hCoeffs->rho32v2 + v * (hCoeffs->rho32v3 + v * (hCoeffs->rho32v4 + v * (hCoeffs->rho32v5
                                                                + v * (hCoeffs->rho32v6 + hCoeffs->rho32v6l * eulerlogxabs
                                                                + (hCoeffs->rho32v8 + hCoeffs->rho32v8l * eulerlogxabs) * v2))))));
            break;
        case 1:
            deltalm = vh3 * (hCoeffs->delta31vh3 + vh3 * (hCoeffs->delta31vh6
                                        + vh * (hCoeffs->delta31vh7 + hCoeffs->delta31vh9 * vh * vh)))
                + hCoeffs->delta31v5 * v * v2 * v2;
            rholm = 1. + v2 * (hCoeffs->rho31v2 + v * (hCoeffs->rho31v3 + v * (hCoeffs->rho31v4
                                                + v * (hCoeffs->rho31v5 + v * (hCoeffs->rho31v6 + hCoeffs->rho31v6l * eulerlogxabs
                                                                + v * (hCoeffs->rho31v7 + (hCoeffs->rho31v8 + hCoeffs->rho31v8l * eulerlogxabs) * v))))));
            auxflm = v * v2 * hCoeffs->f31v3;
            break;
        default:
            return CEV_FAILURE;
            break;
        }
        break;
    case 4:
        switch (m) {
        case 4:
            //RC: This terms are in Eq.A15 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
            deltalm = vh3 * (hCoeffs->delta44vh3 + vh3 * (hCoeffs->delta44vh6 + vh3 * hCoeffs->delta44vh9))
                    + hCoeffs->delta44v5 * v2 * v2 * v;
            //RC: This terms are in Eq.A8 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
            rholm = 1. + v2 * (hCoeffs->rho44v2
                    + v * (hCoeffs->rho44v3 + v * (hCoeffs->rho44v4 + v * (hCoeffs->rho44v5 + v * (hCoeffs->rho44v6 + hCoeffs->rho44v6l *
                    eulerlogxabs + v2 *( hCoeffs->rho44v8 + hCoeffs->rho44v8l * eulerlogxabs +v2 * (hCoeffs->rho44v10 + hCoeffs->rho44v10l * eulerlogxabs) ) )))));
            facEcc = prec_CalculateWaveformFactor_mode44(vars, eta, dm);
            break;
        case 3:
            deltalm = vh3 * (hCoeffs->delta43vh3 + vh * (hCoeffs->delta43vh4
                        + hCoeffs->delta43vh6 * vh * vh));
            rholm = 1. + v * (hCoeffs->rho43v
                        + v * (hCoeffs->rho43v2
                + v2 * (hCoeffs->rho43v4 + v * (hCoeffs->rho43v5
                                + (hCoeffs->rho43v6 + hCoeffs->rho43v6l * eulerlogxabs) * v))));
            auxflm = v * hCoeffs->f43v;
            break;
        case 2:
            deltalm = vh3 * (hCoeffs->delta42vh3 + hCoeffs->delta42vh6 * vh3);
            rholm = 1. + v2 * (hCoeffs->rho42v2
                        + v * (hCoeffs->rho42v3 + v * (hCoeffs->rho42v4 + v * (hCoeffs->rho42v5
                                                    + (hCoeffs->rho42v6 + hCoeffs->rho42v6l * eulerlogxabs) * v))));
            break;
        case 1:
            deltalm = vh3 * (hCoeffs->delta41vh3 + vh * (hCoeffs->delta41vh4
                        + hCoeffs->delta41vh6 * vh * vh));
            rholm = 1. + v * (hCoeffs->rho41v
                        + v * (hCoeffs->rho41v2
                + v2 * (hCoeffs->rho41v4 + v * (hCoeffs->rho41v5
                                + (hCoeffs->rho41v6 + hCoeffs->rho41v6l * eulerlogxabs) * v))));
            auxflm = v * hCoeffs->f41v;
            break;
        default:
            return CEV_FAILURE;
            break;
        }
        break;
    case 5:
        switch (m) {
        case 5:
            //RC: This terms are in Eq.A16 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
            deltalm = vh3 *(hCoeffs->delta55vh3 +vh3*(hCoeffs->delta55vh6 +vh3 *(hCoeffs->delta55vh9)))
                    + hCoeffs->delta55v5 * v2 * v2 * v;
            //RC: This terms are in Eq.A9 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
            rholm = 1. + v2 * (hCoeffs->rho55v2 + v * (hCoeffs->rho55v3 + v * (hCoeffs->rho55v4 + v * (hCoeffs->rho55v5 +
                    v * (hCoeffs->rho55v6 + hCoeffs->rho55v6l * eulerlogxabs +
                    v2 * (hCoeffs->rho55v8 + hCoeffs->rho55v8l * eulerlogxabs +
                    v2 * (hCoeffs->rho55v10 + hCoeffs->rho55v10l * eulerlogxabs )))))));
        //RC: This terms are in Eq.A12 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
            auxflm = v2 * v * (hCoeffs->f55v3 + v * (hCoeffs->f55v4 + v * (hCoeffs->f55v5c) ));
            // facEcc = prec_CalculateWaveformFactor_mode55(vars, eta, dm);
            break;
        case 4:
            deltalm = vh3 * (hCoeffs->delta54vh3 + hCoeffs->delta54vh4 * vh);
            rholm = 1. + v2 * (hCoeffs->rho54v2 + v * (hCoeffs->rho54v3
                            + hCoeffs->rho54v4 * v));
            break;
        case 3:
            deltalm = hCoeffs->delta53vh3 * vh3;
            rholm = 1. + v2 * (hCoeffs->rho53v2
                        + v * (hCoeffs->rho53v3 + v * (hCoeffs->rho53v4 + hCoeffs->rho53v5 * v)));
            break;
        case 2:
            deltalm = vh3 * (hCoeffs->delta52vh3 + hCoeffs->delta52vh4 * vh);
            rholm = 1. + v2 * (hCoeffs->rho52v2 + v * (hCoeffs->rho52v3
                            + hCoeffs->rho52v4 * v));
            break;
        case 1:
            deltalm = hCoeffs->delta51vh3 * vh3;
            rholm = 1. + v2 * (hCoeffs->rho51v2
                        + v * (hCoeffs->rho51v3 + v * (hCoeffs->rho51v4 + hCoeffs->rho51v5 * v)));
            break;
        default:
            return CEV_FAILURE;
            break;
        }
        break;
    case 6:
        switch (m) {
        case 6:
            deltalm = hCoeffs->delta66vh3 * vh3;
            rholm = 1. + v2 * (hCoeffs->rho66v2 + v * (hCoeffs->rho66v3
                            + hCoeffs->rho66v4 * v));
            break;
        case 5:
            deltalm = hCoeffs->delta65vh3 * vh3;
            rholm = 1. + v2 * (hCoeffs->rho65v2 + hCoeffs->rho65v3 * v);
            break;
        case 4:
            deltalm = hCoeffs->delta64vh3 * vh3;
            rholm = 1. + v2 * (hCoeffs->rho64v2 + v * (hCoeffs->rho64v3
                            + hCoeffs->rho64v4 * v));
            break;
        case 3:
            deltalm = hCoeffs->delta63vh3 * vh3;
            rholm = 1. + v2 * (hCoeffs->rho63v2 + hCoeffs->rho63v3 * v);
            break;
        case 2:
            deltalm = hCoeffs->delta62vh3 * vh3;
            rholm = 1. + v2 * (hCoeffs->rho62v2 + v * (hCoeffs->rho62v3
                            + hCoeffs->rho62v4 * v));
            break;
        case 1:
            deltalm = hCoeffs->delta61vh3 * vh3;
            rholm = 1. + v2 * (hCoeffs->rho61v2 + hCoeffs->rho61v3 * v);
            break;
        default:
            return CEV_FAILURE;
            break;
        }
        break;
    case 7:
        switch (m) {
        case 7:
            deltalm = hCoeffs->delta77vh3 * vh3;
            rholm = 1. + v2 * (hCoeffs->rho77v2 + hCoeffs->rho77v3 * v);
            break;
        case 6:
            deltalm = 0.0;
            rholm = 1. + hCoeffs->rho76v2 * v2;
            break;
        case 5:
            deltalm = hCoeffs->delta75vh3 * vh3;
            rholm = 1. + v2 * (hCoeffs->rho75v2 + hCoeffs->rho75v3 * v);
            break;
        case 4:
            deltalm = 0.0;
            rholm = 1. + hCoeffs->rho74v2 * v2;
            break;
        case 3:
            deltalm = hCoeffs->delta73vh3 * vh3;
            rholm = 1. + v2 * (hCoeffs->rho73v2 + hCoeffs->rho73v3 * v);
            break;
        case 2:
            deltalm = 0.0;
            rholm = 1. + hCoeffs->rho72v2 * v2;
            break;
        case 1:
            deltalm = hCoeffs->delta71vh3 * vh3;
            rholm = 1. + v2 * (hCoeffs->rho71v2 + hCoeffs->rho71v3 * v);
            break;
        default:
            return CEV_FAILURE;
            break;
        }
        break;
    case 8:
        switch (m) {
        case 8:
            deltalm = 0.0;
            rholm = 1. + hCoeffs->rho88v2 * v2;
            break;
        case 7:
            deltalm = 0.0;
            rholm = 1. + hCoeffs->rho87v2 * v2;
            break;
        case 6:
            deltalm = 0.0;
            rholm = 1. + hCoeffs->rho86v2 * v2;
            break;
        case 5:
            deltalm = 0.0;
            rholm = 1. + hCoeffs->rho85v2 * v2;
            break;
        case 4:
            deltalm = 0.0;
            rholm = 1. + hCoeffs->rho84v2 * v2;
            break;
        case 3:
            deltalm = 0.0;
            rholm = 1. + hCoeffs->rho83v2 * v2;
            break;
        case 2:
            deltalm = 0.0;
            rholm = 1. + hCoeffs->rho82v2 * v2;
            break;
        case 1:
            deltalm = 0.0;
            rholm = 1. + hCoeffs->rho81v2 * v2;
            break;
        default:
            return CEV_FAILURE;
            break;
        }
        break;
    default:
        return CEV_FAILURE;
        break;
    }

    //debugPK
    //     if (debugPK)
    //     XLAL_PRINT_INFO("rho_%d_%d = %.12e + %.12e \n", l, m, creal(rholm),cimag(rholm));
    // if (debugPK)
    //     XLAL_PRINT_INFO("exp(delta_%d_%d) = %.16e + i%.16e (abs=%e)\n", l, m, creal(cexp(I * deltalm)),
    //             cimag(cexp(I * deltalm)), cabs(cexp(I * deltalm)));
    /* Raise rholm to the lth power */
    rholmPwrl = 1.0;
    i = l;
    while (i--) {
        rholmPwrl *= rholm;
    }
    /*
        * In the equal-mass odd m case, there is no contribution from
        * nonspin terms,  and the only contribution comes from the auxflm
        * term that is proportional to chiA (asymmetric spins). In this
        * case, we must ignore the nonspin terms directly, since the leading
        * term defined by CalculateThisMultipolePrefix in
        * LALSimIMREOBNewtonianMultipole.c is not zero (see comments there).
        */
    if (eta == 0.25 && m % 2) {
        rholmPwrl = auxflm + facEcc;
    } else {
        rholmPwrl += auxflm + facEcc;
    }
    // if (m==1) {
    // 	printf("f21v1 = %.16f f21v3 = %.16f f21v4 = %.16f f21v5 = %.16f\n", hCoeffs->f21v1, hCoeffs->f21v3, hCoeffs->f21v4, hCoeffs->f21v5);
    // }

    // if (r > 0.0 && debugPK) {
    //     XLAL_PRINT_INFO("YP::dynamics variables in waveform: %i, %i, r = %.12e, v = %.12e\n", l, m, r, v);
    //     XLAL_PRINT_INFO("rholm^l = %.16e + %.16e, Tlm = %.16e + i %.16e, \nSlm = %.16e, hNewton = %.16e + i %.16e, delta = %.16e\n", creal(rholmPwrl), cimag(rholmPwrl), creal(Tlm), cimag(Tlm), Slm, creal(hNewton), cimag(hNewton), 0.0);
    // }
    /* Put all factors in Eq. 17 together */
    *hlm = Tlm * cexp(I * deltalm) * Slm * rholmPwrl;
    *hlm *= hNewton;
    // if (r > 8.5 && debugPK) {
    //     XLAL_PRINT_INFO("YP::FullWave: %.16e,%.16e, %.16e\n", creal(*hlm), cimag(*hlm), cabs(*hlm));
    // }
    return CEV_SUCCESS;
}


/* ---------------------------------------------------------------------------- */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/*                                    V2                                        */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/* ---------------------------------------------------------------------------- */
static COMPLEX16 prec_CalculateWaveformFactor_mode22_v2(SEOBPrecWaveformVariables *vars,
    REAL8 eta, REAL8 dm, REAL8 v, REAL8 eulerlogxabs, FacWaveformCoeffs *hCoeffs)
{
    // REAL8 x, x2, x3, x4, x5, x6, x7, x8;
    REAL8 eta2 = eta*eta;
    REAL8 nchia, nchis, lchia, lchis, echia, echis;
    REAL8 nchia2, nchis2, lchia2, lchis2, echia2, echis2;
    REAL8 chi1chi1, chi1chi2, chi2chi2;
    REAL8 Jn, Jl, Je;
    REAL8 Jn2, Jl2, Je2;
    REAL8 x, x2, x3, x4;
    REAL8 y, y2, y3, y4, y5, y6;
    REAL8 z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13, z14, z15, z16, z18, z20, z21, z22, z24, z30;
    nchia = vars->nchia;
    nchis = vars->nchis;
    lchia = vars->lchia;
    lchis = vars->lchis;
    echia = vars->echia;
    echis = vars->echis;
    nchia2 = nchia*nchia;
    nchis2 = nchis*nchis;
    lchia2 = lchia*lchia;
    lchis2 = lchis*lchis;
    echia2 = echia*echia;
    echis2 = echis*echis;
    chi1chi1 = vars->chi1chi1;
    chi1chi2 = vars->chi1chi2;
    chi2chi2 = vars->chi2chi2;
    Jn = vars->Jn;
    Jl = vars->Jl;
    Je = vars->Je;
    Jn2 = Jn*Jn;
    Jl2 = Jl*Jl;
    Je2 = Je*Je;
    x = vars->x;x2 = vars->x2;x3 = vars->x3;x4 = vars->x4;
    y = vars->y;y2 = vars->y2;y3 = vars->y3;y4 = vars->y4;y5 = vars->y5;y6 = vars->y6;
    z2 = vars->z2;z3 = vars->z3;z4 = vars->z4;z5 = vars->z5;z6 = vars->z6;
    z7 = vars->z7;z8 = vars->z8;z9 = vars->z9;z10 = vars->z10;z11 = vars->z11;
    z12 = vars->z12;z13 = vars->z13;z14 = vars->z14;z15 = vars->z15;z16 = vars->z16;
    z18 = vars->z18;z20 = vars->z20;z21 = vars->z21;z22 = vars->z22;z24 = vars->z24;
    z30 = vars->z30;
    REAL8 v2 = v*v;
    COMPLEX16 rho22v2, rho22v3, rho22v4;
    rho22v2 = (-86 + 55*eta + 21*lchia + 21*dm*lchis + I*(21)*nchia + I*(21)*dm*nchis)/(84.);
    rho22v3 = (-2*(dm*echia + echis - echis*eta))/(3.);
    rho22v4 = (-82220 + 21168*echia2 + 42336*dm*echia*echis + 21168*echis2 - 66050*eta - 84672*echia2*eta + 19583*eta2 + 5544*lchia + 
        25326*eta*lchia + 9261*lchia2 - 42336*eta*lchia2 + 5544*dm*lchis + 3150*dm*eta*lchis + 18522*dm*lchia*lchis + 
        42336*dm*eta*lchia*lchis + 9261*lchis2 + 47628*eta*lchis2 - 42336*eta2*lchis2 + I*(37296)*nchia - I*(140994)*eta*nchia - 
        I*(66150)*lchia*nchia + I*(254016)*eta*lchia*nchia - I*(66150)*dm*lchis*nchia + 1323*nchia2 + I*(37296)*dm*nchis - 
        I*(21042)*dm*eta*nchis - I*(66150)*dm*lchia*nchis - I*(66150)*lchis*nchis + I*(10584)*eta*lchis*nchis + 2646*dm*nchia*nchis - 
        42336*dm*eta*nchia*nchis + 1323*nchis2 - 47628*eta*nchis2 + 42336*eta2*nchis2)/(42336.);
    return 1. + v2 * (rho22v2 + v * (rho22v3
                                + v * (rho22v4
                    + v * (hCoeffs->rho22v5 + v * (hCoeffs->rho22v6
                                + hCoeffs->rho22v6l * eulerlogxabs + v * (hCoeffs->rho22v7
                                                        + v * (hCoeffs->rho22v8 + hCoeffs->rho22v8l * eulerlogxabs
                                                            + (hCoeffs->rho22v10 + hCoeffs->rho22v10l * eulerlogxabs) * v2)))))));
}

static COMPLEX16 prec_CalculateWaveformFactor_mode21_v2(SEOBPrecWaveformVariables *vars,
    REAL8 eta, REAL8 dm, REAL8 v, REAL8 eulerlogxabs, FacWaveformCoeffs *hCoeffs)
{
    // REAL8 x, x2, x3, x4, x5, x6, x7, x8;
    REAL8 eta2 = eta*eta;
    REAL8 nchia, nchis, lchia, lchis, echia, echis;
    REAL8 nchia2, nchis2, lchia2, lchis2, echia2, echis2;
    REAL8 chi1chi1, chi1chi2, chi2chi2;
    REAL8 Jn, Jl, Je;
    REAL8 Jn2, Jl2, Je2;
    REAL8 x, x2, x3, x4;
    REAL8 y, y2, y3, y4, y5, y6;
    REAL8 z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13, z14, z15, z16, z18, z20, z21, z22, z24, z30;
    nchia = vars->nchia;
    nchis = vars->nchis;
    lchia = vars->lchia;
    lchis = vars->lchis;
    echia = vars->echia;
    echis = vars->echis;
    nchia2 = nchia*nchia;
    nchis2 = nchis*nchis;
    lchia2 = lchia*lchia;
    lchis2 = lchis*lchis;
    echia2 = echia*echia;
    echis2 = echis*echis;
    chi1chi1 = vars->chi1chi1;
    chi1chi2 = vars->chi1chi2;
    chi2chi2 = vars->chi2chi2;
    Jn = vars->Jn;
    Jl = vars->Jl;
    Je = vars->Je;
    Jn2 = Jn*Jn;
    Jl2 = Jl*Jl;
    Je2 = Je*Je;
    x = vars->x;x2 = vars->x2;x3 = vars->x3;x4 = vars->x4;
    y = vars->y;y2 = vars->y2;y3 = vars->y3;y4 = vars->y4;y5 = vars->y5;y6 = vars->y6;
    z2 = vars->z2;z3 = vars->z3;z4 = vars->z4;z5 = vars->z5;z6 = vars->z6;
    z7 = vars->z7;z8 = vars->z8;z9 = vars->z9;z10 = vars->z10;z11 = vars->z11;
    z12 = vars->z12;z13 = vars->z13;z14 = vars->z14;z15 = vars->z15;z16 = vars->z16;
    z18 = vars->z18;z20 = vars->z20;z21 = vars->z21;z22 = vars->z22;z24 = vars->z24;
    z30 = vars->z30;
    REAL8 v2 = v*v;
    COMPLEX16 f21v1, f21v2, f21v3, auxflm;
    if (dm)
    {
        f21v1 = -3.*echis/2. - 3.*echia/2./dm;
        f21v2 = -6.*I*nchia - 6.*I*nchis/dm + 2.*I*nchis*eta/dm + 6.*lchia + 6.*lchis/dm - 2.*eta*lchis/dm;
        f21v3 = 61.*echis/12. + 15.*I*echis*nchia/2. + 15.*I*echia*nchis/2. + 61.*echia/12./dm + 
            15.*I*echia*nchia/2./dm + 15.*I*echis*nchis/2./dm + 79.*echis*eta/84. + 
            3.*I*echis*nchia*eta + 3.*I*echia*nchis*eta + 383.*echia*eta/84./dm - 
            30.*I*echia*nchia*eta/dm + 6.*I*echis*nchis*eta/dm - 6.*I*echis*nchis*eta2/dm - 
            3.*echis*lchia/2. - 3.*echia*lchia/2./dm + 
            3.*echis*eta*lchia + 6.*echia*eta*lchia/dm - 3.*echia*lchis/2. - 
            3.*echis*lchis/2./dm + 3.*echia*eta*lchis + 6.*echis*eta*lchis/dm - 6.*echis*eta2*lchis/dm;
    } else {
        f21v1 = - 3.*echia/2.;
        f21v2 =  - 6.*I*nchis + 2.*I*nchis*eta + 6.*lchis - 2.*eta*lchis;
        f21v3 = 61.*echia/12. + 15.*I*echia*nchia/2. + 15.*I*echis*nchis/2. + 383.*echia*eta/84. - 
            30.*I*echia*nchia*eta + 6.*I*echis*nchis*eta - 6.*I*echis*nchis*eta2 - 3.*echia*lchia/2. + 
            6.*echia*eta*lchia  - 3.*echis*lchis/2. + 6.*echis*eta*lchis - 6.*echis*eta2*lchis;
    }
    auxflm = v * (f21v1 + 
        v * (f21v2 + 
        v * (f21v3 + 
        v * hCoeffs->f21v4 + 
        v2 * (hCoeffs->f21v5 + v * hCoeffs->f21v6 + v2 * hCoeffs->f21v7c))));

    return auxflm;
}

static COMPLEX16 prec_CalculateWaveformFactor_mode33_v2(SEOBPrecWaveformVariables *vars,
    REAL8 eta, REAL8 dm, REAL8 v, REAL8 vh3, REAL8 eulerlogxabs, FacWaveformCoeffs *hCoeffs)
{
    // REAL8 x, x2, x3, x4, x5, x6, x7, x8;
    REAL8 eta2 = eta*eta;
    REAL8 nchia, nchis, lchia, lchis, echia, echis;
    REAL8 nchia2, nchis2, lchia2, lchis2, echia2, echis2;
    REAL8 chi1chi1, chi1chi2, chi2chi2;
    REAL8 Jn, Jl, Je;
    REAL8 Jn2, Jl2, Je2;
    REAL8 x, x2, x3, x4;
    REAL8 y, y2, y3, y4, y5, y6;
    REAL8 z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13, z14, z15, z16, z18, z20, z21, z22, z24, z30;
    nchia = vars->nchia;
    nchis = vars->nchis;
    lchia = vars->lchia;
    lchis = vars->lchis;
    echia = vars->echia;
    echis = vars->echis;
    nchia2 = nchia*nchia;
    nchis2 = nchis*nchis;
    lchia2 = lchia*lchia;
    lchis2 = lchis*lchis;
    echia2 = echia*echia;
    echis2 = echis*echis;
    chi1chi1 = vars->chi1chi1;
    chi1chi2 = vars->chi1chi2;
    chi2chi2 = vars->chi2chi2;
    Jn = vars->Jn;
    Jl = vars->Jl;
    Je = vars->Je;
    Jn2 = Jn*Jn;
    Jl2 = Jl*Jl;
    Je2 = Je*Je;
    x = vars->x;x2 = vars->x2;x3 = vars->x3;x4 = vars->x4;
    y = vars->y;y2 = vars->y2;y3 = vars->y3;y4 = vars->y4;y5 = vars->y5;y6 = vars->y6;
    z2 = vars->z2;z3 = vars->z3;z4 = vars->z4;z5 = vars->z5;z6 = vars->z6;
    z7 = vars->z7;z8 = vars->z8;z9 = vars->z9;z10 = vars->z10;z11 = vars->z11;
    z12 = vars->z12;z13 = vars->z13;z14 = vars->z14;z15 = vars->z15;z16 = vars->z16;
    z18 = vars->z18;z20 = vars->z20;z21 = vars->z21;z22 = vars->z22;z24 = vars->z24;
    z30 = vars->z30;
    REAL8 v2 = v*v;
    COMPLEX16 f33v2, f33v3, auxflm;
    if (dm)
    {
        f33v2 = -16.*(I*nchis + lchis)*eta/9./dm;
        f33v3 = echis*(-2. + 5.*eta/2.) + echia*(-2. + 19.*eta/2.) / dm;
    } else {
        f33v2 = -16.*(I*nchis + lchis)*eta/9.;
        f33v3 = echia*(-2. + 19.*eta/2.);
    }
    auxflm = v2 * (f33v2 + v * (f33v3 + v * (hCoeffs->f33v4 + v * (hCoeffs->f33v5  + v * hCoeffs->f33v6)))) + _Complex_I * vh3 * vh3 * hCoeffs->f33vh6;

    return auxflm;
}

static COMPLEX16 prec_CalculateWaveformFactor_mode44_v2(SEOBPrecWaveformVariables *vars,
    REAL8 eta, REAL8 dm, REAL8 v, REAL8 eulerlogxabs, FacWaveformCoeffs *hCoeffs)
{
    REAL8 eta2 = eta*eta;
    REAL8 v2, m1Plus3eta;
    v2 = v*v;
    m1Plus3eta = -1. + 3.*eta;
    REAL8 nchia, nchis, lchia, lchis, echia, echis;
    nchia = vars->nchia;
    nchis = vars->nchis;
    lchia = vars->lchia;
    lchis = vars->lchis;
    echia = vars->echia;
    echis = vars->echis;
    COMPLEX16 rho44v2;
    rho44v2 = hCoeffs->rho44v2 + 81.*eta*(-I*nchia + I*nchis*dm - lchia + dm*lchis) / (256.*m1Plus3eta);
    COMPLEX16 rholm;
    rholm = 1. + v2 * (rho44v2
        + v * (hCoeffs->rho44v3 + v * (hCoeffs->rho44v4 + v * (hCoeffs->rho44v5 + v * (hCoeffs->rho44v6 + hCoeffs->rho44v6l *
        eulerlogxabs + v2 *( hCoeffs->rho44v8 + hCoeffs->rho44v8l * eulerlogxabs +v2 * (hCoeffs->rho44v10 + hCoeffs->rho44v10l * eulerlogxabs) ) )))));

    return rholm;
}

/* ---------------------------------------------------------------------------- */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/*                                    V3                                        */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/* ---------------------------------------------------------------------------- */
static COMPLEX16 prec_CalculateWaveformFactor_mode22_v3(SEOBPrecWaveformVariables *vars,
    REAL8 eta, REAL8 dm, REAL8 v, REAL8 eulerlogxabs, FacWaveformCoeffs *hCoeffs, COMPLEX16 *fEcc)
{
    // REAL8 x, x2, x3, x4, x5, x6, x7, x8;
    REAL8 eta2 = eta*eta;
    REAL8 nchia, nchis, lchia, lchis, echia, echis;
    REAL8 nchia2, nchis2, lchia2, lchis2, echia2, echis2;
    nchia = vars->nchia;
    nchis = vars->nchis;
    lchia = vars->lchia;
    lchis = vars->lchis;
    echia = vars->echia;
    echis = vars->echis;
    nchia2 = nchia*nchia;
    nchis2 = nchis*nchis;
    lchia2 = lchia*lchia;
    lchis2 = lchis*lchis;
    echia2 = echia*echia;
    echis2 = echis*echis;
    REAL8 x, x2, x3, x4;
    REAL8 y, y2, y3, y4, y5, y6;
    REAL8 z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13, z14, z15, z16, z18, z20, z21, z22, z24, z30;
    x = vars->x;x2 = vars->x2;x3 = vars->x3;x4 = vars->x4;
    y = vars->y;y2 = vars->y2;y3 = vars->y3;y4 = vars->y4;y5 = vars->y5;y6 = vars->y6;
    z2 = vars->z2;z3 = vars->z3;z4 = vars->z4;z5 = vars->z5;z6 = vars->z6;
    z7 = vars->z7;z8 = vars->z8;z9 = vars->z9;z10 = vars->z10;z11 = vars->z11;
    z12 = vars->z12;z13 = vars->z13;z14 = vars->z14;z15 = vars->z15;z16 = vars->z16;
    z18 = vars->z18;z20 = vars->z20;z21 = vars->z21;z22 = vars->z22;z24 = vars->z24;
    z30 = vars->z30;
    COMPLEX16 fEpn0, fEpn2, fEpn3, fEpn4;
    fEpn0 = -(-1. + 2.*z4 + cpow(y - I*z6, 2.)) / (2.*z4);
    fEpn2 = x2*(-14*(1 + eta) + 21*(-7 + 3*eta + 2*lchia + (2*I)*nchia + 2*dm*(lchis + I*nchis))*z12 - 112*z18 + 2*(86 - 55*eta - 21*lchia - (21*I)*nchia - (21*I)*dm*((-I)*lchis + nchis))*z20 + (53 + 30*eta)*z24 + (-8 + 31*eta)*z20*z16 + y4*(-7*(-3 + eta) + (29 - 10*eta)*z12 - 56*z6) + 56*z6 + I*y3*(-7*(-3 + eta) + (-37 + 41*eta)*z12 - 56*z6)*z6 + I*y*(14*(1 + eta) + (15 - 101*eta)*z12 + (-37 + 41*eta)*z24 - 56*z6)*z6 + y2*(-7 + 21*eta + z12*(74 + 9*eta + 7*z6*(-16 + 3*(1 + eta)*z6)))) / (84.*z16);
    fEpn3 = x3*(-2*dm*echia*(-1 + y2 + 5*z12 - 4*z16 + (4*I)*y*z6) + echis*(2 - eta + (-2 + eta)*y2 + z12*(-10 + 9*eta - 8*(-1 + eta)*z4) + I*(-8 + 5*eta)*y*z6)) / (6.*z10);
    fEpn4 = (x4*(6*y6*(-28*(-3 + eta)*(-3 + eta) + (-1399 + 13*eta + 43*eta2)*z12 + 16*(50 + 11*eta)*z18 + (-277 + 130*eta + 80*eta2)*z24 - 448*(-3 + eta)*z6) + 6*y2*(-224*(-1 + eta2) - 2*(-853 + 614*eta - 168*dm*echia*echis*(1 + eta) + 84*echia2*(-1 + 4*eta) - 269*eta2 + 84*echis2*(-1 - 2*eta + 2*eta2) + 
        126*lchia - 42*eta*lchia - 84*lchia2 + 336*eta*lchia2 + 126*dm*lchis - 42*dm*eta*lchis - 168*dm*lchia*lchis - 168*dm*eta*lchia*lchis - 84*lchis2 - 168*eta*lchis2 + 168*eta2*lchis2 + (126*I)*nchia - (42*I)*eta*nchia + 252*nchia2 - 1008*eta*nchia2 + (126*I)*dm*nchis - (42*I)*dm*eta*nchis + 504*dm*nchia*nchis + 252*nchis2)*z12 + 
        16*(-31 + 51*eta + 42*lchia + 42*dm*lchis + (42*I)*nchia + (42*I)*dm*nchis)*z18 + 4*(-3 + eta)*(-86 + 55*eta + 21*lchia + 21*dm*lchis + (21*I)*nchia + (21*I)*dm*nchis)*z20 + 2*(1595 - 420*echis2*(1 - 2*eta)*(1 - 2*eta) - 1620*eta + 840*dm*echia*echis*(-1 + 2*eta) + 420*echia2*(-1 + 4*eta) - 205*eta2 - 714*lchia - 402*eta*lchia - 420*lchia2 + 
        1680*eta*lchia2 - 714*dm*lchis - 186*dm*eta*lchis - 840*dm*lchia*lchis + 1680*dm*eta*lchia*lchis - 420*lchis2 + 1680*eta*lchis2 - 1680*eta2*lchis2 - (714*I)*nchia - (300*I)*eta*nchia + (252*I)*lchia*nchia - (1008*I)*eta*lchia*nchia + (252*I)*dm*lchis*nchia - (504*I)*dm*eta*lchis*nchia - (714*I)*dm*nchis - (120*I)*dm*eta*nchis + 
        (252*I)*dm*lchia*nchis - (504*I)*dm*eta*lchia*nchis + (252*I)*lchis*nchis - (1008*I)*eta*lchis*nchis + (1008*I)*eta2*lchis*nchis)*z24 + 32*(-86 + 55*eta + 21*lchia + 21*dm*lchis + (21*I)*nchia + (21*I)*dm*nchis)*z20*z6 + 32*(-116 + 33*eta)*z30 + 4*(1 + eta)*(-86 + 55*eta + 21*lchia + 21*dm*lchis + (21*I)*nchia + 
        (21*I)*dm*nchis)*z20*z12 + (2073 + 1075*eta - 1181*eta2)*z20*z16 - 80*(-8 + 31*eta)*z20*z22 + (-71 + 20*eta + 40*eta2)*z24*z24 + 448*(-3 + eta)*z6) - (3*I)*y5*z6*(35*(-3 + eta)*(-3 + eta) - 2*(-664 + 55*eta + 22*eta2)*z12 - 16*(79 + eta)*z18 + (-571 + 748*eta + 389*eta2)*z24 + 560*(-3 + eta)*z6) - I*y*z6*(420*(1 + eta)*(1 + eta) - 
        12*(-419 + eta2*(25 - 504*echis2 - 504*lchis2) - 252*nchia2 + 4*eta*(-137 + 126*echis2 + 126*dm*(echia*echis + lchia*lchis) + 126*lchis2 + 252*nchia2) - 504*dm*nchia*nchis - 252*nchis2)*z12 - 48*(-15 + 227*eta)*z18 + (11717 - 6048*echis2*(1 - 2*eta)*(1 - 2*eta) - 3190*eta + 12096*dm*echia*echis*(-1 + 2*eta) + 
        6048*echia2*(-1 + 4*eta) + 6721*eta2 - 15120*lchia - 504*eta*lchia - 3024*lchia2 + 12096*eta*lchia2 - 15120*dm*lchis - 1512*dm*eta*lchis - 6048*dm*lchia*lchis + 12096*dm*eta*lchia*lchis - 3024*lchis2 + 12096*eta*lchis2 - 12096*eta2*lchis2 - (9072*I)*nchia - (14040*I)*eta*nchia + (6048*I)*lchia*nchia - 
        (24192*I)*eta*lchia*nchia + (6048*I)*dm*lchis*nchia - (12096*I)*dm*eta*lchis*nchia - 3024*nchia2 + 12096*eta*nchia2 - (9072*I)*dm*nchis + (3960*I)*dm*eta*nchis + (6048*I)*dm*lchia*nchis - (12096*I)*dm*eta*lchia*nchis + (6048*I)*lchis*nchis - (24192*I)*eta*lchis*nchis + (24192*I)*eta2*lchis*nchis - 
        6048*dm*nchia*nchis + 12096*dm*eta*nchia*nchis - 3024*nchis2 + 12096*eta*nchis2 - 12096*eta2*nchis2)*z24 + 336*(-37 + 41*eta)*z30 + 6*(2880 - 4811*eta + 640*eta2)*z20*z16 - 3*(127 - 700*eta + 103*eta2)*z24*z24 - 3360*(1 + eta)*z6) + 6*y4*(28*(-3 - 14*eta + 5*eta2) + 
        3*(-567 - 56*echis2*(1 - 2*eta)*(1 - 2*eta) - 183*eta + 112*dm*echia*echis*(-1 + 2*eta) + 56*echia2*(-1 + 4*eta) - 81*eta2 - 56*lchia2 + 224*eta*lchia2 - 112*dm*lchia*lchis + 224*dm*eta*lchia*lchis - 56*lchis2 + 224*eta*lchis2 - 224*eta2*lchis2)*z12 + 32*(20 + 31*eta)*z18 - 
        7*(-163 + 12*eta + 42*eta2)*z24 + (-340 + 193*eta + 17*eta2)*z20*z16 + 896*(1 + eta)*z6) - (6*I)*y3*z6*(-70*(-3 - 2*eta + eta2) + (1177 + 168*echis2*(1 - 2*eta)*(1 - 2*eta) + 898*eta - 336*dm*echia*echis*(-1 + 2*eta) - 168*echia2*(-1 + 4*eta) + 69*eta2 + 168*lchia2 - 672*eta*lchia2 + 
        336*dm*lchia*lchis - 672*dm*eta*lchia*lchis + 168*lchis2 - 672*eta*lchis2 + 672*eta2*lchis2)*z12 + (496 - 1824*eta)*z18 - 2*(-509 + 1459*eta + 281*eta2)*z24 + 56*(-37 + 41*eta)*z30 + (-349 + 724*eta + 143*eta2)*z20*z16 - 280*(5 + eta)*z6) + 6*(112*(1 + eta)*(1 + eta) - 
        28*(-4*eta2*(-4 + 9*echis2 + 9*lchis2) + eta*(-49 + 36*echis2 + 6*lchia + 36*lchis2 + (6*I)*nchia + 72*nchia2 + 6*dm*(6*echia*echis + lchis + 6*lchia*lchis + I*nchis)) + 6*(-13 + lchia + I*nchia - 3*nchia2 + dm*(lchis + I*nchis - 6*nchia*nchis) - 3*nchis2))*z12 + 
        336*(-7 + 2*lchia + (2*I)*nchia + 2*dm*(lchis + I*nchis))*z18 - 8*(1 + eta)*(-86 + 55*eta + 21*lchia + 21*dm*lchis + (21*I)*nchia + (21*I)*dm*nchis)*z20 - 2*(-853 + 620*eta - 504*dm*echia*echis*(1 + 3*eta) + 252*echia2*(-1 + 4*eta) - 733*eta2 + 252*echis2*(-1 - 6*eta + 6*eta2) + 
        462*lchia + 498*eta*lchia - 252*lchia2 + 1008*eta*lchia2 + 462*dm*lchis - 246*dm*eta*lchis - 504*dm*lchia*lchis - 1512*dm*eta*lchia*lchis - 252*lchis2 - 1512*eta*lchis2 + 1512*eta2*lchis2 + (798*I)*nchia - (642*I)*eta*nchia + (1260*I)*lchia*nchia - (5040*I)*eta*lchia*nchia + 
        (1260*I)*dm*lchis*nchia + (504*I)*dm*eta*lchis*nchia + 504*nchia2 - 2016*eta*nchia2 + (798*I)*dm*nchis - (282*I)*dm*eta*nchis + (1260*I)*dm*lchia*nchis + (504*I)*dm*eta*lchia*nchis + (1260*I)*lchis*nchis + (1008*I)*eta*lchis*nchis - (1008*I)*eta2*lchis*nchis + 1008*dm*nchia*nchis + 504*nchis2)*z24 + 
        32*(-86 + 55*eta + 21*lchia + 21*dm*lchis + (21*I)*nchia + (21*I)*dm*nchis)*z20*z6 + 32*(-53 + 33*eta)*z30 - 4*(7 + eta)*(-86 + 55*eta + 21*lchia + 21*dm*lchis + (21*I)*nchia + (21*I)*dm*nchis)*z20*z12 + (632 - 101*eta - 1008*dm*echia*echis*(1 + 2*eta) + 504*echia2*(-1 + 4*eta) - 409*eta2 + 
        504*echis2*(-1 - 4*eta + 4*eta2) + 420*lchia + 1332*eta*lchia - 504*lchia2 + 2016*eta*lchia2 + 420*dm*lchis - 156*dm*eta*lchis - 1008*dm*lchia*lchis - 2016*dm*eta*lchia*lchis - 504*lchis2 - 2016*eta*lchis2 + 2016*eta2*lchis2 + (1092*I)*nchia - (948*I)*eta*nchia + (2520*I)*lchia*nchia - 
        (10080*I)*eta*lchia*nchia + (2520*I)*dm*lchis*nchia + (1008*I)*dm*eta*lchis*nchia + 504*nchia2 - 2016*eta*nchia2 + (1092*I)*dm*nchis - (228*I)*dm*eta*nchis + (2520*I)*dm*lchia*nchis + (1008*I)*dm*eta*lchia*nchis + (2520*I)*lchis*nchis + (2016*I)*eta*lchis*nchis - 
        (2016*I)*eta2*lchis*nchis + 1008*dm*nchia*nchis + 504*nchis2)*z20*z16 - 80*(-8 + 31*eta)*z20*z22 + 4*(1 + eta)*(-86 + 55*eta + 21*lchia + 21*dm*lchis + (21*I)*nchia + (21*I)*dm*nchis)*z22*z22 - 2*(161 - 710*eta + 192*eta2)*z24*z24 - 896*(1 + eta)*z6 + (-8 - 43*eta + 103*eta2)*z30*z30))) / (6048*z14*z14);
    if (PREC_FLAG > 3 || PREC_FLAG==1)
        *fEcc = 0.0;
    else
        *fEcc = fEpn0 + fEpn2 + fEpn3 + fEpn4;
    REAL8 v2 = v*v;
    COMPLEX16 rho22v2, rho22v3, rho22v4;
    rho22v2 = (-86 + 55*eta + 21.*I*nchia + 21.*dm*(I*nchis + lchis))/(84.);
    rho22v3 = (-2.*(dm*echia + echis - echis*eta))/(3.);
    rho22v4 = dm*echia*echis + echia2*(1/2 - 2*eta) + (19583*eta2)/42336 - 
        (121*eta*lchia)/336 + (7*lchia2)/32 - eta*lchia2 + (13*dm*lchis)/28 + 
        (11*dm*eta*lchis)/336 + (7*dm*lchia*lchis)/16 + dm*eta*lchia*lchis + 
        (7*lchis2)/32 + (9*eta*lchis2)/8 - eta2*lchis2 + (I/112)*(eta*(-317 + 672*lchia) - 
        175*(lchia + dm*lchis))*nchia - ((I*(-296*dm + 167*dm*eta + 525*dm*lchia + 525*lchis - 84*eta*lchis) + 
        21*dm*(-1 + 16*eta)*nchia)*nchis)/336 + ((-9*eta)/8 + eta2)*nchis2 + 
        (-82220 + 21168*echis2 - 66050*eta + 19656*lchia + 63*nchia*(592*I + 21*nchia) + 1323*nchis2)/42336;
    if (PREC_FLAG == 4 || PREC_FLAG == 3)
        return 1. + v2 * (rho22v2 + v * (rho22v3
                                    + v * (rho22v4
                        + v * (hCoeffs->rho22v5 + v * (hCoeffs->rho22v6
                                    + hCoeffs->rho22v6l * eulerlogxabs + v * (hCoeffs->rho22v7
                                                            + v * (hCoeffs->rho22v8 + hCoeffs->rho22v8l * eulerlogxabs
                                                                + (hCoeffs->rho22v10 + hCoeffs->rho22v10l * eulerlogxabs) * v2)))))));
    else
        return 1. + v2 * (hCoeffs->rho22v2 + v * (hCoeffs->rho22v3
                                    + v * (hCoeffs->rho22v4
                        + v * (hCoeffs->rho22v5 + v * (hCoeffs->rho22v6
                                    + hCoeffs->rho22v6l * eulerlogxabs + v * (hCoeffs->rho22v7
                                                            + v * (hCoeffs->rho22v8 + hCoeffs->rho22v8l * eulerlogxabs
                                                                + (hCoeffs->rho22v10 + hCoeffs->rho22v10l * eulerlogxabs) * v2)))))));

}

static COMPLEX16 prec_CalculateWaveformFactor_mode21_v3(SEOBPrecWaveformVariables *vars,
    REAL8 eta, REAL8 dm, REAL8 v, REAL8 eulerlogxabs, FacWaveformCoeffs *hCoeffs, COMPLEX16 *fEcc)
{
    // REAL8 x, x2, x3, x4, x5, x6, x7, x8;
    REAL8 eta2 = eta*eta;
    REAL8 nchia, nchis, lchia, lchis, echia, echis;
    REAL8 nchia2, nchis2, lchia2, lchis2, echia2, echis2;
    nchia = vars->nchia;
    nchis = vars->nchis;
    lchia = vars->lchia;
    lchis = vars->lchis;
    echia = vars->echia;
    echis = vars->echis;
    nchia2 = nchia*nchia;
    nchis2 = nchis*nchis;
    lchia2 = lchia*lchia;
    lchis2 = lchis*lchis;
    echia2 = echia*echia;
    echis2 = echis*echis;
    REAL8 x, x2, x3, x4;
    REAL8 y, y2, y3, y4, y5, y6;
    REAL8 z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13, z14, z15, z16, z18, z20, z21, z22, z24, z30;
    x = vars->x;x2 = vars->x2;x3 = vars->x3;x4 = vars->x4;
    y = vars->y;y2 = vars->y2;y3 = vars->y3;y4 = vars->y4;y5 = vars->y5;y6 = vars->y6;
    z2 = vars->z2;z3 = vars->z3;z4 = vars->z4;z5 = vars->z5;z6 = vars->z6;
    z7 = vars->z7;z8 = vars->z8;z9 = vars->z9;z10 = vars->z10;z11 = vars->z11;
    z12 = vars->z12;z13 = vars->z13;z14 = vars->z14;z15 = vars->z15;z16 = vars->z16;
    z18 = vars->z18;z20 = vars->z20;z21 = vars->z21;z22 = vars->z22;z24 = vars->z24;
    z30 = vars->z30;
    COMPLEX16 fEpn1, fEpn2, fEpn3, fEpn4;
    if (dm)
    {
        fEpn1 = -1. + 1./z8;
        fEpn2 = 3. * x * (-1. + z16) * (echia + echis*dm) / (2.*z14*dm);
        fEpn3 = (-(x2*((84*I)*(-3 + eta)*nchis*(-1 + z12)*z12 + dm*y2*(-14*(-3 + eta) + (106 + 19*eta)*z12 - 112*z6) + 
            (3*I)*y*(42*dm*lchia + 14*(3 - 2*eta)*lchis + dm*(-83 + 12*eta)*z12)*z6 + 2*(14*dm*(1 + eta) + 
            (-21*eta*lchis + dm*(97 - 26*eta + (126*I)*nchia))*z12 + (21*eta*lchis + dm*(-55 + 12*eta - (126*I)*nchia))*z24 - 56*dm*z6)))) / (42*dm*z20);
        fEpn4 = (-(x3*(echia*(-294*(1 + eta) + (-2261 + 1889*eta + (1260*I)*(-1 + 4*eta)*nchia - (252*I)*dm*(5 + 2*eta)*nchis)*z12 - 
            42*(1 + eta)*z16 - 12*((21*I)*(-7 + lchia + dm*lchis) - (2*I)*eta*(-52 + 42*lchia + 21*dm*lchis) + 21*(-1 + 4*eta)*nchia + 
            21*dm*(-1 + 2*eta)*nchis)*y*z18 + 168*z22 + 3*(175 - 84*lchia - 84*dm*lchis + eta*(-605 + 336*lchia + 168*dm*lchis))*z24 - 
            21*(7 + eta)*z20*z8 + 2*(427 + 131*eta + 126*lchia - 504*eta*lchia + 126*dm*lchis - 252*dm*eta*lchis - (630*I)*(-1 + 4*eta)*nchia + 
            (126*I)*dm*(5 + 2*eta)*nchis)*z20*z12 + 21*(1 + eta)*z20*z20 + 1176*z6 + 3*y2*(49*(-3 + eta) - I*((-I)*(203 + 141*eta) + 84*(-1 + 4*eta)*nchia + 
            84*dm*(-1 + 2*eta)*nchis)*z12 + 7*(-3 + eta)*z16 + 56*z22 + 7*(1 + eta)*z20*z8 + 392*z6)) + echis*(dm*(eta*(-294 + (673 - (504*I)*nchia)*z12 - 
            42*z16 + (72*I)*(1 + 7*lchia + (7*I)*nchia)*y*z18 + 9*(-55 + 56*lchia)*z24 - 21*z20*z8 + 3*y2*(49 + (-69 - (168*I)*nchia)*z12 + 7*z16 + 7*z20*z8) + 
            (158 - 504*lchia + (504*I)*nchia)*z20*z12 + 21*z20*z20) + 7*(-42 + (-323 - (180*I)*nchia)*z12 - 6*z16 + 36*(7*I - I*lchia + nchia)*y*z18 + 24*z22 + 
            (75 - 36*lchia)*z24 - 21*z20*z8 + 2*(61 + 18*lchia + (90*I)*nchia)*z20*z12 + 3*z20*z20 + 168*z6 + 3*y2*(-21 + (-29 + (12*I)*nchia)*z12 - 3*z16 + 8*z22 + z20*z8 + 56*z6))) - 
            (252*I)*z12*(-(nchis*(-5 + y2 + 5*z20 - I*y*z6 - 4*eta*(1 + y2 - z20 - I*y*z6) + 4*eta2*(1 + y2 - z20 - I*y*z6))) + (1 - 2*eta)*(1 - 2*eta)*lchis*z6*(y + I*z6*(-1 + z8))))))) / (168*dm*z20*z6);
    } else {
        fEpn1 = -1. + 1./z8;
        fEpn2 = 3.*echia*x*(-1. + z16) / (2.*z14);
        fEpn3 = (I*x2*(6*nchis*(-1. + z12)*z12 - 2.*eta*nchis*(-1. + z12)*z12 - 3*lchis*y*z6 + eta*lchis*(-I*z12 + I*z24 + 2.*y*z6)))/z20;
        fEpn4 = (x3*((I*252)*echis*z12*(nchis*(5 - y2 - 5*z20 + (I)*y*z6 + 4*eta*(1 + y2 - z20 - (I)*y*z6) - 4*eta2*(1 + y2 - z20 - (I*1)*y*z6)) + (1. - 2.*eta)*(1.-2.*eta)*lchis*z6*(y + (I)*z6*(-1. + z8))) + 
            echia*(7*(z12*(323 - (I*180)*nchia*(-1 + z20) - 2*(61 + 18*lchia)*z20) + 3*y2*(21 + (29 - (I*12)*nchia)*z12 + 3*z16 - 8*z22 - 56*z6 - z20*z8) + 3*(14 + 2*z16 + (I*12)*(-7 + lchia + (I)*nchia)*y*z18 - z20*z20 - 
            8*z22 - 25*z24 + 12*lchia*z24 - 56*z6 + 7*z20*z8)) + eta*(z12*(-1889 + (I*5040)*nchia*(-1 + z20) + 2*(-131 + 504*lchia)*z20) + 3*(98 + 14*z16 + 16*((I*26) - (I*21)*lchia + 21*nchia)*y*z18 - 7*(z20*z20) + 605*z24 - 
            336*lchia*z24 + 7*z20*z8) + 3*y2*(3*(47 + (I*112)*nchia)*z12 - 7*(7 + z16 + z20*z8))))))/(168.*z20*z6);
    }
    if (PREC_FLAG > 3 || PREC_FLAG==1)
        *fEcc = 0.0;
    else
        *fEcc = fEpn1 + fEpn2 + fEpn3 + fEpn4;
    REAL8 v2 = v*v;
    COMPLEX16 f21v1, f21v2, f21v3, auxflm;
    if (dm)
    {
        f21v1 = -3.*echis/2. - 3.*echia/2./dm;
        f21v2 = -6.*I*nchia + (-6.*I*nchis + 2.*I*nchis*eta + eta*lchis)/dm;
        f21v3 = (61*echis)/12. + (79*echis*eta)/84. + (3*echis*lchia)/2. - 3.*echis*eta*lchia + (3*echia*lchis)/2. - 3.*echia*eta*lchis + ((15.*I)/2)*echis*nchia + (3.*I)*echis*eta*nchia + ((15.*I)/2.)*echia*nchis + (3.*I)*echia*eta*nchis + 
            ((61*echia)/12. + (131*echia*eta)/84. + (3*echia*lchia)/2. - 6.*echia*eta*lchia + (3*echis*lchis)/2. - 6.*echis*eta*lchis + 6.*echis*eta2*lchis + ((15.*I)/2.)*echia*nchia - (30.*I)*echia*eta*nchia + ((15.*I)/2.)*echis*nchis + (6.*I)*echis*eta*nchis - (6.*I)*echis*eta2*nchis)/dm;
    } else {
        f21v1 = -3.*echia/2.;
        f21v2 = (-6.*I*nchis + 2.*I*nchis*eta + eta*lchis);
        f21v3 = (61*echia)/12. + (131*echia*eta)/84. + (3*echia*lchia)/2. - 6.*echia*eta*lchia + (3*echis*lchis)/2. - 6.*echis*eta*lchis + 6.*echis*eta2*lchis + ((15.*I)/2.)*echia*nchia - (30.*I)*echia*eta*nchia + ((15.*I)/2.)*echis*nchis + (6.*I)*echis*eta*nchis - (6.*I)*echis*eta2*nchis;
    }
    if (PREC_FLAG == 4 || PREC_FLAG == 3)
        auxflm = v * (f21v1 + 
            v * (f21v2 + 
            v * (f21v3 + 
            v * hCoeffs->f21v4 + 
            v2 * (hCoeffs->f21v5 + v * hCoeffs->f21v6 + v2 * hCoeffs->f21v7c))));
    else
        auxflm = v * (hCoeffs->f21v1 + 
            v * (0.0 + 
            v * (hCoeffs->f21v3 + 
            v * hCoeffs->f21v4 + 
            v2 * (hCoeffs->f21v5 + v * hCoeffs->f21v6 + v2 * hCoeffs->f21v7c))));

    return auxflm;
}

static COMPLEX16 prec_CalculateWaveformFactor_mode33_v3(SEOBPrecWaveformVariables *vars,
    REAL8 eta, REAL8 dm, REAL8 v, REAL8 vh3, REAL8 eulerlogxabs, FacWaveformCoeffs *hCoeffs, COMPLEX16 *fEcc)
{
    // REAL8 x, x2, x3, x4, x5, x6, x7, x8;
    REAL8 eta2 = eta*eta;
    REAL8 nchia, nchis, lchia, lchis, echia, echis;
    REAL8 nchia2, nchis2, lchia2, lchis2, echia2, echis2;
    nchia = vars->nchia;
    nchis = vars->nchis;
    lchia = vars->lchia;
    lchis = vars->lchis;
    echia = vars->echia;
    echis = vars->echis;
    nchia2 = nchia*nchia;
    nchis2 = nchis*nchis;
    lchia2 = lchia*lchia;
    lchis2 = lchis*lchis;
    echia2 = echia*echia;
    echis2 = echis*echis;
    REAL8 x, x2, x3, x4;
    REAL8 y, y2, y3, y4, y5, y6;
    REAL8 z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13, z14, z15, z16, z18, z20, z21, z22, z24, z30;
    x = vars->x;x2 = vars->x2;x3 = vars->x3;x4 = vars->x4;
    y = vars->y;y2 = vars->y2;y3 = vars->y3;y4 = vars->y4;y5 = vars->y5;y6 = vars->y6;
    z2 = vars->z2;z3 = vars->z3;z4 = vars->z4;z5 = vars->z5;z6 = vars->z6;
    z7 = vars->z7;z8 = vars->z8;z9 = vars->z9;z10 = vars->z10;z11 = vars->z11;
    z12 = vars->z12;z13 = vars->z13;z14 = vars->z14;z15 = vars->z15;z16 = vars->z16;
    z18 = vars->z18;z20 = vars->z20;z21 = vars->z21;z22 = vars->z22;z24 = vars->z24;
    z30 = vars->z30;
    COMPLEX16 fEpn1, fEpn3, fEpn4;
    if (dm)
    {
        fEpn1 = ((-2*I)*y3 + (2*I)*y*(2 + 3*z12) - 6*y2*z6 + 2*(-1 + z12)*z6) / (9.*z6);
        fEpn3 = (2*eta*x2*((25*I)*dm*y - 81*dm*z10 + 42*dm*z6 + 18*((-I)*lchis + nchis)*(y + (4*I)*(-1 + z4)*z6))) / (81.*z6*dm) -
            (I*x2*((-36*I)*(-2 + eta)*y4*z18 + y*(36*(1 + eta) + 218*z12 + 216*z18 + 9*(-19 + 23*eta)*z24 + 9*(5 - 7*eta)*z20*z16 - 144*z6) + (3*I)*z18*(-148 + (9 + 42*eta)*z12 + 2*(-1 + 5*eta)*z24 + 189*z4 - 48*z6) + 3*y5*(3*(-3 + eta) + (-11 + eta)*z12 + 24*z6) - (6*I)*y2*z18*(15 - 18*eta + (-11 + eta)*z12 + 24*z6) - 6*y3*(6*(-1 + eta) + z6*(12 + (11 - 4*eta + 2*(-1 + 5*eta)*z12 - 36*z6)*z6)))) / (162.*z18);
        fEpn4 = (x3*(dm*echis*(I*(-2 + eta)*y3 + I*y*(4 - 2*eta + (-24 + 17*eta)*z12) + 2*(3 - 2*eta)*y2*z6 - (-2 + 3*eta)*(-1 + z12)*z6) + echia*((2*I)*(-1 + 4*eta)*y3 + I*y*(4 - 16*eta + (-24 + 101*eta)*z12) + (6 - 25*eta)*y2*z6 - 2*(-1 + 5*eta)*(-1 + z12)*z6))) / (9.*z12*dm);
    } else {
        fEpn1 = 0.0;
        fEpn3 = (2*eta*x2*((25*I)*dm*y - 81*dm*z10 + 42*dm*z6 + 18*((-I)*lchis + nchis)*(y + (4*I)*(-1 + z4)*z6))) / (81.*z6);
        fEpn4 = (echia*x3*((2*I)*(-1 + 4*eta)*y3 + I*y*(4 - 16*eta + (-24 + 101*eta)*z12) + (6 - 25*eta)*y2*z6 - 2*(-1 + 5*eta)*(-1 + z12)*z6)) / (9.*z12);
    }
    if (PREC_FLAG > 3 || PREC_FLAG==1)
        *fEcc = 0.0;
    else
        *fEcc = fEpn1 + fEpn3 + fEpn4;
    REAL8 v2 = v*v;
    COMPLEX16 f33v2, f33v3, auxflm;
    if (dm)
    {
        f33v2 = -16.*(I*nchis + lchis)*eta/9./dm;
        f33v3 = echis*(-2. + 5.*eta/2.) + echia*(-2. + 19.*eta/2.) / dm;
    } else {
        f33v2 = -16.*(I*nchis + lchis)*eta/9.;
        f33v3 = echia*(-2. + 19.*eta/2.);
    }
    if (PREC_FLAG == 4 || PREC_FLAG == 3)
        auxflm = v2 * (f33v2 + v * (f33v3 + v * (hCoeffs->f33v4 + v * (hCoeffs->f33v5  + v * hCoeffs->f33v6)))) + _Complex_I * vh3 * vh3 * hCoeffs->f33vh6;
    else
        auxflm = v2 * (0.0 + v * (hCoeffs->f33v3 + v * (hCoeffs->f33v4 + v * (hCoeffs->f33v5  + v * hCoeffs->f33v6)))) + _Complex_I * vh3 * vh3 * hCoeffs->f33vh6;

    return auxflm;
}

static COMPLEX16 prec_CalculateWaveformFactor_mode44_v3(SEOBPrecWaveformVariables *vars,
    REAL8 eta, REAL8 dm, REAL8 v, REAL8 eulerlogxabs, FacWaveformCoeffs *hCoeffs, COMPLEX16 *fEcc)
{
    REAL8 eta2 = eta*eta;
    REAL8 v2, m1Plus3eta;
    v2 = v*v;
    m1Plus3eta = -1. + 3.*eta;
    REAL8 nchia, nchis, lchia, lchis, echia, echis;
    nchia = vars->nchia;
    nchis = vars->nchis;
    lchia = vars->lchia;
    lchis = vars->lchis;
    echia = vars->echia;
    echis = vars->echis;
    // print_debug("nchia = %.16e\n", nchia);
    REAL8 x, x2, x3, x4;
    REAL8 y, y2, y3, y4, y5, y6;
    REAL8 z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13, z14, z15, z16, z18, z20, z21, z22, z24, z30;
    x = vars->x;x2 = vars->x2;x3 = vars->x3;x4 = vars->x4;
    y = vars->y;y2 = vars->y2;y3 = vars->y3;y4 = vars->y4;y5 = vars->y5;y6 = vars->y6;
    z2 = vars->z2;z3 = vars->z3;z4 = vars->z4;z5 = vars->z5;z6 = vars->z6;
    z7 = vars->z7;z8 = vars->z8;z9 = vars->z9;z10 = vars->z10;z11 = vars->z11;
    z12 = vars->z12;z13 = vars->z13;z14 = vars->z14;z15 = vars->z15;z16 = vars->z16;
    z18 = vars->z18;z20 = vars->z20;z21 = vars->z21;z22 = vars->z22;z24 = vars->z24;
    z30 = vars->z30;
    COMPLEX16 fEpn2, fEpn4;
    fEpn2 = (7 + 6*y4 + 51*z12 - 18*y2*(1 + 2*z12) + 6*z24 - (24*I)*y3*z6 + (6*I)*y*(9 + 4*z12)*z6 - 64*z8) / (64.*z8);
    fEpn4 = (x2*(-3080*(-1 + 2*eta + 3*eta2) + 20*(1169 + 3315*eta2 - 4*eta*(1123 + 99*lchia + (99*I)*nchia - 99*dm*(lchis + I*nchis)))*z12 - 
        44880*(-1 + 3*eta)*z18 + (-90943 - 182520*eta2 + 60*eta*(6001 + 132*lchia + (132*I)*nchia - 132*dm*(lchis + I*nchis)))*z24 - 
        21120*(-1 + 3*eta)*z30 + 9*(1147 - 10740*eta + 12720*eta2)*z20*z16 + 60*(8 - 106*eta + 183*eta2)*z24*z24 + 12320*(-1 + 3*eta)*z6 - 
        60*y6*(-22*(3 - 10*eta + 3*eta2) + (-74 + 238*eta + 15*eta2)*z12 - 176*(-1 + 3*eta)*z6) + 6*y*z6*((-990*I)*(-1 + 2*eta + 3*eta2) + 
        (12869*I + (13635*I)*eta2 + 20*eta*(-1999*I - (165*I)*lchia + 165*nchia + (165*I)*dm*(lchis + I*nchis)))*z12 - (8800*I)*(-1 + 3*eta)*z18 - 
        I*(6047 - 24130*eta + 22755*eta2)*z24 + (20*I)*(49 - 278*eta + 267*eta2)*z20*z16 + (3960*I)*(-1 + 3*eta)*z6) - (120*I)*y5*z6*(11*(3 - 10*eta + 3*eta2) + 
        (115 - 410*eta + 69*eta2)*z12 + 88*(-1 + 3*eta)*z6) - 60*y4*(22*(7 - 26*eta + 15*eta2) + (-70 + 218*eta - 339*eta2)*z12 - 528*(-1 + 3*eta)*z18 + 
        (172 - 794*eta + 519*eta2)*z24 + 352*(-1 + 3*eta)*z6) + (6*I)*y3*z6*(55*(19 - 74*eta + 51*eta2) + (3803 - 15770*eta + 2625*eta2)*z12 - 
        8800*(-1 + 3*eta)*z18 + 1320*(-1 + 2*eta + 3*eta2)*z24 + 2200*(-1 + 3*eta)*z6) - 2*y2*(-110*(-15 + 2*eta + 129*eta2) + 
        20*(1114 + 1797*eta2 + eta*(-3722 - 99*lchia - (99*I)*nchia + 99*dm*(lchis + I*nchis)))*z12 + 6600*(-1 + 3*eta)*z18 - 3*(11273 - 44540*eta + 22650*eta2)*z24 + 
        10560*(-1 + 3*eta)*z30 + 30*(238 - 926*eta + 321*eta2)*z20*z16 + 9680*(-1 + 3*eta)*z6))) / (42240*(-1 + 3*eta)*z20);
    if (PREC_FLAG > 3 && PREC_FLAG==1)
        *fEcc = 0.0;
    else
        *fEcc = fEpn2 + fEpn4;
    COMPLEX16 rho44v2;
    rho44v2 = hCoeffs->rho44v2 + 81.*eta*(-I*nchia + I*nchis*dm - lchia + dm*lchis) / (256.*m1Plus3eta);
    COMPLEX16 rholm;
    if (PREC_FLAG == 4 || PREC_FLAG == 3)
        rholm = 1. + v2 * (rho44v2
            + v * (hCoeffs->rho44v3 + v * (hCoeffs->rho44v4 + v * (hCoeffs->rho44v5 + v * (hCoeffs->rho44v6 + hCoeffs->rho44v6l *
            eulerlogxabs + v2 *( hCoeffs->rho44v8 + hCoeffs->rho44v8l * eulerlogxabs +v2 * (hCoeffs->rho44v10 + hCoeffs->rho44v10l * eulerlogxabs) ) )))));
    else
        rholm = 1. + v2 * (hCoeffs->rho44v2
            + v * (hCoeffs->rho44v3 + v * (hCoeffs->rho44v4 + v * (hCoeffs->rho44v5 + v * (hCoeffs->rho44v6 + hCoeffs->rho44v6l *
            eulerlogxabs + v2 *( hCoeffs->rho44v8 + hCoeffs->rho44v8l * eulerlogxabs +v2 * (hCoeffs->rho44v10 + hCoeffs->rho44v10l * eulerlogxabs) ) )))));

    return rholm;
}

static COMPLEX16 prec_CalculateWaveformFactor_mode55_v3(SEOBPrecWaveformVariables *vars,
    REAL8 eta, REAL8 dm)
{
    REAL8 eta2 = eta*eta;
    REAL8 v2, m1Plus3eta;
    m1Plus3eta = -1. + 3.*eta;
    REAL8 nchia, nchis, lchia, lchis, echia, echis;
    nchia = vars->nchia;
    nchis = vars->nchis;
    lchia = vars->lchia;
    lchis = vars->lchis;
    echia = vars->echia;
    echis = vars->echis;
    REAL8 x, x2, x3, x4;
    REAL8 y, y2, y3, y4, y5, y6;
    REAL8 z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13, z14, z15, z16, z18, z20, z21, z22, z24, z30;
    x = vars->x;x2 = vars->x2;x3 = vars->x3;x4 = vars->x4;
    y = vars->y;y2 = vars->y2;y3 = vars->y3;y4 = vars->y4;y5 = vars->y5;y6 = vars->y6;
    z2 = vars->z2;z3 = vars->z3;z4 = vars->z4;z5 = vars->z5;z6 = vars->z6;
    z7 = vars->z7;z8 = vars->z8;z9 = vars->z9;z10 = vars->z10;z11 = vars->z11;
    z12 = vars->z12;z13 = vars->z13;z14 = vars->z14;z15 = vars->z15;z16 = vars->z16;
    z18 = vars->z18;z20 = vars->z20;z21 = vars->z21;z22 = vars->z22;z24 = vars->z24;
    z30 = vars->z30;
    COMPLEX16 fEpn3;
    if (PREC_FLAG > 3 && PREC_FLAG==1)
        fEpn3 = 0.0;
    else
        fEpn3 = -1. + ((24*I)*y5 - (48*I)*y3*(2 + 5*z12) + (2*I)*y*(41 + 312*z12 + 60*z24) + 120*y4*z6 - 12*y2*(33 + 20*z12)*z6 + (172 + 429*z12 + 24*z24)*z6) / (625.*z10);
    return fEpn3;
}

/**
 * This function calculates hlm mode factorized-resummed waveform
 * for given dynamical variables.
 * Eq. 17 and the entire Appendix of PRD 86, 024011 (2012) + changes
 * described in ￼the section "Factorized waveforms" of https://dcc.ligo.org/T1400476
 */
INT
prec_EOBGetPrecEccSpinFactorizedWaveform_v2(
            COMPLEX16 * hlm,	/**< OUTPUT, hlm waveforms */
            REAL8Vector * values,	/**< dyanmical variables: \f$(r,\phi,p_r,p_\phi)\f$ */
            REAL8Vector * cartvalues,	/**< dyanmical variables */
            const REAL8 v,	/**< velocity */
            const REAL8 Hreal,	/**< real Hamiltonian */
            const INT l,	/**< l mode index */
            const INT m,	/**< m mode index */
            SEOBPrecWaveformVariables *vars,
            SpinEOBParams * params	/**< Spin EOB parameters */
)
{
            // print_debug("here l,m = (%d, %d)\n", l, m);
    int		debugPK = 0;
    /* Status of function calls */
    INT status;
    INT i;

    REAL8 eta;
    REAL8 r , pp, Omega, v2, Omegav2, vh, vh3, k, hathatk, eulerlogxabs;
    //pr
    REAL8 rcrossp_x, rcrossp_y, rcrossp_z;
    REAL8 Slm, deltalm;
    COMPLEX16 auxflm = 0.0;
    REAL8 hathatksq4, hathatk4pi, Tlmprefac, Tlmprodfac;
    COMPLEX16 Tlm, rholmPwrl,rholm;
    COMPLEX16 hNewton;
    COMPLEX16 facEcc;
    gsl_sf_result	z2;

    /* Non-Keplerian velocity */
    REAL8		vPhi    , vPhi2;

    /* Pre-computed coefficients */

    FacWaveformCoeffs *hCoeffs = params->hCoeffs;

    if (abs(m) > (INT) l) 
    {
        return CEV_FAILURE;
    }
    if (m == 0) 
    {
        return CEV_FAILURE;
    }
    eta = params->eta;

    /* Check our eta was sensible */
    if (eta > 0.25) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Eta seems to be > 0.25 - this isn't allowed!");
        return CEV_FAILURE;
    }
    /*
    * else if ( eta == 0.25 && m % 2 ) { // If m is odd and dM = 0, hLM
    * will be zero memset( hlm, 0, sizeof( COMPLEX16 ) ); return
    * XLAL_SUCCESS; }
    */

    r = values->data[0];
    pp = values->data[3];

    rcrossp_x = cartvalues->data[1] * cartvalues->data[5] - cartvalues->data[2] * cartvalues->data[4];
    rcrossp_y = cartvalues->data[2] * cartvalues->data[3] - cartvalues->data[0] * cartvalues->data[5];
    rcrossp_z = cartvalues->data[0] * cartvalues->data[4] - cartvalues->data[1] * cartvalues->data[3];

    //pp = sqrt(rcrossp_x * rcrossp_x + rcrossp_y * rcrossp_y + rcrossp_z * rcrossp_z);

    v2 = v * v;
    Omega = v2 * v;
    Omegav2 = Omega*v2;
    vh3 = Hreal * Omega;
    vh = cbrt(vh3);
    eulerlogxabs = CST_GAMMA + log(2.0 * (REAL8) m * v);

    FacWaveformCoeffs oldCoeffs = *(params->hCoeffs); //RC: the function XLALSimIMRSpinPrecEOBNonKeplerCoeff is calculating again the coefficients hCoeffs, for the omega. These are different because
    // for the dynamics we are using different coefficients for the 21 mode. I store the old coefficients and the I recover them ofter vPhi is computed.

    /* Calculate the non-Keplerian velocity */
    if (params->alignedSpins) 
    {

        vPhi = XLALSimIMRSpinAlignedEOBNonKeplerCoeff(values->data, params);

        if (IS_REAL8_FAIL_NAN(vPhi)) 
        {
            return CEV_FAILURE;
        }
        vPhi = r * cbrt(vPhi);

        // if (debugPK)
        //     XLAL_PRINT_INFO("In XLALSimIMRSpinEOBFluxCalculateNewtonianMultipole, getting rW = %.12e\n", vPhi);
        vPhi *= Omega;
        vPhi2 = vPhi * vPhi;
    } 
    else 
    {
        vPhi = XLALSimIMRSpinPrecEOBNonKeplerCoeff(cartvalues->data, params);
        if (IS_REAL8_FAIL_NAN(vPhi)) 
        {
            return CEV_FAILURE;
        }
        vPhi = r * cbrt(vPhi);

        // if (debugPK)
        //     XLAL_PRINT_INFO("In XLALSimIMRSpinEOBFluxCalculateNewtonianMultipole, getting rW = %.12e\n",
        //             vPhi);
        vPhi *= Omega;
        vPhi2 = vPhi * vPhi;
    }
    *hCoeffs = oldCoeffs; //RC: Here I recover the old coefficients
    /*
        * Calculate the newtonian multipole, 1st term in Eq. 17, given by
        * Eq. A1
        */
    // if (debugPK) {
    //     XLAL_PRINT_INFO("\nValues inside XLALSimIMRSpinEOBFluxGetSpinFactorizedWaveform:\n");
    //     for (i = 0; i < 11; i++)
    //         XLAL_PRINT_INFO("values[%d] = %.12e\n", i, cartvalues->data[i]);

    //     XLAL_PRINT_INFO("Calculating hNewton, with v = %.12e, vPhi = %.12e, r = %.12e, Phi = %.12e, l = %d, m = %d\n",
    //             v, vPhi, r, values->data[1], (UINT4) l, (UINT4) m);
    // }
    status = XLALSimIMRSpinEOBCalculateNewtonianMultipole(&hNewton,
            vPhi2, r, cartvalues->data[12] + cartvalues->data[13], (UINT) l, m, params);
                                
    if (status!=CEV_SUCCESS) 
    {
        return CEV_FAILURE;
    }
    /* Calculate the source term, 2nd term in Eq. 17, given by Eq. A5, Hreal is given by Eq.5 and Heff is in Eq.2 */
    if (((l + m) % 2) == 0) 
    {
        Slm = (Hreal * Hreal - 1.) / (2. * eta) + 1.;
    } 
    else 
    {
        Slm = v * pp;
        //Slm = v * sqrt(rcrossp_x * rcrossp_x + rcrossp_y * rcrossp_y + rcrossp_z * rcrossp_z);
    }
    // if (debugPK)
    //     XLAL_PRINT_INFO("In XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform: Hreal = %e, Slm = %e, eta = %e\n", Hreal, Slm, eta);

    /*
    * Calculate the Tail term, 3rd term in Eq. 17,
    * given by Eq. A6, and Eq. (42) of
    * http://arxiv.org/pdf/1212.4357.pdf (or PRD 87 084035 (2013))
    */
    k = m * Omega;
    hathatk = Hreal * k;

    gsl_sf_result	lnr1, arg1;
    status = gsl_sf_lngamma_complex_e(l + 1.0, -2.0 * hathatk, &lnr1, &arg1);
    if (status != GSL_SUCCESS) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Error in GSL function");
        return CEV_FAILURE;
    }
    status = gsl_sf_fact_e(l, &z2);
    if (status != GSL_SUCCESS) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Error in GSL function");
        return CEV_FAILURE;
    }
    Tlm = cexp((lnr1.val + CST_PI * hathatk) + I * (
            arg1.val + 2.0 * hathatk * log(4.0 * k / sqrt(CST_E))));
    Tlm /= z2.val;


    // if (debugPK){
    //     hathatksq4 = 4. * hathatk * hathatk;
    //     hathatk4pi = 4. * LAL_PI * hathatk;
    //     /* Calculating the prefactor of Tlm, outside the multiple product */
    //     Tlmprefac = sqrt(hathatk4pi / (1. - exp(-hathatk4pi))) / z2.val;

    //     /* Calculating the multiple product factor */
    //     for (Tlmprodfac = 1., i = 1; i <= l; i++)
    //         Tlmprodfac *= (hathatksq4 + (REAL8) i * i);

    //     REAL8		Tlmold;
    //     Tlmold = Tlmprefac * sqrt(Tlmprodfac);
    //     XLAL_PRINT_INFO("Tlm = %e + i%e, |Tlm| = %.16e (should be %.16e)\n", creal(Tlm), cimag(Tlm), cabs(Tlm), Tlmold);
    // }

    /* Calculate the residue phase and amplitude terms */
    /*
        * deltalm is the 4th term in Eq. 17, delta 22 given by Eq. A15,
        * others
        */
    /* rholm is the 5th term in Eq. 17, given by Eqs. A8 - A14 */
    /*
        * auxflm is a special part of the 5th term in Eq. 17, given by Eq.
        * A15
        */
    /*
        * Actual values of the coefficients are defined in the another function
        * see file LALSimIMRSpinEOBFactorizedWaveformCoefficientsPrec.c
        */
    facEcc = 0.0;
    REAL8 dm;
    if (1. - 4.*eta < 0.0) dm = 0.0;
    else dm = sqrt(1. - 4.*eta);
    switch (l) {
    case 2:
        switch (abs(m)) {
        case 2:
            deltalm = vh3 * (hCoeffs->delta22vh3 + vh3 * (hCoeffs->delta22vh6
                    + vh * vh * (hCoeffs->delta22vh9 * vh)))
                + Omega*(hCoeffs->delta22v5 * v2 + Omega*(hCoeffs->delta22v6  + hCoeffs->delta22v8 *v2));
            // rholm = 1. + v2 * (hCoeffs->rho22v2 + v * (hCoeffs->rho22v3
            //                     + v * (hCoeffs->rho22v4
            //         + v * (hCoeffs->rho22v5 + v * (hCoeffs->rho22v6
            //                     + hCoeffs->rho22v6l * eulerlogxabs + v * (hCoeffs->rho22v7
            //                                             + v * (hCoeffs->rho22v8 + hCoeffs->rho22v8l * eulerlogxabs
            //                                                 + (hCoeffs->rho22v10 + hCoeffs->rho22v10l * eulerlogxabs) * v2)))))));
            rholm = prec_CalculateWaveformFactor_mode22_v3(vars, eta, dm, v, eulerlogxabs, hCoeffs, &facEcc);
            //FIXME
                // if (debugPK){
                //     XLAL_PRINT_INFO("PK:: rho22v2 = %.12e, rho22v3 = %.12e, rho22v4 = %.12e,\n rho22v5 = %.16e, rho22v6 = %.16e, rho22v6LOG = %.16e, \n rho22v7 = %.12e, rho22v8 = %.16e, rho22v8LOG = %.16e, \n rho22v10 = %.16e, rho22v10LOG = %.16e\n, rho22v6 = %.12e, rho22v8 = %.12e, rho22v10 = %.12e\n",
                //         hCoeffs->rho22v2, hCoeffs->rho22v3, hCoeffs->rho22v4,
                //         hCoeffs->rho22v5, hCoeffs->rho22v6, hCoeffs->rho22v6l,
                //         hCoeffs->rho22v7, hCoeffs->rho22v8, hCoeffs->rho22v8l,
                //         hCoeffs->rho22v10, hCoeffs->rho22v10l,
                //         hCoeffs->rho22v6 + hCoeffs->rho22v6l * eulerlogxabs,
                //         hCoeffs->rho22v8 + hCoeffs->rho22v8l * eulerlogxabs,
                //         hCoeffs->rho22v10 + hCoeffs->rho22v10l * eulerlogxabs);}
            break;
        case 1:
            {
                deltalm = vh3 * (hCoeffs->delta21vh3 + vh3 * (hCoeffs->delta21vh6
                                            + vh * (hCoeffs->delta21vh7 + (hCoeffs->delta21vh9) * vh * vh)))
                    + Omegav2*(hCoeffs->delta21v5  + hCoeffs->delta21v7 * v2);
                rholm = 1. + v * (hCoeffs->rho21v1
                            + v * (hCoeffs->rho21v2 + v * (hCoeffs->rho21v3 + v * (hCoeffs->rho21v4
                                                        + v * (hCoeffs->rho21v5 + v * (hCoeffs->rho21v6 + hCoeffs->rho21v6l * eulerlogxabs
                                                                    + v * (hCoeffs->rho21v7 + hCoeffs->rho21v7l * eulerlogxabs
                                                                            + v * (hCoeffs->rho21v8 + hCoeffs->rho21v8l * eulerlogxabs
                                                                                + (hCoeffs->rho21v10 + hCoeffs->rho21v10l * eulerlogxabs) * v2))))))));
                // auxflm = v * (hCoeffs->f21v1 + v2 * (hCoeffs->f21v3 + v * hCoeffs->f21v4 + v2 * (hCoeffs->f21v5 + v * hCoeffs->f21v6 + v2 * hCoeffs->f21v7c)));
                auxflm = prec_CalculateWaveformFactor_mode21_v3(vars, eta, dm, v, eulerlogxabs, hCoeffs, &facEcc);
            }
            break;
        default:
            return CEV_FAILURE;
            break;
        }
        break;
    case 3:
        switch (m) {
        case 3:
            deltalm = vh3 * (hCoeffs->delta33vh3 + vh3 * (hCoeffs->delta33vh6 + hCoeffs->delta33vh9 * vh3))
                + hCoeffs->delta33v5 * v * v2 * v2 + hCoeffs->delta33v7 * v2 * v2 * v2 * v;	
            //R.C: delta33v7 is set to 0, whoever is adding it here as a coefficient is evil, TODO: double check that is zero and then remove it
            //RC: This terms are in Eq.A6 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
            rholm = 1. + v2 * (hCoeffs->rho33v2 + v * (hCoeffs->rho33v3 + v * (hCoeffs->rho33v4 + v * (hCoeffs->rho33v5 + v * (hCoeffs->rho33v6 +
                    hCoeffs->rho33v6l * eulerlogxabs + v * (hCoeffs->rho33v7 + v * (hCoeffs->rho33v8 + hCoeffs->rho33v8l * eulerlogxabs +
                    v2*(hCoeffs->rho33v10 + hCoeffs->rho33v10l*eulerlogxabs))))))));
            //RC: This terms are in Eq.A10 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
            // auxflm = v * (v2 * (hCoeffs->f33v3 + v * (hCoeffs->f33v4 + v * (hCoeffs->f33v5  + v * hCoeffs->f33v6)))) + _Complex_I * vh3 * vh3 * hCoeffs->f33vh6;
            auxflm = prec_CalculateWaveformFactor_mode33_v3(vars, eta, dm, v, vh3, eulerlogxabs, hCoeffs, &facEcc);
            break;
        case 2:
            deltalm = vh3 * (hCoeffs->delta32vh3 + vh * (hCoeffs->delta32vh4 + vh * vh * (hCoeffs->delta32vh6
                            + hCoeffs->delta32vh9 * vh3)));
            rholm = 1. + v * (hCoeffs->rho32v
                        + v * (hCoeffs->rho32v2 + v * (hCoeffs->rho32v3 + v * (hCoeffs->rho32v4 + v * (hCoeffs->rho32v5
                                                                + v * (hCoeffs->rho32v6 + hCoeffs->rho32v6l * eulerlogxabs
                                                                + (hCoeffs->rho32v8 + hCoeffs->rho32v8l * eulerlogxabs) * v2))))));
            break;
        case 1:
            deltalm = vh3 * (hCoeffs->delta31vh3 + vh3 * (hCoeffs->delta31vh6
                                        + vh * (hCoeffs->delta31vh7 + hCoeffs->delta31vh9 * vh * vh)))
                + hCoeffs->delta31v5 * v * v2 * v2;
            rholm = 1. + v2 * (hCoeffs->rho31v2 + v * (hCoeffs->rho31v3 + v * (hCoeffs->rho31v4
                                                + v * (hCoeffs->rho31v5 + v * (hCoeffs->rho31v6 + hCoeffs->rho31v6l * eulerlogxabs
                                                                + v * (hCoeffs->rho31v7 + (hCoeffs->rho31v8 + hCoeffs->rho31v8l * eulerlogxabs) * v))))));
            auxflm = v * v2 * hCoeffs->f31v3;
            break;
        default:
            return CEV_FAILURE;
            break;
        }
        break;
    case 4:
        switch (m) {
        case 4:
            //RC: This terms are in Eq.A15 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
            deltalm = vh3 * (hCoeffs->delta44vh3 + vh3 * (hCoeffs->delta44vh6 + vh3 * hCoeffs->delta44vh9))
                    + hCoeffs->delta44v5 * v2 * v2 * v;
            //RC: This terms are in Eq.A8 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
            // rholm = 1. + v2 * (hCoeffs->rho44v2
            //         + v * (hCoeffs->rho44v3 + v * (hCoeffs->rho44v4 + v * (hCoeffs->rho44v5 + v * (hCoeffs->rho44v6 + hCoeffs->rho44v6l *
            //         eulerlogxabs + v2 *( hCoeffs->rho44v8 + hCoeffs->rho44v8l * eulerlogxabs +v2 * (hCoeffs->rho44v10 + hCoeffs->rho44v10l * eulerlogxabs) ) )))));
            rholm = prec_CalculateWaveformFactor_mode44_v3(vars, eta, dm, v, eulerlogxabs, hCoeffs, &facEcc);
            // facEcc = prec_CalculateWaveformFactor_mode44(vars, eta, dm);
            break;
        case 3:
            deltalm = vh3 * (hCoeffs->delta43vh3 + vh * (hCoeffs->delta43vh4
                        + hCoeffs->delta43vh6 * vh * vh));
            rholm = 1. + v * (hCoeffs->rho43v
                        + v * (hCoeffs->rho43v2
                + v2 * (hCoeffs->rho43v4 + v * (hCoeffs->rho43v5
                                + (hCoeffs->rho43v6 + hCoeffs->rho43v6l * eulerlogxabs) * v))));
            auxflm = v * hCoeffs->f43v;
            break;
        case 2:
            deltalm = vh3 * (hCoeffs->delta42vh3 + hCoeffs->delta42vh6 * vh3);
            rholm = 1. + v2 * (hCoeffs->rho42v2
                        + v * (hCoeffs->rho42v3 + v * (hCoeffs->rho42v4 + v * (hCoeffs->rho42v5
                                                    + (hCoeffs->rho42v6 + hCoeffs->rho42v6l * eulerlogxabs) * v))));
            break;
        case 1:
            deltalm = vh3 * (hCoeffs->delta41vh3 + vh * (hCoeffs->delta41vh4
                        + hCoeffs->delta41vh6 * vh * vh));
            rholm = 1. + v * (hCoeffs->rho41v
                        + v * (hCoeffs->rho41v2
                + v2 * (hCoeffs->rho41v4 + v * (hCoeffs->rho41v5
                                + (hCoeffs->rho41v6 + hCoeffs->rho41v6l * eulerlogxabs) * v))));
            auxflm = v * hCoeffs->f41v;
            break;
        default:
            return CEV_FAILURE;
            break;
        }
        break;
    case 5:
        switch (m) {
        case 5:
            //RC: This terms are in Eq.A16 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
            deltalm = vh3 *(hCoeffs->delta55vh3 +vh3*(hCoeffs->delta55vh6 +vh3 *(hCoeffs->delta55vh9)))
                    + hCoeffs->delta55v5 * v2 * v2 * v;
            //RC: This terms are in Eq.A9 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
            rholm = 1. + v2 * (hCoeffs->rho55v2 + v * (hCoeffs->rho55v3 + v * (hCoeffs->rho55v4 + v * (hCoeffs->rho55v5 +
                    v * (hCoeffs->rho55v6 + hCoeffs->rho55v6l * eulerlogxabs +
                    v2 * (hCoeffs->rho55v8 + hCoeffs->rho55v8l * eulerlogxabs +
                    v2 * (hCoeffs->rho55v10 + hCoeffs->rho55v10l * eulerlogxabs )))))));
        //RC: This terms are in Eq.A12 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
            auxflm = v2 * v * (hCoeffs->f55v3 + v * (hCoeffs->f55v4 + v * (hCoeffs->f55v5c) ));
            facEcc = prec_CalculateWaveformFactor_mode55_v3(vars, eta, dm);
            break;
        case 4:
            deltalm = vh3 * (hCoeffs->delta54vh3 + hCoeffs->delta54vh4 * vh);
            rholm = 1. + v2 * (hCoeffs->rho54v2 + v * (hCoeffs->rho54v3
                            + hCoeffs->rho54v4 * v));
            break;
        case 3:
            deltalm = hCoeffs->delta53vh3 * vh3;
            rholm = 1. + v2 * (hCoeffs->rho53v2
                        + v * (hCoeffs->rho53v3 + v * (hCoeffs->rho53v4 + hCoeffs->rho53v5 * v)));
            break;
        case 2:
            deltalm = vh3 * (hCoeffs->delta52vh3 + hCoeffs->delta52vh4 * vh);
            rholm = 1. + v2 * (hCoeffs->rho52v2 + v * (hCoeffs->rho52v3
                            + hCoeffs->rho52v4 * v));
            break;
        case 1:
            deltalm = hCoeffs->delta51vh3 * vh3;
            rholm = 1. + v2 * (hCoeffs->rho51v2
                        + v * (hCoeffs->rho51v3 + v * (hCoeffs->rho51v4 + hCoeffs->rho51v5 * v)));
            break;
        default:
            return CEV_FAILURE;
            break;
        }
        break;
    case 6:
        switch (m) {
        case 6:
            deltalm = hCoeffs->delta66vh3 * vh3;
            rholm = 1. + v2 * (hCoeffs->rho66v2 + v * (hCoeffs->rho66v3
                            + hCoeffs->rho66v4 * v));
            break;
        case 5:
            deltalm = hCoeffs->delta65vh3 * vh3;
            rholm = 1. + v2 * (hCoeffs->rho65v2 + hCoeffs->rho65v3 * v);
            break;
        case 4:
            deltalm = hCoeffs->delta64vh3 * vh3;
            rholm = 1. + v2 * (hCoeffs->rho64v2 + v * (hCoeffs->rho64v3
                            + hCoeffs->rho64v4 * v));
            break;
        case 3:
            deltalm = hCoeffs->delta63vh3 * vh3;
            rholm = 1. + v2 * (hCoeffs->rho63v2 + hCoeffs->rho63v3 * v);
            break;
        case 2:
            deltalm = hCoeffs->delta62vh3 * vh3;
            rholm = 1. + v2 * (hCoeffs->rho62v2 + v * (hCoeffs->rho62v3
                            + hCoeffs->rho62v4 * v));
            break;
        case 1:
            deltalm = hCoeffs->delta61vh3 * vh3;
            rholm = 1. + v2 * (hCoeffs->rho61v2 + hCoeffs->rho61v3 * v);
            break;
        default:
            return CEV_FAILURE;
            break;
        }
        break;
    case 7:
        switch (m) {
        case 7:
            deltalm = hCoeffs->delta77vh3 * vh3;
            rholm = 1. + v2 * (hCoeffs->rho77v2 + hCoeffs->rho77v3 * v);
            break;
        case 6:
            deltalm = 0.0;
            rholm = 1. + hCoeffs->rho76v2 * v2;
            break;
        case 5:
            deltalm = hCoeffs->delta75vh3 * vh3;
            rholm = 1. + v2 * (hCoeffs->rho75v2 + hCoeffs->rho75v3 * v);
            break;
        case 4:
            deltalm = 0.0;
            rholm = 1. + hCoeffs->rho74v2 * v2;
            break;
        case 3:
            deltalm = hCoeffs->delta73vh3 * vh3;
            rholm = 1. + v2 * (hCoeffs->rho73v2 + hCoeffs->rho73v3 * v);
            break;
        case 2:
            deltalm = 0.0;
            rholm = 1. + hCoeffs->rho72v2 * v2;
            break;
        case 1:
            deltalm = hCoeffs->delta71vh3 * vh3;
            rholm = 1. + v2 * (hCoeffs->rho71v2 + hCoeffs->rho71v3 * v);
            break;
        default:
            return CEV_FAILURE;
            break;
        }
        break;
    case 8:
        switch (m) {
        case 8:
            deltalm = 0.0;
            rholm = 1. + hCoeffs->rho88v2 * v2;
            break;
        case 7:
            deltalm = 0.0;
            rholm = 1. + hCoeffs->rho87v2 * v2;
            break;
        case 6:
            deltalm = 0.0;
            rholm = 1. + hCoeffs->rho86v2 * v2;
            break;
        case 5:
            deltalm = 0.0;
            rholm = 1. + hCoeffs->rho85v2 * v2;
            break;
        case 4:
            deltalm = 0.0;
            rholm = 1. + hCoeffs->rho84v2 * v2;
            break;
        case 3:
            deltalm = 0.0;
            rholm = 1. + hCoeffs->rho83v2 * v2;
            break;
        case 2:
            deltalm = 0.0;
            rholm = 1. + hCoeffs->rho82v2 * v2;
            break;
        case 1:
            deltalm = 0.0;
            rholm = 1. + hCoeffs->rho81v2 * v2;
            break;
        default:
            return CEV_FAILURE;
            break;
        }
        break;
    default:
        return CEV_FAILURE;
        break;
    }
        // if (isnan(creal(facEcc)) || isnan(cimag(facEcc)))
        //     print_debug("(%d,%d)this is nan\n", l, m);
    //debugPK
    //     if (debugPK)
    //     XLAL_PRINT_INFO("rho_%d_%d = %.12e + %.12e \n", l, m, creal(rholm),cimag(rholm));
    // if (debugPK)
    //     XLAL_PRINT_INFO("exp(delta_%d_%d) = %.16e + i%.16e (abs=%e)\n", l, m, creal(cexp(I * deltalm)),
    //             cimag(cexp(I * deltalm)), cabs(cexp(I * deltalm)));
    /* Raise rholm to the lth power */
    rholmPwrl = 1.0;
    i = l;
    while (i--) {
        rholmPwrl *= rholm;
    }
    /*
        * In the equal-mass odd m case, there is no contribution from
        * nonspin terms,  and the only contribution comes from the auxflm
        * term that is proportional to chiA (asymmetric spins). In this
        * case, we must ignore the nonspin terms directly, since the leading
        * term defined by CalculateThisMultipolePrefix in
        * LALSimIMREOBNewtonianMultipole.c is not zero (see comments there).
        */
    if (eta == 0.25 && m % 2) {
        rholmPwrl = auxflm + facEcc;
    } else {
        rholmPwrl += auxflm + facEcc;
    }
    // if (m==1) {
    // 	printf("f21v1 = %.16f f21v3 = %.16f f21v4 = %.16f f21v5 = %.16f\n", hCoeffs->f21v1, hCoeffs->f21v3, hCoeffs->f21v4, hCoeffs->f21v5);
    // }

    // if (r > 0.0 && debugPK) {
    //     XLAL_PRINT_INFO("YP::dynamics variables in waveform: %i, %i, r = %.12e, v = %.12e\n", l, m, r, v);
    //     XLAL_PRINT_INFO("rholm^l = %.16e + %.16e, Tlm = %.16e + i %.16e, \nSlm = %.16e, hNewton = %.16e + i %.16e, delta = %.16e\n", creal(rholmPwrl), cimag(rholmPwrl), creal(Tlm), cimag(Tlm), Slm, creal(hNewton), cimag(hNewton), 0.0);
    // }
    /* Put all factors in Eq. 17 together */
    *hlm = Tlm * cexp(I * deltalm) * Slm * rholmPwrl;
    *hlm *= hNewton;
    // if (r > 8.5 && debugPK) {
    //     XLAL_PRINT_INFO("YP::FullWave: %.16e,%.16e, %.16e\n", creal(*hlm), cimag(*hlm), cabs(*hlm));
    // }
    return CEV_SUCCESS;
}


INT prec_CalculateFactorizedWaveformCorrection(
    COMPLEX16 * rholm,
    COMPLEX16 * flm,
    REAL8Vector *values,
    REAL8Vector *cartvalues,
    const REAL8 v,
    const REAL8 Hreal,
    const INT l,
    const INT m,
    SEOBPrecWaveformVariables *vars,
    SpinEOBParams *params
)
{
    /* Status of function calls */
    INT status;
    INT i;

    REAL8 eta;
    REAL8 r , pp, Omega, v2, Omegav2, vh, vh3, k, hathatk, eulerlogxabs;
    //pr
    REAL8 rcrossp_x, rcrossp_y, rcrossp_z;
    REAL8 Slm, deltalm;
    COMPLEX16 auxflm = 0.0;
    REAL8 hathatksq4, hathatk4pi, Tlmprefac, Tlmprodfac;
    COMPLEX16 Tlm, rholmPwrl;
    COMPLEX16 hNewton;
    COMPLEX16 facEcc = 0.0;
    gsl_sf_result	z2;

    /* Non-Keplerian velocity */
    REAL8		vPhi    , vPhi2;

    /* Pre-computed coefficients */

    FacWaveformCoeffs *hCoeffs = params->hCoeffs;

    if (abs(m) > (INT) l) 
    {
        return CEV_FAILURE;
    }
    if (m == 0) 
    {
        return CEV_FAILURE;
    }
    eta = params->eta;

    /* Check our eta was sensible */
    if (eta > 0.25) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Eta seems to be > 0.25 - this isn't allowed!");
        return CEV_FAILURE;
    }
    /*
    * else if ( eta == 0.25 && m % 2 ) { // If m is odd and dM = 0, hLM
    * will be zero memset( hlm, 0, sizeof( COMPLEX16 ) ); return
    * XLAL_SUCCESS; }
    */

    r = values->data[0];
    pp = values->data[3];

    rcrossp_x = cartvalues->data[1] * cartvalues->data[5] - cartvalues->data[2] * cartvalues->data[4];
    rcrossp_y = cartvalues->data[2] * cartvalues->data[3] - cartvalues->data[0] * cartvalues->data[5];
    rcrossp_z = cartvalues->data[0] * cartvalues->data[4] - cartvalues->data[1] * cartvalues->data[3];

    //pp = sqrt(rcrossp_x * rcrossp_x + rcrossp_y * rcrossp_y + rcrossp_z * rcrossp_z);

    v2 = v * v;
    Omega = v2 * v;
    Omegav2 = Omega*v2;
    vh3 = Hreal * Omega;
    vh = cbrt(vh3);
    eulerlogxabs = CST_GAMMA + log(2.0 * (REAL8) m * v);

    FacWaveformCoeffs oldCoeffs = *(params->hCoeffs); //RC: the function XLALSimIMRSpinPrecEOBNonKeplerCoeff is calculating again the coefficients hCoeffs, for the omega. These are different because
    /* Calculate the non-Keplerian velocity */
    if (params->alignedSpins) 
    {

        vPhi = XLALSimIMRSpinAlignedEOBNonKeplerCoeff(values->data, params);

        if (IS_REAL8_FAIL_NAN(vPhi)) 
        {
            return CEV_FAILURE;
        }
        vPhi = r * cbrt(vPhi);

        // if (debugPK)
        //     XLAL_PRINT_INFO("In XLALSimIMRSpinEOBFluxCalculateNewtonianMultipole, getting rW = %.12e\n", vPhi);
        vPhi *= Omega;
        vPhi2 = vPhi * vPhi;
    } 
    else 
    {
        vPhi = XLALSimIMRSpinPrecEOBNonKeplerCoeff(cartvalues->data, params);
        if (IS_REAL8_FAIL_NAN(vPhi)) 
        {
            return CEV_FAILURE;
        }
        vPhi = r * cbrt(vPhi);

        // if (debugPK)
        //     XLAL_PRINT_INFO("In XLALSimIMRSpinEOBFluxCalculateNewtonianMultipole, getting rW = %.12e\n",
        //             vPhi);
        vPhi *= Omega;
        vPhi2 = vPhi * vPhi;
    }
    *hCoeffs = oldCoeffs; //RC: Here I recover the old coefficients
    REAL8 dm;
    if (1. - 4.*eta < 0.0) dm = 0.0;
    else dm = sqrt(1. - 4.*eta);
    COMPLEX16 val_rholm;
    switch (l) {
    case 2:
        switch (abs(m)) {
        case 2:
            val_rholm = prec_CalculateWaveformFactor_mode22_v3(vars, eta, dm, v, eulerlogxabs, hCoeffs, &facEcc);
            break;
        case 1:
            {
                val_rholm = 1. + v * (hCoeffs->rho21v1
                            + v * (hCoeffs->rho21v2 + v * (hCoeffs->rho21v3 + v * (hCoeffs->rho21v4
                                                        + v * (hCoeffs->rho21v5 + v * (hCoeffs->rho21v6 + hCoeffs->rho21v6l * eulerlogxabs
                                                                    + v * (hCoeffs->rho21v7 + hCoeffs->rho21v7l * eulerlogxabs
                                                                            + v * (hCoeffs->rho21v8 + hCoeffs->rho21v8l * eulerlogxabs
                                                                                + (hCoeffs->rho21v10 + hCoeffs->rho21v10l * eulerlogxabs) * v2))))))));
                auxflm = prec_CalculateWaveformFactor_mode21_v3(vars, eta, dm, v, eulerlogxabs, hCoeffs, &facEcc);
            }
            break;
        default:
            return CEV_FAILURE;
            break;
        }
        break;
    case 3:
        switch (m) {
        case 3:
            val_rholm = 1. + v2 * (hCoeffs->rho33v2 + v * (hCoeffs->rho33v3 + v * (hCoeffs->rho33v4 + v * (hCoeffs->rho33v5 + v * (hCoeffs->rho33v6 +
                    hCoeffs->rho33v6l * eulerlogxabs + v * (hCoeffs->rho33v7 + v * (hCoeffs->rho33v8 + hCoeffs->rho33v8l * eulerlogxabs +
                    v2*(hCoeffs->rho33v10 + hCoeffs->rho33v10l*eulerlogxabs))))))));
            auxflm = prec_CalculateWaveformFactor_mode33_v3(vars, eta, dm, v, vh3, eulerlogxabs, hCoeffs, &facEcc);
            break;
        case 2:
            val_rholm = 1. + v * (hCoeffs->rho32v
                        + v * (hCoeffs->rho32v2 + v * (hCoeffs->rho32v3 + v * (hCoeffs->rho32v4 + v * (hCoeffs->rho32v5
                                                                + v * (hCoeffs->rho32v6 + hCoeffs->rho32v6l * eulerlogxabs
                                                                + (hCoeffs->rho32v8 + hCoeffs->rho32v8l * eulerlogxabs) * v2))))));
            break;
        case 1:
            val_rholm = 1. + v2 * (hCoeffs->rho31v2 + v * (hCoeffs->rho31v3 + v * (hCoeffs->rho31v4
                                                + v * (hCoeffs->rho31v5 + v * (hCoeffs->rho31v6 + hCoeffs->rho31v6l * eulerlogxabs
                                                                + v * (hCoeffs->rho31v7 + (hCoeffs->rho31v8 + hCoeffs->rho31v8l * eulerlogxabs) * v))))));
            auxflm = v * v2 * hCoeffs->f31v3;
            break;
        default:
            return CEV_FAILURE;
            break;
        }
        break;
    case 4:
        switch (m) {
        case 4:
            val_rholm = prec_CalculateWaveformFactor_mode44_v3(vars, eta, dm, v, eulerlogxabs, hCoeffs, &facEcc);
            break;
        case 3:
            val_rholm = 1. + v * (hCoeffs->rho43v
                        + v * (hCoeffs->rho43v2
                + v2 * (hCoeffs->rho43v4 + v * (hCoeffs->rho43v5
                                + (hCoeffs->rho43v6 + hCoeffs->rho43v6l * eulerlogxabs) * v))));
            auxflm = v * hCoeffs->f43v;
            break;
        case 2:
            val_rholm = 1. + v2 * (hCoeffs->rho42v2
                        + v * (hCoeffs->rho42v3 + v * (hCoeffs->rho42v4 + v * (hCoeffs->rho42v5
                                                    + (hCoeffs->rho42v6 + hCoeffs->rho42v6l * eulerlogxabs) * v))));
            break;
        case 1:
            val_rholm = 1. + v * (hCoeffs->rho41v
                        + v * (hCoeffs->rho41v2
                + v2 * (hCoeffs->rho41v4 + v * (hCoeffs->rho41v5
                                + (hCoeffs->rho41v6 + hCoeffs->rho41v6l * eulerlogxabs) * v))));
            auxflm = v * hCoeffs->f41v;
            break;
        default:
            return CEV_FAILURE;
            break;
        }
        break;
    case 5:
        switch (m) {
        case 5:
            val_rholm = 1. + v2 * (hCoeffs->rho55v2 + v * (hCoeffs->rho55v3 + v * (hCoeffs->rho55v4 + v * (hCoeffs->rho55v5 +
                    v * (hCoeffs->rho55v6 + hCoeffs->rho55v6l * eulerlogxabs +
                    v2 * (hCoeffs->rho55v8 + hCoeffs->rho55v8l * eulerlogxabs +
                    v2 * (hCoeffs->rho55v10 + hCoeffs->rho55v10l * eulerlogxabs )))))));
            auxflm = v2 * v * (hCoeffs->f55v3 + v * (hCoeffs->f55v4 + v * (hCoeffs->f55v5c) ));
            facEcc = prec_CalculateWaveformFactor_mode55_v3(vars, eta, dm);
            break;
        case 4:
            val_rholm = 1. + v2 * (hCoeffs->rho54v2 + v * (hCoeffs->rho54v3
                            + hCoeffs->rho54v4 * v));
            break;
        case 3:
            val_rholm = 1. + v2 * (hCoeffs->rho53v2
                        + v * (hCoeffs->rho53v3 + v * (hCoeffs->rho53v4 + hCoeffs->rho53v5 * v)));
            break;
        case 2:
            val_rholm = 1. + v2 * (hCoeffs->rho52v2 + v * (hCoeffs->rho52v3
                            + hCoeffs->rho52v4 * v));
            break;
        case 1:
            val_rholm = 1. + v2 * (hCoeffs->rho51v2
                        + v * (hCoeffs->rho51v3 + v * (hCoeffs->rho51v4 + hCoeffs->rho51v5 * v)));
            break;
        default:
            return CEV_FAILURE;
            break;
        }
        break;
    case 6:
        switch (m) {
        case 6:
            val_rholm = 1. + v2 * (hCoeffs->rho66v2 + v * (hCoeffs->rho66v3
                            + hCoeffs->rho66v4 * v));
            break;
        case 5:
            val_rholm = 1. + v2 * (hCoeffs->rho65v2 + hCoeffs->rho65v3 * v);
            break;
        case 4:
            val_rholm = 1. + v2 * (hCoeffs->rho64v2 + v * (hCoeffs->rho64v3
                            + hCoeffs->rho64v4 * v));
            break;
        case 3:
            val_rholm = 1. + v2 * (hCoeffs->rho63v2 + hCoeffs->rho63v3 * v);
            break;
        case 2:
            val_rholm = 1. + v2 * (hCoeffs->rho62v2 + v * (hCoeffs->rho62v3
                            + hCoeffs->rho62v4 * v));
            break;
        case 1:
            val_rholm = 1. + v2 * (hCoeffs->rho61v2 + hCoeffs->rho61v3 * v);
            break;
        default:
            return CEV_FAILURE;
            break;
        }
        break;
    case 7:
        switch (m) {
        case 7:
            val_rholm = 1. + v2 * (hCoeffs->rho77v2 + hCoeffs->rho77v3 * v);
            break;
        case 6:
            val_rholm = 1. + hCoeffs->rho76v2 * v2;
            break;
        case 5:
            val_rholm = 1. + v2 * (hCoeffs->rho75v2 + hCoeffs->rho75v3 * v);
            break;
        case 4:
            val_rholm = 1. + hCoeffs->rho74v2 * v2;
            break;
        case 3:
            val_rholm = 1. + v2 * (hCoeffs->rho73v2 + hCoeffs->rho73v3 * v);
            break;
        case 2:
            val_rholm = 1. + hCoeffs->rho72v2 * v2;
            break;
        case 1:
            val_rholm = 1. + v2 * (hCoeffs->rho71v2 + hCoeffs->rho71v3 * v);
            break;
        default:
            return CEV_FAILURE;
            break;
        }
        break;
    case 8:
        switch (m) {
        case 8:
            val_rholm = 1. + hCoeffs->rho88v2 * v2;
            break;
        case 7:
            val_rholm = 1. + hCoeffs->rho87v2 * v2;
            break;
        case 6:
            val_rholm = 1. + hCoeffs->rho86v2 * v2;
            break;
        case 5:
            val_rholm = 1. + hCoeffs->rho85v2 * v2;
            break;
        case 4:
            val_rholm = 1. + hCoeffs->rho84v2 * v2;
            break;
        case 3:
            val_rholm = 1. + hCoeffs->rho83v2 * v2;
            break;
        case 2:
            val_rholm = 1. + hCoeffs->rho82v2 * v2;
            break;
        case 1:
            val_rholm = 1. + hCoeffs->rho81v2 * v2;
            break;
        default:
            return CEV_FAILURE;
            break;
        }
        break;
    default:
        return CEV_FAILURE;
        break;
    }

    //debugPK
    //     if (debugPK)
    //     XLAL_PRINT_INFO("rho_%d_%d = %.12e + %.12e \n", l, m, creal(rholm),cimag(rholm));
    // if (debugPK)
    //     XLAL_PRINT_INFO("exp(delta_%d_%d) = %.16e + i%.16e (abs=%e)\n", l, m, creal(cexp(I * deltalm)),
    //             cimag(cexp(I * deltalm)), cabs(cexp(I * deltalm)));
    /* Raise rholm to the lth power */
    rholmPwrl = 1.0;
    i = l;
    while (i--) {
        rholmPwrl *= val_rholm;
    }
    /*
    * In the equal-mass odd m case, there is no contribution from
    * nonspin terms,  and the only contribution comes from the auxflm
    * term that is proportional to chiA (asymmetric spins). In this
    * case, we must ignore the nonspin terms directly, since the leading
    * term defined by CalculateThisMultipolePrefix in
    * LALSimIMREOBNewtonianMultipole.c is not zero (see comments there).
    */
    if (eta == 0.25 && m % 2) {
        rholmPwrl = auxflm + facEcc;
    } else {
        rholmPwrl += auxflm + facEcc;
    }
    *rholm = val_rholm;
    *flm = rholmPwrl;
    // print_debug("(%d, %d), rhopwrl - flm = %.5e + i%.5e\n, facEcc = %.5e + i%.5e, auxflm = %.5e + i%.5e", l, m,
    //     creal(rholmPwrl - val_rholm*val_rholm), cimag(rholmPwrl - val_rholm*val_rholm),
    //     creal(facEcc), cimag(facEcc),
    //     creal(auxflm), cimag(auxflm));
    return CEV_SUCCESS;
}