/**
 * Writer: Xiaolin.liu
 * xiaolin.liu@mail.bnu.edu.cn
 *
 * This module contains basic functions for  calculation.
 * Functions list:
 * Kernel:
 * 20xx.xx.xx, LOC
 **/

#include "newFactorizedWaveform.h"
#include "myLog.h"
#include "pFactorizedWaveform.h"
#include "pHamiltonian.h"
#include <gsl/gsl_sf_gamma.h>

INT CalculateFacWaveformAmpResV2(FacWaveformCoeffs *const hCoeffs, const REAL8 eta, const INT l, const INT m,
                                 const REAL8 v, const REAL8 vPhi, const REAL8 vh3, const REAL8 eulerlogxabs,
                                 const REAL8 x0, const REAL8 x1, const REAL8 x2, COMPLEX16 *facAmpRes)
{
    int i, Use_V1 = 0;
    REAL8 f1, f2;
    REAL8 ff00, ff10, ff01, ff11, ff20, ff02, ff40, ff31, ff22, ff13, ff04;
    REAL8 ff60, ff51, ff42, ff33, ff24, ff15, ff06, ff30, ff21, ff12, ff03;
    REAL8 ff50, ff41, ff32, ff23, ff14, ff05;
    COMPLEX16 rholm, rholmPwrl;
    REAL8 rhoHigh = 0.0, rhoLow;
    COMPLEX16 A0lm = 1.0, A1lm = 0.0, A2lm = 0.0, A3lm = 0.0, A4lm = 0.0;
    COMPLEX16 Flm0 = 0.0, Flm1 = 0.0, Flm3 = 0.0;
    REAL8 v2, vPhi2;
    COMPLEX16 rho2PN;
    COMPLEX16 auxflm = 0.0;
    REAL8 g2;
    REAL8 x0_2, x0_3, x0_4;
    x0_2 = x0 * x0;
    x0_3 = x0_2 * x0;
    x0_4 = x0_2 * x0_2;
    v2 = v * v;
    vPhi2 = vPhi * vPhi;
    f1 = x1 / x0;
    f2 = x2 / x0;
    ff00 = 1;
    ff10 = f1;
    ff01 = f2;
    ff11 = f1 * f2;
    ff20 = f1 * f1;
    ff02 = f2 * f2;
    ff30 = ff20 * f1;
    ff21 = ff20 * f2;
    ff12 = ff02 * f1;
    ff03 = ff02 * f2;
    ff40 = ff20 * ff20;
    ff31 = ff20 * ff11;
    ff22 = ff11 * ff11;
    ff13 = ff11 * ff02;
    ff04 = ff02 * ff02;
    ff50 = ff40 * f1;
    ff41 = ff31 * f1;
    ff32 = ff22 * f1;
    ff23 = ff22 * f2;
    ff14 = ff13 * f2;
    ff05 = ff04 * f2;
    ff60 = ff40 * ff20;
    ff51 = ff40 * ff11;
    ff42 = ff40 * ff02;
    ff33 = ff22 * ff11;
    ff24 = ff20 * ff04;
    ff15 = ff11 * ff04;
    ff06 = ff02 * ff04;
    g2 = pow(f2, -4. / 3.);
    switch (l)
    {
    case 2:
        switch (ABS(m))
        {
        case 2:
            rhoHigh =
                v2 * v2 * v *
                (hCoeffs->rho22v5 +
                 v * (hCoeffs->rho22v6 + hCoeffs->rho22v6l * eulerlogxabs +
                      v * (hCoeffs->rho22v7 + v * (hCoeffs->rho22v8 + hCoeffs->rho22v8l * eulerlogxabs +
                                                   v2 * (hCoeffs->rho22v10 + hCoeffs->rho22v10l * eulerlogxabs)))));
            rhoLow = 1. + v2 * (hCoeffs->rho22v2 + v * (hCoeffs->rho22v3 + v * (hCoeffs->rho22v4)));

            A0lm = hCoeffs->h22T0ff00 * ff00 + I * hCoeffs->h22T0ff11 * ff11 + hCoeffs->h22T0ff02 * ff02 +
                   hCoeffs->h22T0ff20 * ff20;
            A2lm = hCoeffs->h22T2ff00 * ff00 + I * hCoeffs->h22T2ff11 * ff11 + hCoeffs->h22T2ff02 * ff02 +
                   hCoeffs->h22T2ff20 * ff20 + hCoeffs->h22T2ff40 * ff40 + I * hCoeffs->h22T2ff31 * ff31 +
                   hCoeffs->h22T2ff22 * ff22 + I * hCoeffs->h22T2ff13 * ff13 + hCoeffs->h22T2ff04 * ff04;
            // A3lm = hCoeffs->h22T3ff01 * ff01 + I*hCoeffs->h22T3ff10 * ff10;
            A3lm = hCoeffs->h22T3ff01 * ff01 + I * hCoeffs->h22T3ff10 * ff10 + hCoeffs->h22T3ff03 * ff03 +
                   I * hCoeffs->h22T3ff30 * ff30 + hCoeffs->h22T3ff21 * ff21 + I * hCoeffs->h22T3ff12 * ff12;
            A4lm = hCoeffs->h22T4ff00 * ff00 + hCoeffs->h22T4ff11 * ff11 + hCoeffs->h22T4ff02 * ff02 +
                   hCoeffs->h22T4ff20 * ff20 + hCoeffs->h22T4ff40 * ff40 + I * hCoeffs->h22T4ff31 * ff31 +
                   hCoeffs->h22T4ff22 * ff22 + I * hCoeffs->h22T4ff13 * ff13 + hCoeffs->h22T4ff04 * ff04 +
                   hCoeffs->h22T4ff60 * ff60 + I * hCoeffs->h22T4ff51 * ff51 + hCoeffs->h22T4ff42 * ff42 +
                   I * hCoeffs->h22T4ff33 * ff33 + hCoeffs->h22T4ff24 * ff24 + I * hCoeffs->h22T4ff15 * ff15 +
                   hCoeffs->h22T4ff60 * ff60;
            break;
        case 1:
            rhoHigh = v2 * v2 *
                      (hCoeffs->rho21v4 +
                       v * (hCoeffs->rho21v5 +
                            v * (hCoeffs->rho21v6 + hCoeffs->rho21v6l * eulerlogxabs +
                                 v * (hCoeffs->rho21v7 + hCoeffs->rho21v7l * eulerlogxabs +
                                      v * (hCoeffs->rho21v8 + hCoeffs->rho21v8l * eulerlogxabs +
                                           (hCoeffs->rho21v10 + hCoeffs->rho21v10l * eulerlogxabs) * v2)))));
            rhoLow = 1. + v * (hCoeffs->rho21v1 + v * (hCoeffs->rho21v2 + v * (hCoeffs->rho21v3)));
            A0lm = g2 * hCoeffs->h21T0ff01 * ff01;
            A2lm = g2 * (hCoeffs->h21T2ff21 * ff21 + hCoeffs->h21T2ff12 * I * ff12 + hCoeffs->h21T2ff10 * I * ff10 +
                         hCoeffs->h21T2ff03 * ff03 + hCoeffs->h21T2ff01 * ff01);
            // (v*hCoeffs->f21v1 + v2*v*hCoeffs->f21v3)*(vPhi*vPhi2/x0_3) +
            Flm0 = v2 * v2 * (hCoeffs->f21v4 + v * (hCoeffs->f21v5 + v * (hCoeffs->f21v6))) * (vPhi * vPhi2 / x0_3);
            // Flm0 = v2*v*(hCoeffs->f21v3 +
            //        v*(hCoeffs->f21v4 +
            //     v * (hCoeffs->f21v5 +
            //     v * (hCoeffs->f21v6)))) * (vPhi*vPhi2/x0_3);

            Flm1 = g2 * hCoeffs->h21T1ff00 * ff00;
            // A1lm = Flm1;
            Flm3 = g2 * (hCoeffs->h21T3ff20 * ff20 + hCoeffs->h21T3ff11 * I * ff11 + hCoeffs->h21T3ff02 * ff02 +
                         hCoeffs->h21T3ff00 * ff00);
            // A3lm = Flm3;
            auxflm = Flm0 + Flm1 * x0 + Flm3 * x0_3;
            // auxflm = Flm0;
            // if (eta == 0.25 && m % 2)
            // {
            // auxflm = Flm1 * x0 + Flm0;
            // auxflm = Flm0 + Flm1 * x0;
            // auxflm += pow(x0, 7) * (hCoeffs->f21v7cEff00 +
            //     I*hCoeffs->f21v7cEff10*ff10 + I*hCoeffs->f21v7cEff11*ff11+
            //     hCoeffs->f21v7cEff01 * ff01 + hCoeffs->f21v7cEff02*ff02);
            auxflm += pow(x0, 12) * hCoeffs->f21v7cEff00;
            // auxflm = Flm0;
            // auxflm = g2 * hCoeffs->f21v1 * x0 +
            //          v2 * v * (hCoeffs->f21v3 +
            //          v * (hCoeffs->f21v4 +
            //          v * (hCoeffs->f21v5 +
            //          v * (hCoeffs->f21v6 +
            //          v * (hCoeffs->f21v7c))))) * (vPhi*vPhi2/x0_3);
            // print_debug("f21v7c = %.16f\n", hCoeffs->f21v7c);
            // auxflm = Flm0 * cexp( x0* (Flm1/Flm0 - x0*cpow(Flm1/Flm0,2)/2. +
            // x0_2*(cpow(Flm1/Flm0,3)/3. + Flm3/Flm0) ) );
            // }
            // else
            // {
            //     A1lm = Flm1;
            //     A3lm = Flm3;
            //     auxflm = Flm0;
            // }
            // print_debug("Flm0 = %g + i%g, Flm1 = %g + i%g, Flm3 = %g + i%g\n",
            // Flm0, Flm1, Flm3); print_debug("v = %g, vPhi = %g, x0 = %g\n", v,
            // vPhi, x0); print_debug("f21v4 = %g, f21v5 = %g, f21v6 = %g, f21v7
            // = %g\n",
            //     hCoeffs->f21v4, hCoeffs->f21v5, hCoeffs->f21v6,
            //     hCoeffs->f21v7c);
            break;
        default:
            PRINT_LOG_INFO(LOG_CRITICAL,
                           "Currently mode (%d, %d) is unsupported in new "
                           "factorized waveform.",
                           l, m);
            return CEV_FAILURE;
        }
        break;
    case 3:
        switch (ABS(m))
        {
        case 3:
            rhoHigh = v2 * v2 *
                      (hCoeffs->rho33v4 +
                       v * (hCoeffs->rho33v5 +
                            v * (hCoeffs->rho33v6 + hCoeffs->rho33v6l * eulerlogxabs +
                                 v * (hCoeffs->rho33v7 +
                                      v * (hCoeffs->rho33v8 + hCoeffs->rho33v8l * eulerlogxabs +
                                           v2 * (hCoeffs->rho33v10 + hCoeffs->rho33v10l * eulerlogxabs))))));
            rhoLow = 1. + v2 * (hCoeffs->rho33v2 + v * hCoeffs->rho33v3);
            A0lm = I * hCoeffs->h33T0ff30 * ff30 + hCoeffs->h33T0ff21 * ff21 + I * hCoeffs->h33T0ff12 * ff12 +
                   I * hCoeffs->h33T0ff10 * ff10 + hCoeffs->h33T0ff03 * ff03 + hCoeffs->h33T0ff01 * ff01;

            A2lm = hCoeffs->h33T2ff50 * I * ff50 + hCoeffs->h33T2ff41 * ff41 + hCoeffs->h33T2ff32 * I * ff32 +
                   hCoeffs->h33T2ff30 * I * ff30 + hCoeffs->h33T2ff23 * ff23 + hCoeffs->h33T2ff21 * ff21 +
                   hCoeffs->h33T2ff14 * I * ff14 + hCoeffs->h33T2ff12 * I * ff12 + hCoeffs->h33T2ff10 * I * ff10 +
                   hCoeffs->h33T2ff05 * ff05 + hCoeffs->h33T2ff03 * ff03 + hCoeffs->h33T2ff01 * ff01;
            Flm0 = (v2 * v2 * (hCoeffs->f33v4 + v * (hCoeffs->f33v5 + v * (hCoeffs->f33v6))) +
                    I * vh3 * vh3 * hCoeffs->f33vh6) *
                   (vPhi * vPhi2 / x0_3);

            Flm3 = (hCoeffs->h33T3ff40 * ff40 + hCoeffs->h33T3ff31 * I * ff31 + hCoeffs->h33T3ff22 * ff22 +
                    hCoeffs->h33T3ff13 * I * ff13 + hCoeffs->h33T3ff20 * ff20 + hCoeffs->h33T3ff11 * I * ff11 +
                    hCoeffs->h33T3ff02 * ff02 + hCoeffs->h33T3ff00 * ff00);
            if (eta == 0.25 && m % 2)
            {
                auxflm = Flm0 + Flm3 * x0_3;
            }
            else
            {
                A3lm = Flm3;
                auxflm = Flm0;
            }
            break;
        case 2:
            rhoHigh = v2 * v *
                      (hCoeffs->rho32v3 +
                       v * (hCoeffs->rho32v4 +
                            v * (hCoeffs->rho32v5 + v * (hCoeffs->rho32v6 + hCoeffs->rho32v6l * eulerlogxabs +
                                                         (hCoeffs->rho32v8 + hCoeffs->rho32v8l * eulerlogxabs) * v2))));
            rhoLow = 1. + v * (hCoeffs->rho32v + v * (hCoeffs->rho32v2));
            A0lm = g2 * (hCoeffs->h32T0ff11 * I * ff11 + hCoeffs->h32T0ff02 * ff02);
            A1lm = g2 * (hCoeffs->h32T1ff10 * I * ff10 + hCoeffs->h32T1ff01 * ff01);
            A2lm = g2 * (hCoeffs->h32T2ff31 * I * ff31 + hCoeffs->h32T2ff22 * ff22 + hCoeffs->h32T2ff13 * I * ff13 +
                         hCoeffs->h32T2ff11 * I * ff11 + hCoeffs->h32T2ff04 * ff04 + hCoeffs->h32T2ff02 * ff02);
            break;
        case 1:
            Use_V1 = 1;
            rhoHigh = v2 * v2 *
                      (hCoeffs->rho31v4 +
                       v * (hCoeffs->rho31v5 +
                            v * (hCoeffs->rho31v6 + hCoeffs->rho31v6l * eulerlogxabs +
                                 v * (hCoeffs->rho31v7 + (hCoeffs->rho31v8 + hCoeffs->rho31v8l * eulerlogxabs) * v))));
            rhoLow = 1. + v2 * (hCoeffs->rho31v2 + v * (hCoeffs->rho31v3));
            A0lm = hCoeffs->h31T0ff30 * I * ff30 + hCoeffs->h31T0ff21 * ff21 + hCoeffs->h31T0ff12 * I * ff12 +
                   hCoeffs->h31T0ff10 * I * ff10 + hCoeffs->h31T0ff03 * ff03 + hCoeffs->h31T0ff01 * ff01;
            A2lm = hCoeffs->h31T2ff50 * I * ff50 + hCoeffs->h31T2ff41 * ff41 + hCoeffs->h31T2ff32 * I * ff32 +
                   hCoeffs->h31T2ff30 * I * ff30 + hCoeffs->h31T2ff23 * ff23 + hCoeffs->h31T2ff21 * ff21 +
                   hCoeffs->h31T2ff14 * I * ff14 + hCoeffs->h31T2ff12 * I * ff12 + hCoeffs->h31T2ff10 * I * ff10 +
                   hCoeffs->h31T2ff05 * ff05 + hCoeffs->h31T2ff03 * ff03 + hCoeffs->h31T2ff01 * ff01;
            auxflm = (hCoeffs->h31T3ff20 * ff20 + hCoeffs->h31T3ff11 * I * ff11 + hCoeffs->h31T3ff02 * ff02 +
                      hCoeffs->h31T3ff00 * ff00) *
                     x0_3;
            break;
        default:
            PRINT_LOG_INFO(LOG_CRITICAL,
                           "Currently mode (%d, %d) is unsupported in new "
                           "factorized waveform.",
                           l, m);
            return CEV_FAILURE;
        }
        break;
    case 4:
        switch (ABS(m))
        {
        case 4:
            rhoHigh = v * v2 *
                      (hCoeffs->rho44v3 +
                       v * (hCoeffs->rho44v4 +
                            v * (hCoeffs->rho44v5 +
                                 v * ((hCoeffs->rho44v6 + hCoeffs->rho44v6l * eulerlogxabs) +
                                      v2 * ((hCoeffs->rho44v8 + hCoeffs->rho44v8l * eulerlogxabs) +
                                            v2 * ((hCoeffs->rho44v10 + hCoeffs->rho44v10l * eulerlogxabs)))))));
            rhoLow = 1. + v2 * hCoeffs->rho44v2;
            A0lm = hCoeffs->h44T0ff40 * ff40 + hCoeffs->h44T0ff31 * I * ff31 + hCoeffs->h44T0ff22 * ff22 +
                   hCoeffs->h44T0ff20 * ff20 + hCoeffs->h44T0ff13 * I * ff13 + hCoeffs->h44T0ff11 * I * ff11 +
                   hCoeffs->h44T0ff04 * ff04 + hCoeffs->h44T0ff02 * ff02 + hCoeffs->h44T0ff00 * ff00;
            A2lm = hCoeffs->h44T2ff60 * ff60 + hCoeffs->h44T2ff51 * I * ff51 + hCoeffs->h44T2ff42 * ff42 +
                   hCoeffs->h44T2ff40 * ff40 + hCoeffs->h44T2ff31 * I * ff31 + hCoeffs->h44T2ff24 * ff24 +
                   hCoeffs->h44T2ff22 * ff22 + hCoeffs->h44T2ff20 * ff20 + hCoeffs->h44T2ff15 * I * ff15 +
                   hCoeffs->h44T2ff13 * I * ff13 + hCoeffs->h44T2ff11 * I * ff11 + hCoeffs->h44T2ff06 * ff06 +
                   hCoeffs->h44T2ff04 * ff04 + hCoeffs->h44T2ff02 * ff02 + hCoeffs->h44T2ff00 * ff00;

            break;
        case 3:
            rhoHigh = v2 * (hCoeffs->rho43v2 +
                            v2 * (hCoeffs->rho43v4 +
                                  v * (hCoeffs->rho43v5 + (hCoeffs->rho43v6 + hCoeffs->rho43v6l * eulerlogxabs) * v)));
            rhoLow = 1.;
            A0lm = g2 * (hCoeffs->h43T0ff21 * ff21 + hCoeffs->h43T0ff12 * I * ff12 + hCoeffs->h43T0ff03 * ff03 +
                         hCoeffs->h43T0ff01 * ff01);

            auxflm = x0 * g2 *
                     (hCoeffs->h43T1ff20 * ff20 + hCoeffs->h43T1ff11 * I * ff11 + hCoeffs->h43T1ff02 * ff02 +
                      hCoeffs->h43T1ff00 * ff00);
            break;
        case 2:
            Use_V1 = 1;
            rhoHigh = v2 * v *
                      (hCoeffs->rho42v3 +
                       v * (hCoeffs->rho42v4 +
                            v * (hCoeffs->rho42v5 + (hCoeffs->rho42v6 + hCoeffs->rho42v6l * eulerlogxabs) * v)));
            rhoLow = 1. + v2 * hCoeffs->rho42v2;
            A0lm = hCoeffs->h42T0ff40 * ff40 + hCoeffs->h42T0ff31 * I * ff31 + hCoeffs->h42T0ff20 * ff20 +
                   hCoeffs->h42T0ff13 * I * ff13 + hCoeffs->h42T0ff11 * I * ff11 + hCoeffs->h42T0ff04 * ff04 +
                   hCoeffs->h42T0ff02 * ff02 + hCoeffs->h42T0ff00 * ff00;
            A2lm = hCoeffs->h42T2ff60 * ff60 + hCoeffs->h42T2ff51 * I * ff51 + hCoeffs->h42T2ff42 * ff42 +
                   hCoeffs->h42T2ff40 * ff40 + hCoeffs->h42T2ff33 * I * ff33 + hCoeffs->h42T2ff31 * I * ff31 +
                   hCoeffs->h42T2ff24 * ff24 + hCoeffs->h42T2ff22 * ff22 + hCoeffs->h42T2ff20 * ff20 +
                   hCoeffs->h42T2ff15 * I * ff15 + hCoeffs->h42T2ff13 * I * ff13 + hCoeffs->h42T2ff11 * I * ff11 +
                   hCoeffs->h42T2ff06 * ff06 + hCoeffs->h42T2ff04 * ff04 + hCoeffs->h42T2ff02 * ff02 +
                   hCoeffs->h42T2ff00 * ff00;
            break;
        case 1:
            Use_V1 = 1;
            rhoHigh = v2 * (hCoeffs->rho41v2 +
                            v2 * (hCoeffs->rho41v4 +
                                  v * (hCoeffs->rho41v5 + (hCoeffs->rho41v6 + hCoeffs->rho41v6l * eulerlogxabs) * v)));
            rhoLow = 1.;
            A0lm = g2 * (hCoeffs->h41T0ff21 * ff21 + hCoeffs->h41T0ff12 * I * ff12 + hCoeffs->h41T0ff03 * ff03 +
                         hCoeffs->h41T0ff01 * ff01);

            auxflm = x0 * g2 *
                     (hCoeffs->h41T1ff20 * ff20 + hCoeffs->h41T1ff11 * I * ff11 + hCoeffs->h41T1ff02 * ff02 +
                      hCoeffs->h41T1ff00 * ff00);
            break;
        default:
            PRINT_LOG_INFO(LOG_CRITICAL,
                           "Currently mode (%d, %d) is unsupported in new "
                           "factorized waveform.",
                           l, m);
            return CEV_FAILURE;
        }
        break;

    case 5:
        Use_V1 = 1;
        switch (ABS(m))
        {
        case 5:
            // RC: This terms are in Eq.A9 in
            // https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028
            // [arXiv:1803.10701]
            rhoLow =
                1. + v2 * (hCoeffs->rho55v2 +
                           v * (hCoeffs->rho55v3 +
                                v * (hCoeffs->rho55v4 +
                                     v * (hCoeffs->rho55v5 +
                                          v * (hCoeffs->rho55v6 + hCoeffs->rho55v6l * eulerlogxabs +
                                               v2 * (hCoeffs->rho55v8 + hCoeffs->rho55v8l * eulerlogxabs +
                                                     v2 * (hCoeffs->rho55v10 + hCoeffs->rho55v10l * eulerlogxabs)))))));
            // RC: This terms are in Eq.A12 in
            // https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028
            // [arXiv:1803.10701]
            auxflm = v2 * v * (hCoeffs->f55v3 + v * (hCoeffs->f55v4));
            break;
        case 4:
            rhoLow = 1. + v2 * (hCoeffs->rho54v2 + v * (hCoeffs->rho54v3 + hCoeffs->rho54v4 * v));
            break;
        case 3:
            rhoLow =
                1. + v2 * (hCoeffs->rho53v2 + v * (hCoeffs->rho53v3 + v * (hCoeffs->rho53v4 + hCoeffs->rho53v5 * v)));
            break;
        case 2:
            rhoLow = 1. + v2 * (hCoeffs->rho52v2 + v * (hCoeffs->rho52v3 + hCoeffs->rho52v4 * v));
            break;
        case 1:
            rhoLow =
                1. + v2 * (hCoeffs->rho51v2 + v * (hCoeffs->rho51v3 + v * (hCoeffs->rho51v4 + hCoeffs->rho51v5 * v)));
            break;
        default:
            PRINT_LOG_INFO(LOG_CRITICAL,
                           "Currently mode (%d, %d) is unsupported in new "
                           "factorized waveform.",
                           l, m);
            return CEV_FAILURE;
        }
        break;
    case 6:
        Use_V1 = 1;
        switch (ABS(m))
        {
        case 6:
            rhoLow = 1. + v2 * (hCoeffs->rho66v2 + v * (hCoeffs->rho66v3 + hCoeffs->rho66v4 * v));
            break;
        case 5:
            rhoLow = 1. + v2 * (hCoeffs->rho65v2 + hCoeffs->rho65v3 * v);
            break;
        case 4:
            rhoLow = 1. + v2 * (hCoeffs->rho64v2 + v * (hCoeffs->rho64v3 + hCoeffs->rho64v4 * v));
            break;
        case 3:
            rhoLow = 1. + v2 * (hCoeffs->rho63v2 + hCoeffs->rho63v3 * v);
            break;
        case 2:
            rhoLow = 1. + v2 * (hCoeffs->rho62v2 + v * (hCoeffs->rho62v3 + hCoeffs->rho62v4 * v));
            break;
        case 1:
            rhoLow = 1. + v2 * (hCoeffs->rho61v2 + hCoeffs->rho61v3 * v);
            break;
        default:
            PRINT_LOG_INFO(LOG_CRITICAL,
                           "Currently mode (%d, %d) is unsupported in new "
                           "factorized waveform.",
                           l, m);
            return CEV_FAILURE;
        }
        break;
    case 7:
        Use_V1 = 1;
        switch (ABS(m))
        {
        case 7:
            rhoLow = 1. + v2 * (hCoeffs->rho77v2 + hCoeffs->rho77v3 * v);
            break;
        case 6:
            rhoLow = 1. + hCoeffs->rho76v2 * v2;
            break;
        case 5:
            rhoLow = 1. + v2 * (hCoeffs->rho75v2 + hCoeffs->rho75v3 * v);
            break;
        case 4:
            rhoLow = 1. + hCoeffs->rho74v2 * v2;
            break;
        case 3:
            rhoLow = 1. + v2 * (hCoeffs->rho73v2 + hCoeffs->rho73v3 * v);
            break;
        case 2:
            rhoLow = 1. + hCoeffs->rho72v2 * v2;
            break;
        case 1:
            rhoLow = 1. + v2 * (hCoeffs->rho71v2 + hCoeffs->rho71v3 * v);
            break;
        default:
            PRINT_LOG_INFO(LOG_CRITICAL,
                           "Currently mode (%d, %d) is unsupported in new "
                           "factorized waveform.",
                           l, m);
            return CEV_FAILURE;
        }
        break;
    case 8:
        Use_V1 = 1;
        switch (ABS(m))
        {
        case 8:
            rhoLow = 1. + hCoeffs->rho88v2 * v2;
            break;
        case 7:
            rhoLow = 1. + hCoeffs->rho87v2 * v2;
            break;
        case 6:
            rhoLow = 1. + hCoeffs->rho86v2 * v2;
            break;
        case 5:
            rhoLow = 1. + hCoeffs->rho85v2 * v2;
            break;
        case 4:
            rhoLow = 1. + hCoeffs->rho84v2 * v2;
            break;
        case 3:
            rhoLow = 1. + hCoeffs->rho83v2 * v2;
            break;
        case 2:
            rhoLow = 1. + hCoeffs->rho82v2 * v2;
            break;
        case 1:
            rhoLow = 1. + hCoeffs->rho81v2 * v2;
            break;
        default:
            PRINT_LOG_INFO(LOG_CRITICAL,
                           "Currently mode (%d, %d) is unsupported in new "
                           "factorized waveform.",
                           l, m);
            return CEV_FAILURE;
        }
        break;
    default:
        PRINT_LOG_INFO(LOG_CRITICAL, "Currently mode (%d, %d) is unsupported in new factorized waveform.", l, m);
        return CEV_FAILURE;
    }
    if (Use_V1)
    {
        rholmPwrl = 1.0;
        rholm = rhoLow + rhoHigh;
        i = l;
        while (i--)
            rholmPwrl *= rholm;
    }
    else
    {
        COMPLEX16 T0, T1, T2, T3, T4;
        REAL8 LL = (REAL8)l;
        T0 = 1.;
        T1 = -A1lm / A0lm;
        T2 = cpow(A1lm / A0lm, 2.) - A2lm / A0lm;
        T3 = -cpow(A1lm / A0lm, 3) + 2. * A1lm * A2lm / A0lm / A0lm - A3lm / A0lm;
        T4 = cpow(A1lm / A0lm, 4.) - 3. * A1lm * A1lm * A2lm / cpow(A0lm, 3.) + cpow(A2lm / A0lm, 2) +
             2. * A1lm * A3lm / A0lm - A4lm / A0lm;
        rho2PN = A0lm / (T0 + T1 * x0 + T2 * x0_2 + T3 * x0_3 + T4 * x0_4);
        rholm = rhoLow + rhoHigh;
        rholmPwrl = 1.0;
        REAL8 rholmPwrlSub = 1.0;
        i = l;
        while (i--)
        {
            rholmPwrl *= rholm;
            rholmPwrlSub *= rhoLow;
        }
        // rholmPwrl = rho2PN;
        rholmPwrl = rho2PN + (rholmPwrl - rholmPwrlSub) * pow(vPhi / x0, l + (l + m) % 2);
        // rholmPwrl = rholmPwrl * pow(vPhi/x0, l + (l + m) % 2);
        // if (isnan(rho2PN))
        // print_debug("(%d,%d), h21T2ff21 = %g, h21T2ff12 = %g, h21T2ff10 = %g,
        // h21T2ff03 = %g, h21T2ff01 = %g + i %g\n", l, m,
        //         hCoeffs->h21T2ff21, hCoeffs->h21T2ff12,
        //         hCoeffs->h21T2ff10, hCoeffs->h21T2ff03,
        //         creal(hCoeffs->h21T2ff01), cimag(hCoeffs->h21T2ff01));
    }
    if (eta == 0.25 && m % 2)
    {
        rholmPwrl = auxflm;
    }
    else
    {
        rholmPwrl += auxflm;
    }
    *facAmpRes = rholmPwrl;
    return CEV_SUCCESS;
}

/**
 * This function calculates hlm mode factorized-resummed waveform
 * for given dynamical variables.
 * Eq. 17 and the entire Appendix of the paper.
 */
INT XLALSimIMRSpinEOBGetSpinFactorizedWaveformV2(COMPLEX16 *ret,
                                                 /**< OUTPUT, hlm waveforms */
                                                 REAL8Vector *values,
                                                 /**< dyanmical variables */
                                                 const REAL8 v,
                                                 /**< velocity */
                                                 const REAL8 dr,
                                                 /**< radial velocity */
                                                 const REAL8 ncrv,
                                                 /**< nVec Cross vVec */
                                                 const REAL8 Hreal,
                                                 /**< real Hamiltonian */
                                                 const INT4 l,
                                                 /**< l mode index */
                                                 const INT4 m,
                                                 /**< m mode index */
                                                 SpinEOBParams *params,
                                                 /**< Spin EOB parameters */
                                                 INT ret_type
                                                 /**< ret_type: 0:full, 1:amp, 2:amp_res */
)
{
    /* Status of function calls */
    INT4 status;
    INT4 i;
    INT4 use_hm = 0;

    REAL8 eta;
    REAL8 r, pp, Omega, v2, vh, vh3, k, hathatk, eulerlogxabs; // pr
    REAL8 Slm, deltalm, rholm;
    COMPLEX16 auxflm = 0.0;
    COMPLEX16 Tlm, rholmPwrl;
    COMPLEX16 fTe = 0.0;
    REAL8 hathatksq4, hathatk4pi, Tlmprefac, Tlmprodfac;
    COMPLEX16 hNewton;
    gsl_sf_result lnr1, arg1, z2;

    /* Non-Keplerian velocity */
    REAL8 vPhi, vPhi2;

    /* Pre-computed coefficients */
    FacWaveformCoeffs *hCoeffs = params->hCoeffs;
    use_hm = params->use_hm;
    eta = params->eta;

    r = values->data[0];
    // pr    = values->data[2];
    pp = values->data[3];

    v2 = v * v;
    Omega = v2 * v;
    vh3 = Hreal * Omega;
    vh = cbrt(vh3);
    eulerlogxabs = CST_GAMMA + log(2.0 * (REAL8)m * v);

    /* Calculate the non-Keplerian velocity */
    if (params->alignedSpins)
    {
        // YP: !!!!! SEOBNRv3devel temporary change !!!!!
        vPhi = XLALSimIMRSpinAlignedEOBNonKeplerCoeff(values->data, params);

        if (IS_REAL8_FAIL_NAN(vPhi))
        {
            return CEV_FAILURE;
        }

        vPhi = r * cbrt(vPhi);
        vPhi *= Omega;
        vPhi2 = vPhi * vPhi;
    }
    else
    {
        vPhi = v;
        vPhi2 = v2;
    }
    /* Calculate the newtonian multipole, 1st term in Eq. 17, given by Eq. A1 */
    // YP: !!!!! SEOBNRv3devel temporary change !!!!!
    if (ret_type <= 0)
        status = XLALSimIMRSpinEOBCalculateNewtonianMultipole(&hNewton, vPhi2, r, values->data[1], (UINT)l, m, params);
    else
        status =
            XLALSimIMRSpinEOBFluxCalculateNewtonianMultipole(&hNewton, vPhi2, r, values->data[1], (UINT)l, m, params);
    // YP: !!!!! SEOBNRv3devel temporary change !!!!!
    if (status != CEV_SUCCESS)
    {
        return CEV_FAILURE;
    }
    /* Calculate the source term, 2nd term in Eq. 17, given by Eq. A5 */
    if (((l + m) % 2) == 0)
    {
        Slm = (Hreal * Hreal - 1.) / (2. * eta) + 1.;
    }
    else
    {
        Slm = v * pp;
    }
    // printf( "Hreal = %e, Slm = %e, eta = %e\n", Hreal, Slm, eta );

    /* Calculate the Tail term, 3rd term in Eq. 17, given by Eq. A6 */
    k = m * Omega;
    hathatk = Hreal * k;
    if (ret_type <= 0)
    {
        status = gsl_sf_lngamma_complex_e(l + 1.0, -2.0 * hathatk, &lnr1, &arg1);
        if (status != GSL_SUCCESS)
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "Error in GSL function: %s", gsl_strerror(status));
            return CEV_FAILURE;
        }
        status = gsl_sf_fact_e(l, &z2);
        if (status != GSL_SUCCESS)
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "Error in GSL function: %s", gsl_strerror(status));
            return CEV_FAILURE;
        }
        Tlm = cexp((lnr1.val + CST_PI * hathatk) + I * (arg1.val + 2.0 * hathatk * log(4.0 * k / sqrt(CST_E))));
        Tlm /= z2.val;
        switch (l)
        {
        case 2:
            switch (abs(m))
            {
            case 2:
                deltalm =
                    vh3 * (hCoeffs->delta22vh3 + vh3 * (hCoeffs->delta22vh6 + vh * vh * (hCoeffs->delta22vh9 * vh))) +
                    hCoeffs->delta22v5 * v * v2 * v2 + hCoeffs->delta22v6 * v2 * v2 * v2 +
                    hCoeffs->delta22v8 * v2 * v2 * v2 * v2;
                break;
            case 1:
                deltalm =
                    vh3 * (hCoeffs->delta21vh3 +
                           vh3 * (hCoeffs->delta21vh6 + vh * (hCoeffs->delta21vh7 + vh * vh * (hCoeffs->delta21vh9)))) +
                    hCoeffs->delta21v5 * v * v2 * v2 + hCoeffs->delta21v7 * v2 * v2 * v2 * v;
                break;
            default:
                return CEV_FAILURE;
                break;
            }
            break;
        case 3:
            switch (m)
            {
            case 3:
                deltalm = vh3 * (hCoeffs->delta33vh3 + vh3 * (hCoeffs->delta33vh6 + hCoeffs->delta33vh9 * vh3)) +
                          hCoeffs->delta33v5 * v * v2 * v2;
                break;
            case 2:
                deltalm =
                    vh3 * (hCoeffs->delta32vh3 +
                           vh * (hCoeffs->delta32vh4 + vh * vh * (hCoeffs->delta32vh6 + hCoeffs->delta32vh9 * vh3)));
                break;
            case 1:
                deltalm =
                    vh3 * (hCoeffs->delta31vh3 +
                           vh3 * (hCoeffs->delta31vh6 + vh * (hCoeffs->delta31vh7 + vh * vh * hCoeffs->delta31vh9))) +
                    hCoeffs->delta31v5 * v * v2 * v2;
                break;
            default:
                return CEV_FAILURE;
                break;
            }
            break;
        case 4:
            switch (m)
            {
            case 4:
                if (use_hm)
                {
                    // RC: This terms are in Eq.A15 in
                    // https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028
                    deltalm = vh3 * (hCoeffs->delta44vh3 + vh3 * (hCoeffs->delta44vh6 + vh3 * hCoeffs->delta44vh9)) +
                              hCoeffs->delta44v5 * v2 * v2 * v;
                }
                else
                {
                    deltalm =
                        vh3 * (hCoeffs->delta44vh3 + vh3 * hCoeffs->delta44vh6) + hCoeffs->delta44v5 * v2 * v2 * v;
                }
                break;
            case 3:
                deltalm = vh3 * (hCoeffs->delta43vh3 + vh * (hCoeffs->delta43vh4 + vh * vh * hCoeffs->delta43vh6));
                break;
            case 2:
                deltalm = vh3 * (hCoeffs->delta42vh3 + hCoeffs->delta42vh6 * vh3);
                break;
            case 1:
                deltalm = vh3 * (hCoeffs->delta41vh3 + vh * (hCoeffs->delta41vh4 + vh * vh * hCoeffs->delta41vh6));
                break;
            default:
                return CEV_FAILURE;
                break;
            }
            break;
        case 5:
            switch (m)
            {
            case 5:
                if (use_hm)
                {
                    // RC: This terms are in Eq.A16 in
                    // https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028
                    deltalm = vh3 * (hCoeffs->delta55vh3 + vh3 * (hCoeffs->delta55vh6 + vh3 * (hCoeffs->delta55vh9))) +
                              hCoeffs->delta55v5 * v2 * v2 * v;
                }
                else
                {
                    deltalm = hCoeffs->delta55vh3 * vh3 + hCoeffs->delta55v5 * v2 * v2 * v;
                }
                break;
            case 4:
                deltalm = vh3 * (hCoeffs->delta54vh3 + hCoeffs->delta54vh4 * vh);
                break;
            case 3:
                deltalm = hCoeffs->delta53vh3 * vh3;
                break;
            case 2:
                deltalm = vh3 * (hCoeffs->delta52vh3 + hCoeffs->delta52vh4 * vh);
                break;
            case 1:
                deltalm = hCoeffs->delta51vh3 * vh3;
                break;
            default:
                return CEV_FAILURE;
                break;
            }
            break;
        case 6:
            switch (m)
            {
            case 6:
                deltalm = hCoeffs->delta66vh3 * vh3;
                break;
            case 5:
                deltalm = hCoeffs->delta65vh3 * vh3;
                break;
            case 4:
                deltalm = hCoeffs->delta64vh3 * vh3;
                break;
            case 3:
                deltalm = hCoeffs->delta63vh3 * vh3;
                break;
            case 2:
                deltalm = hCoeffs->delta62vh3 * vh3;
                break;
            case 1:
                deltalm = hCoeffs->delta61vh3 * vh3;
                break;
            default:
                return CEV_FAILURE;
                break;
            }
            break;
        case 7:
            switch (m)
            {
            case 7:
                deltalm = hCoeffs->delta77vh3 * vh3;
                break;
            case 5:
                deltalm = hCoeffs->delta75vh3 * vh3;
                break;
            case 3:
                deltalm = hCoeffs->delta73vh3 * vh3;
                break;
            case 1:
                deltalm = hCoeffs->delta71vh3 * vh3;
                break;
            case 6:
            case 4:
            case 2:
                deltalm = 0.0;
                break;
            default:
                return CEV_FAILURE;
                break;
            }
            break;
        case 8:
            deltalm = 0.0;
            break;
        default:
            return CEV_FAILURE;
            break;
        }
    }
    else
    {
        hathatksq4 = 4. * hathatk * hathatk;
        hathatk4pi = 4. * CST_PI * hathatk;
        status = gsl_sf_fact_e(l, &z2);
        if (status != GSL_SUCCESS)
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "Error in GSL function: %s", gsl_strerror(status));
            return CEV_FAILURE;
        }
        /* Calculating the prefactor of Tlm, outside the multiple product */
        Tlmprefac = sqrt(hathatk4pi / (1. - exp(-hathatk4pi))) / z2.val;

        /* Calculating the multiple product factor */
        for (Tlmprodfac = 1., i = 1; i <= l; i++)
        {
            Tlmprodfac *= (hathatksq4 + (REAL8)i * i);
        }

        Tlm = Tlmprefac * sqrt(Tlmprodfac);
        deltalm = 0.0;
    }
    if (CalculateFacWaveformAmpResV2(hCoeffs, eta, l, m, v, vPhi, vh3, eulerlogxabs, 1 / sqrt(r), dr, ncrv,
                                     &rholmPwrl) != CEV_SUCCESS)
        return CEV_FAILURE;
    if (ret_type > 1)
    {
        *ret = rholmPwrl;
        return CEV_SUCCESS;
    }
    /* Put all factors in Eq. 17 together */
    *ret = (Tlm * cexp(I * deltalm) + fTe) * Slm * (rholmPwrl);
    *ret *= hNewton;
    /*if (r > 8.5)
    {
        printf("YP::FullWave: Reh = %.16e, Imh = %.16e, hAmp = %.16e, hPhi =
    %.16e\n",creal(*hlm),cimag(*hlm),cabs(*hlm),carg(*hlm));
    } */
    return CEV_SUCCESS;
}

REAL8
XLALInspiralSpinFactorizedFluxV2(REAL8Vector *values, /**< dynamical variables */
                                 EOBNonQCCoeffs *nqcCoeffs,
                                 /**< pre-computed NQC coefficients */
                                 const REAL8 omega, /**< orbital frequency */
                                 const REAL8 dr,    /**< radial velocity */
                                 const REAL8 ncrv,  /**< nVec cross vVec */
                                 SpinEOBParams *ak, /**< physical parameters */
                                 const REAL8 H,     /**< real Hamiltonian */
                                 const INT lMax     /**< upper limit of the summation over l */
)
{
    REAL8 flux = 0.0;
    REAL8 v;
    REAL8 omegaSq;
    COMPLEX16 hLM, hNQC;
    INT l, m;

    /* Omegs is the derivative of phi */
    omegaSq = omega * omega;

    v = cbrt(omega);

    //  printf( "v = %.16e\n", v );
    for (l = 2; l <= lMax; l++)
    {
        for (m = 1; m <= l; m++)
        {
            if (XLALSimIMRSpinEOBGetSpinFactorizedWaveformV2(&hLM, values, v, dr, ncrv, H, l, m, ak, 1) != CEV_SUCCESS)
            {
                return REAL8_FAIL_NAN;
            }

            /* For the 2,2 mode, we apply NQC correction to the flux */
            // if (l == 2 && m == 2)
            // {
            //     XLALSimIMREOBNonQCCorrection (&hNQC, values, omega,
            //     nqcCoeffs);
            //     /* Eq. 16 */
            //     hLM *= hNQC;
            // }
            // printf( "l = %d, m = %d, mag(hLM) = %.17e, omega = %.16e\n", l, m,
            // sqrt(creal(hLM)*creal(hLM)+cimag(hLM)*cimag(hLM)), omega );
            /* Eq. 13 */
            flux += (REAL8)(m * m) * omegaSq * (creal(hLM) * creal(hLM) + cimag(hLM) * cimag(hLM));
        }
    }
    return flux * CST_1_PI / 8.0;
}

/*--------------------------------------------------------------*/
/*                                                              */
/*                                                              */
/*                                                              */
/*                            PREC                              */
/*                                                              */
/*                                                              */
/*                                                              */
/*--------------------------------------------------------------*/
INT XLALSimIMRSpinEOBFluxGetPrecSpinFactorizedWaveformV2(
    COMPLEX16 *hlmTab,                                   /**< OUTPUT, hlm waveforms */
    REAL8Vector *values,                                 /**< dyanmical variables: \f$(r,\phi,p_r,p_\phi)\f$ */
    REAL8Vector *cartvalues,                             /**< dyanmical variables */
    const REAL8 v,                                       /**< velocity */
    const REAL8 dr, const REAL8 ncrv, const REAL8 Hreal, /**< real Hamiltonian */
    const INT4 lMax,                                     /**< maximum l mode to compute, compute 0 < m <= lMax */
    SpinEOBParams *params                                /**< Spin EOB parameters */
)
{
    // int		debugPK = 0;
    const REAL8 vPhiKepler = params->alignedSpins ? XLALSimIMRSpinAlignedEOBNonKeplerCoeff(values->data, params)
                                                  : XLALSimIMRSpinPrecEOBNonKeplerCoeff(cartvalues->data, params);
    if (IS_REAL8_FAIL_NAN(vPhiKepler))
    {
        return CEV_FAILURE;
    }
    INT l, m;
    for (l = 2; l <= lMax; l++)
    {
        for (m = 1; m <= l; m++)
        {
            COMPLEX16 *hlm = &hlmTab[l * (lMax + 1) + m];
            /* Status of function calls */
            INT status;
            INT i;

            REAL8 eta;
            REAL8 r, pp, Omega, v2, vh3, /* vh, vh3, */ k, hathatk, eulerlogxabs;
            // pr
            REAL8 rcrossp_x, rcrossp_y, rcrossp_z;
            REAL8 Slm, rholm;
            COMPLEX16 rholmPwrl;
            REAL8 auxflm = 0.0;
            REAL8 hathatksq4, hathatk4pi, Tlmprefac, Tlmprodfac;
            REAL8 Tlm;
            COMPLEX16 hNewton;
            gsl_sf_result z2;

            /* Non-Keplerian velocity */
            REAL8 vPhi, vPhi2;

            /* Pre-computed coefficients */
            FacWaveformCoeffs *hCoeffs = params->hCoeffs;

            eta = params->eta;

            /*
             * else if ( eta == 0.25 && m % 2 ) { // If m is odd and dM = 0, hLM
             * will be zero memset( hlm, 0, sizeof( COMPLEX16 ) ); return
             * XLAL_SUCCESS; }
             */

            // r = sqrt(values->data[0] * values->data[0] + values->data[1] *
            // values->data[1] + values->data[2] * values->data[2]); pr =
            // values->data[2];
            r = values->data[0];
            pp = values->data[3];

            rcrossp_x = cartvalues->data[1] * cartvalues->data[5] - cartvalues->data[2] * cartvalues->data[4];
            rcrossp_y = cartvalues->data[2] * cartvalues->data[3] - cartvalues->data[0] * cartvalues->data[5];
            rcrossp_z = cartvalues->data[0] * cartvalues->data[4] - cartvalues->data[1] * cartvalues->data[3];

            // pp = sqrt(rcrossp_x * rcrossp_x + rcrossp_y * rcrossp_y +
            // rcrossp_z * rcrossp_z);

            v2 = v * v;
            Omega = v2 * v;
            vh3 = Hreal * Omega;
            // vh = cbrt(vh3);
            eulerlogxabs = CST_GAMMA + log(2.0 * (REAL8)m * v);

            /* Calculate the non-Keplerian velocity */
            vPhi = vPhiKepler;

            vPhi = r * cbrt(vPhi);

            // if (debugPK)
            // 	XLAL_PRINT_INFO("In
            // XLALSimIMRSpinEOBFluxCalculateNewtonianMultipole, getting rW =
            // %.12e\n", vPhi);
            vPhi *= Omega;
            vPhi2 = vPhi * vPhi;

            /*
             * Calculate the newtonian multipole, 1st term in Eq. 17, given by
             * Eq. A1
             */
            // debugPK
            //  if (debugPK) {
            //  	XLAL_PRINT_INFO("\nValues inside
            //  XLALSimIMRSpinEOBFluxGetPrecSpinFactorizedWaveform:\n"); 	for (i
            //  = 0; i
            //  < 14; i++) 		XLAL_PRINT_INFO("values[%d] = %.12e\n",
            //  i, values->data[i]);

            // 	XLAL_PRINT_INFO("Calculating hNewton, with v = %.12e, vPhi =
            // %.12e, r = %.12e, Phi = %.12e, l = %d, m = %d\n", 	       v, vPhi,
            // r, values->data[1], (UINT4) l, (UINT4) m);
            // }
            status = XLALSimIMRSpinEOBFluxCalculateNewtonianMultipole(&hNewton, vPhi2, r, values->data[1], (UINT)l, m,
                                                                      params);
            if (status == CEV_FAILURE)
            {
                return CEV_FAILURE;
            }
            /* Calculate the source term, 2nd term in Eq. 17, given by Eq. A5 */
            if (((l + m) % 2) == 0)
            {
                Slm = (Hreal * Hreal - 1.) / (2. * eta) + 1.;
            }
            else
            {
                Slm = v * pp;
                // Slm = v * sqrt(rcrossp_x * rcrossp_x + rcrossp_y * rcrossp_y +
                // rcrossp_z * rcrossp_z);
            }
            // if (debugPK)
            // 	XLAL_PRINT_INFO("In
            // XLALSimIMRSpinEOBFluxGetSpinFactorizedWaveform: Hreal = %e, Slm =
            // %e, eta = %e\n", Hreal, Slm, eta);

            /*
             * Calculate the absolute value of the Tail term, 3rd term in Eq. 17,
             * given by Eq. A6, and Eq. (42) of
             * http://arxiv.org/pdf/1212.4357.pdf
             */
            k = m * Omega;
            hathatk = Hreal * k;
            hathatksq4 = 4. * hathatk * hathatk;
            hathatk4pi = 4. * CST_PI * hathatk;
            /*
             * gsl_sf_result lnr1, arg1; XLAL_CALLGSL( status =
             * gsl_sf_lngamma_complex_e( l+1.0, -2.0*hathatk, &lnr1, &arg1 ) );
             * if (status != GSL_SUCCESS) { XLALPrintError("XLAL Error - %s:
             * Error in GSL function\n", __func__ ); XLAL_ERROR( XLAL_EFUNC ); }
             */
            status = gsl_sf_fact_e(l, &z2);
            if (status != GSL_SUCCESS)
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "Error in GSL function");
                return CEV_FAILURE;
            }
            /*
             * COMPLEX16 Tlmold; Tlmold = cexp( ( lnr1.val + LAL_PI * hathatk ) +
             * I * ( arg1.val + 2.0 * hathatk * log(4.0*k/sqrt(LAL_E)) ) );
             * Tlmold /= z2.val;
             */
            /* Calculating the prefactor of Tlm, outside the multiple product */
            Tlmprefac = sqrt(hathatk4pi / (1. - exp(-hathatk4pi))) / z2.val;

            /* Calculating the multiple product factor */
            for (Tlmprodfac = 1., i = 1; i <= l; i++)
            {
                Tlmprodfac *= (hathatksq4 + (REAL8)i * i);
            }

            Tlm = Tlmprefac * sqrt(Tlmprodfac);

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
             * Actual values of the coefficients are defined in the next function
             * of this file
             */
#if 0
			switch (l) {
			case 2:
				switch (abs(m)) {
				case 2:
					rholm = 1. + v2 * (hCoeffs->rho22v2 + v * (hCoeffs->rho22v3
								     + v * (hCoeffs->rho22v4
					     + v * (hCoeffs->rho22v5 + v * (hCoeffs->rho22v6
									    + hCoeffs->rho22v6l * eulerlogxabs + v * (hCoeffs->rho22v7
														      + v * (hCoeffs->rho22v8 + hCoeffs->rho22v8l * eulerlogxabs
															     + (hCoeffs->rho22v10 + hCoeffs->rho22v10l * eulerlogxabs) * v2)))))));
					//FIXME
						// if (debugPK){
                        //     XLAL_PRINT_INFO("PK:: rho22v2 = %.12e, rho22v3 = %.12e, rho22v4 = %.12e,\n rho22v5 = %.16e, rho22v6 = %.16e, rho22v6LOG = %.16e, \n rho22v7 = %.12e, rho22v8 = %.16e, rho22v8LOG = %.16e, \n rho22v10 = %.16e, rho22v10LOG = %.16e\n, rho22v6 = %.12e, rho22v8 = %.12e, rho22v10 = %.12e\n",
						//        hCoeffs->rho22v2, hCoeffs->rho22v3, hCoeffs->rho22v4,
						//        hCoeffs->rho22v5, hCoeffs->rho22v6, hCoeffs->rho22v6l,
						//        hCoeffs->rho22v7, hCoeffs->rho22v8, hCoeffs->rho22v8l,
						//        hCoeffs->rho22v10, hCoeffs->rho22v10l,
						//        hCoeffs->rho22v6 + hCoeffs->rho22v6l * eulerlogxabs,
						//        hCoeffs->rho22v8 + hCoeffs->rho22v8l * eulerlogxabs,
						//        hCoeffs->rho22v10 + hCoeffs->rho22v10l * eulerlogxabs);}
					break;
				case 1:
					{
						rholm = 1. + v * (hCoeffs->rho21v1
								  + v * (hCoeffs->rho21v2 + v * (hCoeffs->rho21v3 + v * (hCoeffs->rho21v4
															 + v * (hCoeffs->rho21v5 + v * (hCoeffs->rho21v6 + hCoeffs->rho21v6l * eulerlogxabs
																			+ v * (hCoeffs->rho21v7 + hCoeffs->rho21v7l * eulerlogxabs
																			       + v * (hCoeffs->rho21v8 + hCoeffs->rho21v8l * eulerlogxabs
																				      + (hCoeffs->rho21v10 + hCoeffs->rho21v10l * eulerlogxabs) * v2))))))));
						auxflm = v * hCoeffs->f21v1 + v2 * v * hCoeffs->f21v3;
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
					rholm = 1. + v2 * (hCoeffs->rho33v2 + v * (hCoeffs->rho33v3 + v * (hCoeffs->rho33v4
													   + v * (hCoeffs->rho33v5 + v * (hCoeffs->rho33v6 + hCoeffs->rho33v6l * eulerlogxabs
																	  + v * (hCoeffs->rho33v7 + (hCoeffs->rho33v8 + hCoeffs->rho33v8l * eulerlogxabs) * v))))));
					auxflm = v * v2 * hCoeffs->f33v3;
					break;
				case 2:
					rholm = 1. + v * (hCoeffs->rho32v
							  + v * (hCoeffs->rho32v2 + v * (hCoeffs->rho32v3 + v * (hCoeffs->rho32v4 + v * (hCoeffs->rho32v5
																	 + v * (hCoeffs->rho32v6 + hCoeffs->rho32v6l * eulerlogxabs
																		+ (hCoeffs->rho32v8 + hCoeffs->rho32v8l * eulerlogxabs) * v2))))));
					break;
				case 1:
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

					rholm = 1. + v2 * (hCoeffs->rho44v2
					     + v * (hCoeffs->rho44v3 + v * (hCoeffs->rho44v4
						 + v * (hCoeffs->rho44v5 + (hCoeffs->rho44v6
						+ hCoeffs->rho44v6l * eulerlogxabs) * v))));
					break;
				case 3:
					rholm = 1. + v * (hCoeffs->rho43v
							  + v * (hCoeffs->rho43v2
					    + v2 * (hCoeffs->rho43v4 + v * (hCoeffs->rho43v5
									    + (hCoeffs->rho43v6 + hCoeffs->rho43v6l * eulerlogxabs) * v))));
					auxflm = v * hCoeffs->f43v;
					break;
				case 2:
					rholm = 1. + v2 * (hCoeffs->rho42v2
							   + v * (hCoeffs->rho42v3 + v * (hCoeffs->rho42v4 + v * (hCoeffs->rho42v5
														  + (hCoeffs->rho42v6 + hCoeffs->rho42v6l * eulerlogxabs) * v))));
					break;
				case 1:
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
					rholm = 1. + v2 * (hCoeffs->rho55v2
					     + v * (hCoeffs->rho55v3 + v * (hCoeffs->rho55v4
					 + v * (hCoeffs->rho55v5 + hCoeffs->rho55v6 * v))));
					break;
				case 4:
					rholm = 1. + v2 * (hCoeffs->rho54v2 + v * (hCoeffs->rho54v3
								   + hCoeffs->rho54v4 * v));
					break;
				case 3:
					rholm = 1. + v2 * (hCoeffs->rho53v2
							   + v * (hCoeffs->rho53v3 + v * (hCoeffs->rho53v4 + hCoeffs->rho53v5 * v)));
					break;
				case 2:
					rholm = 1. + v2 * (hCoeffs->rho52v2 + v * (hCoeffs->rho52v3
								   + hCoeffs->rho52v4 * v));
					break;
				case 1:
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
					rholm = 1. + v2 * (hCoeffs->rho66v2 + v * (hCoeffs->rho66v3
								   + hCoeffs->rho66v4 * v));
					break;
				case 5:
					rholm = 1. + v2 * (hCoeffs->rho65v2 + hCoeffs->rho65v3 * v);
					break;
				case 4:
					rholm = 1. + v2 * (hCoeffs->rho64v2 + v * (hCoeffs->rho64v3
								   + hCoeffs->rho64v4 * v));
					break;
				case 3:
					rholm = 1. + v2 * (hCoeffs->rho63v2 + hCoeffs->rho63v3 * v);
					break;
				case 2:
					rholm = 1. + v2 * (hCoeffs->rho62v2 + v * (hCoeffs->rho62v3
								   + hCoeffs->rho62v4 * v));
					break;
				case 1:
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
					rholm = 1. + v2 * (hCoeffs->rho77v2 + hCoeffs->rho77v3 * v);
					break;
				case 6:
					rholm = 1. + hCoeffs->rho76v2 * v2;
					break;
				case 5:
					rholm = 1. + v2 * (hCoeffs->rho75v2 + hCoeffs->rho75v3 * v);
					break;
				case 4:
					rholm = 1. + hCoeffs->rho74v2 * v2;
					break;
				case 3:
					rholm = 1. + v2 * (hCoeffs->rho73v2 + hCoeffs->rho73v3 * v);
					break;
				case 2:
					rholm = 1. + hCoeffs->rho72v2 * v2;
					break;
				case 1:
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
					rholm = 1. + hCoeffs->rho88v2 * v2;
					break;
				case 7:
					rholm = 1. + hCoeffs->rho87v2 * v2;
					break;
				case 6:
					rholm = 1. + hCoeffs->rho86v2 * v2;
					break;
				case 5:
					rholm = 1. + hCoeffs->rho85v2 * v2;
					break;
				case 4:
					rholm = 1. + hCoeffs->rho84v2 * v2;
					break;
				case 3:
					rholm = 1. + hCoeffs->rho83v2 * v2;
					break;
				case 2:
					rholm = 1. + hCoeffs->rho82v2 * v2;
					break;
				case 1:
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
				// if (debugPK)
				// XLAL_PRINT_INFO("rho_%d_%d = %.12e \n", l, m, rholm);
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
				rholmPwrl = auxflm;
			} else {
				rholmPwrl += auxflm;
			}
#endif
            if (CalculateFacWaveformAmpResV2(hCoeffs, eta, l, m, v, vPhi, vh3, eulerlogxabs, 1 / sqrt(r), dr, ncrv,
                                             &rholmPwrl) != CEV_SUCCESS)
                return CEV_FAILURE;
            // if (r > 0.0 && debugPK) {
            // 	XLAL_PRINT_INFO("YP::dynamics variables in waveform: %i, %i, r
            // =
            // %.12e, v = %.12e\n", l, m, r, v); 	XLAL_PRINT_INFO("rholm^l =
            // %.16e, Tlm = %.16e + i %.16e, \nSlm = %.16e, hNewton = %.16e + i
            // %.16e, delta =
            // %.16e\n", rholmPwrl, Tlm, 0.0, Slm, creal(hNewton),
            // cimag(hNewton), 0.0);
            // }
            /* Put all factors in Eq. 17 together */
            *hlm = Tlm * Slm * cabs(rholmPwrl);
            *hlm *= hNewton;
            /*
             * if (r > 8.5) { XLAL_PRINT_INFO("YP::FullWave: %.16e,%.16e,
             * %.16e\n",hlm->re,hlm->im,sqrt(hlm->re*hlm->re+hlm->im*hlm->im)); }
             */
        }
    }
    return CEV_SUCCESS;
}

/**
 * This function calculates hlm mode factorized-resummed waveform
 * for given dynamical variables.
 * Eq. 17 and the entire Appendix of PRD 86, 024011 (2012) + changes
 * described in ￼the section "Factorized waveforms" of
 * https://dcc.ligo.org/T1400476
 */
INT XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveformV2(
    COMPLEX16 *hlm,                       /**< OUTPUT, hlm waveforms */
    REAL8Vector *values,                  /**< dyanmical variables: \f$(r,\phi,p_r,p_\phi)\f$ */
    REAL8Vector *cartvalues,              /**< dyanmical variables */
    const REAL8 v,                        /**< velocity */
    const REAL8 dr,                       /**< radial velocity */
    const REAL8 ncrv,                     /**< angular velocity */
    const REAL8 prDot, const REAL8 Hreal, /**< real Hamiltonian */
    const INT l,                          /**< l mode index */
    const INT m,                          /**< m mode index */
    SpinEOBParams *params                 /**< Spin EOB parameters */
)
{
    int debugPK = 0;
    /* Status of function calls */
    INT status;
    INT i;

    REAL8 eta;
    REAL8 r, pp, Omega, v2, Omegav2, vh, vh3, k, hathatk, eulerlogxabs;
    // pr
    REAL8 rcrossp_x, rcrossp_y, rcrossp_z;
    REAL8 Slm, deltalm;
    COMPLEX16 auxflm = 0.0;
    REAL8 hathatksq4, hathatk4pi, Tlmprefac, Tlmprodfac;
    COMPLEX16 Tlm, rholmPwrl, rholm;
    COMPLEX16 fTe = 0.0;
    COMPLEX16 hNewton;
    gsl_sf_result z2;

    /* Non-Keplerian velocity */
    REAL8 vPhi, vPhi2;

    /* Pre-computed coefficients */

    FacWaveformCoeffs *hCoeffs = params->hCoeffs;

    if (abs(m) > (INT)l)
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

    // pp = sqrt(rcrossp_x * rcrossp_x + rcrossp_y * rcrossp_y + rcrossp_z *
    // rcrossp_z);

    v2 = v * v;
    Omega = v2 * v;
    Omegav2 = Omega * v2;
    vh3 = Hreal * Omega;
    vh = cbrt(vh3);
    eulerlogxabs = CST_GAMMA + log(2.0 * (REAL8)m * v);

    FacWaveformCoeffs oldCoeffs = *(params->hCoeffs); // RC: the function XLALSimIMRSpinPrecEOBNonKeplerCoeff
                                                      // is calculating again the coefficients hCoeffs, for
                                                      // the omega. These are different because
    // for the dynamics we are using different coefficients for the 21 mode. I
    // store the old coefficients and the I recover them ofter vPhi is computed.

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
        //     XLAL_PRINT_INFO("In
        //     XLALSimIMRSpinEOBFluxCalculateNewtonianMultipole, getting rW =
        //     %.12e\n", vPhi);
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
        //     XLAL_PRINT_INFO("In
        //     XLALSimIMRSpinEOBFluxCalculateNewtonianMultipole, getting rW =
        //     %.12e\n",
        //             vPhi);
        vPhi *= Omega;
        vPhi2 = vPhi * vPhi;
    }
    *hCoeffs = oldCoeffs; // RC: Here I recover the old coefficients
    /*
     * Calculate the newtonian multipole, 1st term in Eq. 17, given by
     * Eq. A1
     */
    // if (debugPK) {
    //     XLAL_PRINT_INFO("\nValues inside
    //     XLALSimIMRSpinEOBFluxGetSpinFactorizedWaveform:\n"); for (i = 0; i <
    //     11; i++)
    //         XLAL_PRINT_INFO("values[%d] = %.12e\n", i, cartvalues->data[i]);

    //     XLAL_PRINT_INFO("Calculating hNewton, with v = %.12e, vPhi = %.12e, r
    //     =
    //     %.12e, Phi = %.12e, l = %d, m = %d\n",
    //             v, vPhi, r, values->data[1], (UINT4) l, (UINT4) m);
    // }
    status = XLALSimIMRSpinEOBCalculateNewtonianMultipole(
        &hNewton, vPhi2, r, cartvalues->data[12] + cartvalues->data[13], (UINT)l, m, params);

    if (status != CEV_SUCCESS)
    {
        return CEV_FAILURE;
    }
    /* Calculate the source term, 2nd term in Eq. 17, given by Eq. A5, Hreal is
     * given by Eq.5 and Heff is in Eq.2 */
    if (((l + m) % 2) == 0)
    {
        Slm = (Hreal * Hreal - 1.) / (2. * eta) + 1.;
    }
    else
    {
        Slm = v * pp;
        // Slm = v * sqrt(rcrossp_x * rcrossp_x + rcrossp_y * rcrossp_y +
        // rcrossp_z
        // * rcrossp_z);
    }
    // if (debugPK)
    //     XLAL_PRINT_INFO("In XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform:
    //     Hreal = %e, Slm = %e, eta = %e\n", Hreal, Slm, eta);

    /*
     * Calculate the Tail term, 3rd term in Eq. 17,
     * given by Eq. A6, and Eq. (42) of
     * http://arxiv.org/pdf/1212.4357.pdf (or PRD 87 084035 (2013))
     */
    k = m * Omega;
    hathatk = Hreal * k;

    gsl_sf_result lnr1, arg1;
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
    Tlm = cexp((lnr1.val + CST_PI * hathatk) + I * (arg1.val + 2.0 * hathatk * log(4.0 * k / sqrt(CST_E))));
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
    //     XLAL_PRINT_INFO("Tlm = %e + i%e, |Tlm| = %.16e (should be %.16e)\n",
    //     creal(Tlm), cimag(Tlm), cabs(Tlm), Tlmold);
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
    REAL8 pl = sqrt(GET_MAX(r * prDot + 1. / r, 0.0));
    switch (l)
    {
    case 2:
        switch (abs(m))
        {
        case 2:
            deltalm = vh3 * (hCoeffs->delta22vh3 + vh3 * (hCoeffs->delta22vh6 + vh * vh * (hCoeffs->delta22vh9 * vh))) +
                      Omega * (hCoeffs->delta22v5 * v2 + Omega * (hCoeffs->delta22v6 + hCoeffs->delta22v8 * v2));
            // FIXME
            //  if (debugPK){
            //      XLAL_PRINT_INFO("PK:: rho22v2 = %.12e, rho22v3 = %.12e,
            //      rho22v4 =
            //      %.12e,\n rho22v5 = %.16e, rho22v6 = %.16e, rho22v6LOG =
            //      %.16e, \n rho22v7 = %.12e, rho22v8 = %.16e, rho22v8LOG =
            //      %.16e, \n rho22v10 = %.16e, rho22v10LOG = %.16e\n, rho22v6 =
            //      %.12e, rho22v8 = %.12e, rho22v10 = %.12e\n",
            //          hCoeffs->rho22v2, hCoeffs->rho22v3, hCoeffs->rho22v4,
            //          hCoeffs->rho22v5, hCoeffs->rho22v6, hCoeffs->rho22v6l,
            //          hCoeffs->rho22v7, hCoeffs->rho22v8, hCoeffs->rho22v8l,
            //          hCoeffs->rho22v10, hCoeffs->rho22v10l,
            //          hCoeffs->rho22v6 + hCoeffs->rho22v6l * eulerlogxabs,
            //          hCoeffs->rho22v8 + hCoeffs->rho22v8l * eulerlogxabs,
            //          hCoeffs->rho22v10 + hCoeffs->rho22v10l * eulerlogxabs);}
            //     fTe = (CST_PI*(-7.*I*dr*(-15*pow(prDot,4)*pow(r,8) +
            //     180*pow(prDot,5)*pow(r,10) +
            //   210*pow(prDot,2)*pow(r,4)*(12 + 7*r*pow(dr,2)) +
            //   6*prDot*pow(r,3)*pow(dr,2)*(20 + 69*r*pow(dr,2)) -
            //   10*pow(prDot,3)*pow(r,6)*(96 + 137*r*pow(dr,2)) -
            //   3*(3840 + 280*r*pow(dr,2) + 89*pow(r,2)*pow(dr,4))) +
            //     2*pl*(-3633*pow(prDot,5)*pow(r,10) +
            //     3395*pow(prDot,6)*pow(r,12) -
            //   210*pow(prDot,3)*pow(r,6)*(14 + r*pow(dr,2)) +
            //   1260*pow(prDot,4)*pow(r,8)*(3 + 2*r*pow(dr,2)) +
            //   105*prDot*pow(r,3)*pow(dr,2)*(84 + 31*r*pow(dr,2)) -
            //   105*pow(prDot,2)*pow(r,5)*pow(dr,2)*(48 + 47*r*pow(dr,2)) +
            //   60*(672 + 8*pow(r,2)*pow(dr,4) +
            //   7*pow(r,3)*pow(dr,6)))))/(40320.*r); fTe = fTe - I*((omega*(-11
            //   + 12*CST_GAMMA - 6.*I*CST_PI + 36*CST_LN2
            //   + 36*Log(v)))/3.;
            // fTe = fTe - I * (Omega*(-11 + 12*CST_GAMMA - 6.*I*CST_PI +
            // 36*CST_LN2 + 12*log(Omega)))/3.;
            break;
        case 1: {
            deltalm =
                vh3 * (hCoeffs->delta21vh3 +
                       vh3 * (hCoeffs->delta21vh6 + vh * (hCoeffs->delta21vh7 + (hCoeffs->delta21vh9) * vh * vh))) +
                Omegav2 * (hCoeffs->delta21v5 + hCoeffs->delta21v7 * v2);
        }
        break;
        default:
            return CEV_FAILURE;
            break;
        }
        break;
    case 3:
        switch (m)
        {
        case 3:
            deltalm = vh3 * (hCoeffs->delta33vh3 + vh3 * (hCoeffs->delta33vh6 + hCoeffs->delta33vh9 * vh3)) +
                      hCoeffs->delta33v5 * v * v2 * v2 + hCoeffs->delta33v7 * v2 * v2 * v2 * v;
            // R.C: delta33v7 is set to 0, whoever is adding it here as a
            // coefficient is evil, TODO: double check that is zero and then
            // remove it RC: This terms are in Eq.A6 in
            // https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028
            // [arXiv:1803.10701]
            break;
        case 2:
            deltalm = vh3 * (hCoeffs->delta32vh3 +
                             vh * (hCoeffs->delta32vh4 + vh * vh * (hCoeffs->delta32vh6 + hCoeffs->delta32vh9 * vh3)));
            break;
        case 1:
            deltalm = vh3 * (hCoeffs->delta31vh3 +
                             vh3 * (hCoeffs->delta31vh6 + vh * (hCoeffs->delta31vh7 + hCoeffs->delta31vh9 * vh * vh))) +
                      hCoeffs->delta31v5 * v * v2 * v2;
            break;
        default:
            return CEV_FAILURE;
            break;
        }
        break;
    case 4:
        switch (m)
        {
        case 4:
            // RC: This terms are in Eq.A15 in
            // https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028
            // [arXiv:1803.10701]
            deltalm = vh3 * (hCoeffs->delta44vh3 + vh3 * (hCoeffs->delta44vh6 + vh3 * hCoeffs->delta44vh9)) +
                      hCoeffs->delta44v5 * v2 * v2 * v;
            // RC: This terms are in Eq.A8 in
            // https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028
            // [arXiv:1803.10701]
            break;
        case 3:
            deltalm = vh3 * (hCoeffs->delta43vh3 + vh * (hCoeffs->delta43vh4 + hCoeffs->delta43vh6 * vh * vh));
            break;
        case 2:
            deltalm = vh3 * (hCoeffs->delta42vh3 + hCoeffs->delta42vh6 * vh3);
            break;
        case 1:
            deltalm = vh3 * (hCoeffs->delta41vh3 + vh * (hCoeffs->delta41vh4 + hCoeffs->delta41vh6 * vh * vh));
            break;
        default:
            return CEV_FAILURE;
            break;
        }
        break;
    case 5:
        switch (m)
        {
        case 5:
            // RC: This terms are in Eq.A16 in
            // https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028
            // [arXiv:1803.10701]
            deltalm = vh3 * (hCoeffs->delta55vh3 + vh3 * (hCoeffs->delta55vh6 + vh3 * (hCoeffs->delta55vh9))) +
                      hCoeffs->delta55v5 * v2 * v2 * v;
            break;
        case 4:
            deltalm = vh3 * (hCoeffs->delta54vh3 + hCoeffs->delta54vh4 * vh);
            break;
        case 3:
            deltalm = hCoeffs->delta53vh3 * vh3;
            break;
        case 2:
            deltalm = vh3 * (hCoeffs->delta52vh3 + hCoeffs->delta52vh4 * vh);
            break;
        case 1:
            deltalm = hCoeffs->delta51vh3 * vh3;
            break;
        default:
            return CEV_FAILURE;
            break;
        }
        break;
    case 6:
        switch (m)
        {
        case 6:
            deltalm = hCoeffs->delta66vh3 * vh3;
            break;
        case 5:
            deltalm = hCoeffs->delta65vh3 * vh3;
            break;
        case 4:
            deltalm = hCoeffs->delta64vh3 * vh3;
            break;
        case 3:
            deltalm = hCoeffs->delta63vh3 * vh3;
            break;
        case 2:
            deltalm = hCoeffs->delta62vh3 * vh3;
            break;
        case 1:
            deltalm = hCoeffs->delta61vh3 * vh3;
            break;
        default:
            return CEV_FAILURE;
            break;
        }
        break;
    case 7:
        switch (m)
        {
        case 7:
            deltalm = hCoeffs->delta77vh3 * vh3;
            break;
        case 6:
            deltalm = 0.0;
            break;
        case 5:
            deltalm = hCoeffs->delta75vh3 * vh3;
            break;
        case 4:
            deltalm = 0.0;
            break;
        case 3:
            deltalm = hCoeffs->delta73vh3 * vh3;
            break;
        case 2:
            deltalm = 0.0;
            break;
        case 1:
            deltalm = hCoeffs->delta71vh3 * vh3;
            break;
        default:
            return CEV_FAILURE;
            break;
        }
        break;
    case 8:
        switch (m)
        {
        case 8:
        case 7:
        case 6:
        case 5:
        case 4:
        case 3:
        case 2:
        case 1:
            deltalm = 0.0;
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

    if (CalculateFacWaveformAmpResV2(hCoeffs, eta, l, m, v, vPhi, vh3, eulerlogxabs, 1 / sqrt(r), dr, ncrv,
                                     &rholmPwrl) != CEV_SUCCESS)
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Failure in CalculateFacWaveformAmpResV2");
        return CEV_FAILURE;
    }
    // if (r < 2.39988)
    //     print_debug("rholmPwrl = %.3e + i%.3e\n", creal(rholmPwrl),
    //     cimag(rholmPwrl));
    // if (m==1) {
    // 	printf("f21v1 = %.16f f21v3 = %.16f f21v4 = %.16f f21v5 = %.16f\n",
    // hCoeffs->f21v1, hCoeffs->f21v3, hCoeffs->f21v4, hCoeffs->f21v5);
    // }

    // if (r > 0.0 && debugPK) {
    //     XLAL_PRINT_INFO("YP::dynamics variables in waveform: %i, %i, r =
    //     %.12e, v = %.12e\n", l, m, r, v); XLAL_PRINT_INFO("rholm^l = %.16e +
    //     %.16e, Tlm = %.16e + i %.16e, \nSlm = %.16e, hNewton = %.16e + i
    //     %.16e, delta = %.16e\n", creal(rholmPwrl), cimag(rholmPwrl),
    //     creal(Tlm), cimag(Tlm), Slm, creal(hNewton), cimag(hNewton), 0.0);
    // }
    /* Put all factors in Eq. 17 together */
    *hlm = (fTe + Tlm * cexp(I * deltalm)) * Slm * (rholmPwrl);
    *hlm *= hNewton;
    // if (r > 8.5 && debugPK) {
    //     XLAL_PRINT_INFO("YP::FullWave: %.16e,%.16e, %.16e\n", creal(*hlm),
    //     cimag(*hlm), cabs(*hlm));
    // }
    return CEV_SUCCESS;
}

INT XLALSimIMRSpinEOBGetSASpinFactorizedWaveformV2(
    COMPLEX16 *hlm,       /**< OUTPUT, hlm waveforms */
    REAL8Vector *values,  /**< dyanmical variables: \f$(r,\phi,p_r,p_\phi)\f$ */
    const REAL8 v,        /**< velocity */
    const REAL8 dr,       /**< radial velocity */
    const REAL8 ncrv,     /**< angular velocity */
    const REAL8 prDot,    /** unused */
    const REAL8 Hreal,    /**< real Hamiltonian */
    const INT l,          /**< l mode index */
    const INT m,          /**< m mode index */
    SpinEOBParams *params /**< Spin EOB parameters */
)
{
    int debugPK = 0;
    /* Status of function calls */
    INT status;
    INT i;

    REAL8 eta;
    REAL8 r, pp, Omega, v2, Omegav2, vh, vh3, k, hathatk, eulerlogxabs;
    // pr
    //  REAL8 rcrossp_x, rcrossp_y, rcrossp_z;
    REAL8 Slm, deltalm;
    COMPLEX16 auxflm = 0.0;
    REAL8 hathatksq4, hathatk4pi, Tlmprefac, Tlmprodfac;
    COMPLEX16 Tlm, rholmPwrl, rholm;
    COMPLEX16 fTe = 0.0;
    COMPLEX16 hNewton;
    gsl_sf_result z2;

    /* Non-Keplerian velocity */
    REAL8 vPhi, vPhi2;

    /* Pre-computed coefficients */

    FacWaveformCoeffs *hCoeffs = params->hCoeffs;

    if (abs(m) > (INT)l)
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

    // rcrossp_x = cartvalues->data[1] * cartvalues->data[5] -
    // cartvalues->data[2]
    // * cartvalues->data[4]; rcrossp_y = cartvalues->data[2] *
    // cartvalues->data[3] - cartvalues->data[0] * cartvalues->data[5]; rcrossp_z
    // = cartvalues->data[0] * cartvalues->data[4] - cartvalues->data[1] *
    // cartvalues->data[3];

    // pp = sqrt(rcrossp_x * rcrossp_x + rcrossp_y * rcrossp_y + rcrossp_z *
    // rcrossp_z);

    v2 = v * v;
    Omega = v2 * v;
    Omegav2 = Omega * v2;
    vh3 = Hreal * Omega;
    vh = cbrt(vh3);
    eulerlogxabs = CST_GAMMA + log(2.0 * (REAL8)m * v);

    FacWaveformCoeffs oldCoeffs = *(params->hCoeffs); // RC: the function XLALSimIMRSpinPrecEOBNonKeplerCoeff
                                                      // is calculating again the coefficients hCoeffs, for
                                                      // the omega. These are different because
    // for the dynamics we are using different coefficients for the 21 mode. I
    // store the old coefficients and the I recover them ofter vPhi is computed.

    /* Calculate the non-Keplerian velocity */
    vPhi = XLALSimIMRSpinAlignedEOBNonKeplerCoeff(values->data, params);

    if (IS_REAL8_FAIL_NAN(vPhi))
    {
        return CEV_FAILURE;
    }
    vPhi = r * cbrt(vPhi);

    // if (debugPK)
    //     XLAL_PRINT_INFO("In XLALSimIMRSpinEOBFluxCalculateNewtonianMultipole,
    //     getting rW = %.12e\n", vPhi);
    vPhi *= Omega;
    vPhi2 = vPhi * vPhi;
    // *hCoeffs = oldCoeffs; //RC: Here I recover the old coefficients
    /*
     * Calculate the newtonian multipole, 1st term in Eq. 17, given by
     * Eq. A1
     */
    // if (debugPK) {
    //     XLAL_PRINT_INFO("\nValues inside
    //     XLALSimIMRSpinEOBFluxGetSpinFactorizedWaveform:\n"); for (i = 0; i <
    //     11; i++)
    //         XLAL_PRINT_INFO("values[%d] = %.12e\n", i, cartvalues->data[i]);

    //     XLAL_PRINT_INFO("Calculating hNewton, with v = %.12e, vPhi = %.12e, r
    //     =
    //     %.12e, Phi = %.12e, l = %d, m = %d\n",
    //             v, vPhi, r, values->data[1], (UINT4) l, (UINT4) m);
    // }
    status = XLALSimIMRSpinEOBCalculateNewtonianMultipole(&hNewton, vPhi2, r, values->data[1], (UINT)l, m, params);

    if (status != CEV_SUCCESS)
    {
        return CEV_FAILURE;
    }
    /* Calculate the source term, 2nd term in Eq. 17, given by Eq. A5, Hreal is
     * given by Eq.5 and Heff is in Eq.2 */
    if (((l + m) % 2) == 0)
    {
        Slm = (Hreal * Hreal - 1.) / (2. * eta) + 1.;
    }
    else
    {
        Slm = v * pp;
        // Slm = v * sqrt(rcrossp_x * rcrossp_x + rcrossp_y * rcrossp_y +
        // rcrossp_z
        // * rcrossp_z);
    }
    // if (debugPK)
    //     XLAL_PRINT_INFO("In XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform:
    //     Hreal = %e, Slm = %e, eta = %e\n", Hreal, Slm, eta);

    /*
     * Calculate the Tail term, 3rd term in Eq. 17,
     * given by Eq. A6, and Eq. (42) of
     * http://arxiv.org/pdf/1212.4357.pdf (or PRD 87 084035 (2013))
     */
    k = m * Omega;
    hathatk = Hreal * k;

    gsl_sf_result lnr1, arg1;
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
    Tlm = cexp((lnr1.val + CST_PI * hathatk) + I * (arg1.val + 2.0 * hathatk * log(4.0 * k / sqrt(CST_E))));
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
    //     XLAL_PRINT_INFO("Tlm = %e + i%e, |Tlm| = %.16e (should be %.16e)\n",
    //     creal(Tlm), cimag(Tlm), cabs(Tlm), Tlmold);
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
    REAL8 pl = sqrt(GET_MAX(r * prDot + 1. / r, 0.0));
    switch (l)
    {
    case 2:
        switch (abs(m))
        {
        case 2:
            deltalm = vh3 * (hCoeffs->delta22vh3 + vh3 * (hCoeffs->delta22vh6 + vh * vh * (hCoeffs->delta22vh9 * vh))) +
                      Omega * (hCoeffs->delta22v5 * v2 + Omega * (hCoeffs->delta22v6 + hCoeffs->delta22v8 * v2));
            // FIXME
            //  if (debugPK){
            //      XLAL_PRINT_INFO("PK:: rho22v2 = %.12e, rho22v3 = %.12e,
            //      rho22v4 =
            //      %.12e,\n rho22v5 = %.16e, rho22v6 = %.16e, rho22v6LOG =
            //      %.16e, \n rho22v7 = %.12e, rho22v8 = %.16e, rho22v8LOG =
            //      %.16e, \n rho22v10 = %.16e, rho22v10LOG = %.16e\n, rho22v6 =
            //      %.12e, rho22v8 = %.12e, rho22v10 = %.12e\n",
            //          hCoeffs->rho22v2, hCoeffs->rho22v3, hCoeffs->rho22v4,
            //          hCoeffs->rho22v5, hCoeffs->rho22v6, hCoeffs->rho22v6l,
            //          hCoeffs->rho22v7, hCoeffs->rho22v8, hCoeffs->rho22v8l,
            //          hCoeffs->rho22v10, hCoeffs->rho22v10l,
            //          hCoeffs->rho22v6 + hCoeffs->rho22v6l * eulerlogxabs,
            //          hCoeffs->rho22v8 + hCoeffs->rho22v8l * eulerlogxabs,
            //          hCoeffs->rho22v10 + hCoeffs->rho22v10l * eulerlogxabs);}
            //     fTe = (CST_PI*(-7.*I*dr*(-15*pow(prDot,4)*pow(r,8) +
            //     180*pow(prDot,5)*pow(r,10) +
            //   210*pow(prDot,2)*pow(r,4)*(12 + 7*r*pow(dr,2)) +
            //   6*prDot*pow(r,3)*pow(dr,2)*(20 + 69*r*pow(dr,2)) -
            //   10*pow(prDot,3)*pow(r,6)*(96 + 137*r*pow(dr,2)) -
            //   3*(3840 + 280*r*pow(dr,2) + 89*pow(r,2)*pow(dr,4))) +
            //     2*pl*(-3633*pow(prDot,5)*pow(r,10) +
            //     3395*pow(prDot,6)*pow(r,12) -
            //   210*pow(prDot,3)*pow(r,6)*(14 + r*pow(dr,2)) +
            //   1260*pow(prDot,4)*pow(r,8)*(3 + 2*r*pow(dr,2)) +
            //   105*prDot*pow(r,3)*pow(dr,2)*(84 + 31*r*pow(dr,2)) -
            //   105*pow(prDot,2)*pow(r,5)*pow(dr,2)*(48 + 47*r*pow(dr,2)) +
            //   60*(672 + 8*pow(r,2)*pow(dr,4) +
            //   7*pow(r,3)*pow(dr,6)))))/(40320.*r); fTe = fTe - I*((omega*(-11
            //   + 12*CST_GAMMA - 6.*I*CST_PI + 36*CST_LN2
            //   + 36*Log(v)))/3.;
            // fTe = fTe - I * (Omega*(-11 + 12*CST_GAMMA - 6.*I*CST_PI +
            // 36*CST_LN2 + 12*log(Omega)))/3.;
            break;
        case 1: {
            deltalm =
                vh3 * (hCoeffs->delta21vh3 +
                       vh3 * (hCoeffs->delta21vh6 + vh * (hCoeffs->delta21vh7 + (hCoeffs->delta21vh9) * vh * vh))) +
                Omegav2 * (hCoeffs->delta21v5 + hCoeffs->delta21v7 * v2);
        }
        break;
        default:
            return CEV_FAILURE;
            break;
        }
        break;
    case 3:
        switch (m)
        {
        case 3:
            deltalm = vh3 * (hCoeffs->delta33vh3 + vh3 * (hCoeffs->delta33vh6 + hCoeffs->delta33vh9 * vh3)) +
                      hCoeffs->delta33v5 * v * v2 * v2 + hCoeffs->delta33v7 * v2 * v2 * v2 * v;
            // R.C: delta33v7 is set to 0, whoever is adding it here as a
            // coefficient is evil, TODO: double check that is zero and then
            // remove it RC: This terms are in Eq.A6 in
            // https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028
            // [arXiv:1803.10701]
            break;
        case 2:
            deltalm = vh3 * (hCoeffs->delta32vh3 +
                             vh * (hCoeffs->delta32vh4 + vh * vh * (hCoeffs->delta32vh6 + hCoeffs->delta32vh9 * vh3)));
            break;
        case 1:
            deltalm = vh3 * (hCoeffs->delta31vh3 +
                             vh3 * (hCoeffs->delta31vh6 + vh * (hCoeffs->delta31vh7 + hCoeffs->delta31vh9 * vh * vh))) +
                      hCoeffs->delta31v5 * v * v2 * v2;
            break;
        default:
            return CEV_FAILURE;
            break;
        }
        break;
    case 4:
        switch (m)
        {
        case 4:
            // RC: This terms are in Eq.A15 in
            // https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028
            // [arXiv:1803.10701]
            deltalm = vh3 * (hCoeffs->delta44vh3 + vh3 * (hCoeffs->delta44vh6 + vh3 * hCoeffs->delta44vh9)) +
                      hCoeffs->delta44v5 * v2 * v2 * v;
            // RC: This terms are in Eq.A8 in
            // https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028
            // [arXiv:1803.10701]
            break;
        case 3:
            deltalm = vh3 * (hCoeffs->delta43vh3 + vh * (hCoeffs->delta43vh4 + hCoeffs->delta43vh6 * vh * vh));
            break;
        case 2:
            deltalm = vh3 * (hCoeffs->delta42vh3 + hCoeffs->delta42vh6 * vh3);
            break;
        case 1:
            deltalm = vh3 * (hCoeffs->delta41vh3 + vh * (hCoeffs->delta41vh4 + hCoeffs->delta41vh6 * vh * vh));
            break;
        default:
            return CEV_FAILURE;
            break;
        }
        break;
    case 5:
        switch (m)
        {
        case 5:
            // RC: This terms are in Eq.A16 in
            // https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028
            // [arXiv:1803.10701]
            deltalm = vh3 * (hCoeffs->delta55vh3 + vh3 * (hCoeffs->delta55vh6 + vh3 * (hCoeffs->delta55vh9))) +
                      hCoeffs->delta55v5 * v2 * v2 * v;
            break;
        case 4:
            deltalm = vh3 * (hCoeffs->delta54vh3 + hCoeffs->delta54vh4 * vh);
            break;
        case 3:
            deltalm = hCoeffs->delta53vh3 * vh3;
            break;
        case 2:
            deltalm = vh3 * (hCoeffs->delta52vh3 + hCoeffs->delta52vh4 * vh);
            break;
        case 1:
            deltalm = hCoeffs->delta51vh3 * vh3;
            break;
        default:
            return CEV_FAILURE;
            break;
        }
        break;
    case 6:
        switch (m)
        {
        case 6:
            deltalm = hCoeffs->delta66vh3 * vh3;
            break;
        case 5:
            deltalm = hCoeffs->delta65vh3 * vh3;
            break;
        case 4:
            deltalm = hCoeffs->delta64vh3 * vh3;
            break;
        case 3:
            deltalm = hCoeffs->delta63vh3 * vh3;
            break;
        case 2:
            deltalm = hCoeffs->delta62vh3 * vh3;
            break;
        case 1:
            deltalm = hCoeffs->delta61vh3 * vh3;
            break;
        default:
            return CEV_FAILURE;
            break;
        }
        break;
    case 7:
        switch (m)
        {
        case 7:
            deltalm = hCoeffs->delta77vh3 * vh3;
            break;
        case 6:
            deltalm = 0.0;
            break;
        case 5:
            deltalm = hCoeffs->delta75vh3 * vh3;
            break;
        case 4:
            deltalm = 0.0;
            break;
        case 3:
            deltalm = hCoeffs->delta73vh3 * vh3;
            break;
        case 2:
            deltalm = 0.0;
            break;
        case 1:
            deltalm = hCoeffs->delta71vh3 * vh3;
            break;
        default:
            return CEV_FAILURE;
            break;
        }
        break;
    case 8:
        switch (m)
        {
        case 8:
        case 7:
        case 6:
        case 5:
        case 4:
        case 3:
        case 2:
        case 1:
            deltalm = 0.0;
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

    if (CalculateFacWaveformAmpResV2(hCoeffs, eta, l, m, v, vPhi, vh3, eulerlogxabs, 1 / sqrt(r), dr, ncrv,
                                     &rholmPwrl) != CEV_SUCCESS)
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Failure in CalculateFacWaveformAmpResV2");
        return CEV_FAILURE;
    }
    // if (r < 2.39988)
    //     print_debug("rholmPwrl = %.3e + i%.3e\n", creal(rholmPwrl),
    //     cimag(rholmPwrl));
    // if (m==1) {
    // 	printf("f21v1 = %.16f f21v3 = %.16f f21v4 = %.16f f21v5 = %.16f\n",
    // hCoeffs->f21v1, hCoeffs->f21v3, hCoeffs->f21v4, hCoeffs->f21v5);
    // }

    // if (r > 0.0 && debugPK) {
    //     XLAL_PRINT_INFO("YP::dynamics variables in waveform: %i, %i, r =
    //     %.12e, v = %.12e\n", l, m, r, v); XLAL_PRINT_INFO("rholm^l = %.16e +
    //     %.16e, Tlm = %.16e + i %.16e, \nSlm = %.16e, hNewton = %.16e + i
    //     %.16e, delta = %.16e\n", creal(rholmPwrl), cimag(rholmPwrl),
    //     creal(Tlm), cimag(Tlm), Slm, creal(hNewton), cimag(hNewton), 0.0);
    // }
    /* Put all factors in Eq. 17 together */
    *hlm = (fTe + Tlm * cexp(I * deltalm)) * Slm * (rholmPwrl);
    *hlm *= hNewton;
    // if (r > 8.5 && debugPK) {
    //     XLAL_PRINT_INFO("YP::FullWave: %.16e,%.16e, %.16e\n", creal(*hlm),
    //     cimag(*hlm), cabs(*hlm));
    // }
    return CEV_SUCCESS;
}

INT XLALSimIMRSpinEOBGetAmplitudeResidualPrecV2(
    COMPLEX16 *rholmpwrl, REAL8Vector *values, /**< dyanmical variables: \f$(r,\phi,p_r,p_\phi)\f$ */
    REAL8Vector *cartvalues,                   /**< dyanmical variables */
    const REAL8 v, const REAL8 dr, const REAL8 ncrv, const REAL8 Hreal, const INT modeL, const INT modeM,
    SpinEOBParams *params)
{
    INT4 i = 0;
    REAL8 eta;
    REAL8 eulerlogxabs;
    REAL8 rholm;
    REAL8 v2 = v * v, Omega, vh3, r;
    COMPLEX16 auxflm = 0.0;
    COMPLEX16 rholmPwrl;
    FacWaveformCoeffs *hCoeffs = params->hCoeffs;
    eulerlogxabs = CST_GAMMA + log(2.0 * (REAL8)modeM * v);
    REAL8 vPhi;

    if (abs(modeM) > (INT4)modeL)
    {
        return CEV_FAILURE;
    }
    if (modeM == 0)
    {
        return CEV_FAILURE;
    }
    eta = params->eta;

    /* Check our eta was sensible */
    if (eta > 0.25 && eta < 0.25 + 1e-4)
    {
        eta = 0.25;
    }
    if (eta > 0.25)
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Eta seems to be > 0.25 - this isn't allowed!");
        return CEV_FAILURE;
    }
    r = values->data[0];
    Omega = v2 * v;
    vh3 = Hreal * Omega;
    // vh = cbrt(vh3);

    vPhi = params->alignedSpins ? XLALSimIMRSpinAlignedEOBNonKeplerCoeff(values->data, params)
                                : XLALSimIMRSpinPrecEOBNonKeplerCoeff(cartvalues->data, params);
    /* Calculate the non-Keplerian velocity */

    vPhi = r * cbrt(vPhi);

    // if (debugPK)
    // 	XLAL_PRINT_INFO("In XLALSimIMRSpinEOBFluxCalculateNewtonianMultipole,
    // getting rW = %.12e\n", vPhi);
    vPhi *= Omega;

#if 0
    switch (modeL)
    {
        case 2:
    switch (abs (modeM))
    {
    case 2:
    rholm =
        1. + v2 * (hCoeffs->rho22v2 +
            v * (hCoeffs->rho22v3 +
                v * (hCoeffs->rho22v4 +
                v * (hCoeffs->rho22v5 +
                    v * (hCoeffs->rho22v6 +
                    hCoeffs->rho22v6l * eulerlogxabs +
                    v * (hCoeffs->rho22v7 +
                        v * (hCoeffs->rho22v8 +
                            hCoeffs->rho22v8l *
                            eulerlogxabs +
                            (hCoeffs->rho22v10 +
                            hCoeffs->rho22v10l *
                            eulerlogxabs) *
                            v2)))))));
    break;
    case 1:
    {
        rholm =
        1. + v * (hCoeffs->rho21v1 +
            v * (hCoeffs->rho21v2 +
                v * (hCoeffs->rho21v3 +
                v * (hCoeffs->rho21v4 +
                    v * (hCoeffs->rho21v5 +
                        v * (hCoeffs->rho21v6 +
                        hCoeffs->rho21v6l *
                        eulerlogxabs +
                        v * (hCoeffs->rho21v7 +
                            hCoeffs->rho21v7l *
                            eulerlogxabs +
                            v * (hCoeffs->rho21v8 +
                            hCoeffs->rho21v8l *
                            eulerlogxabs +
                            (hCoeffs->
                                rho21v10 +
                                hCoeffs->
                                rho21v10l *
                                eulerlogxabs) *
                            v2))))))));
                auxflm = v * (hCoeffs->f21v1 + v2 * (hCoeffs->f21v3 + v * hCoeffs->f21v4 + v2 * (hCoeffs->f21v5 + v * hCoeffs->f21v6)));
            }
    break;
    default:
        return CEV_FAILURE;
    break;
    }
    break;
    case 3:
    switch (modeM)
    {
    case 3:
        rholm =
            1. + v2 * (hCoeffs->rho33v2 +
                        v * (hCoeffs->rho33v3 +
                            v * (hCoeffs->rho33v4 +
                                v * (hCoeffs->rho33v5 +
                                        v * (hCoeffs->rho33v6 +
                                            hCoeffs->rho33v6l * eulerlogxabs +
                                            v * (hCoeffs->rho33v7 +
                                                v * (hCoeffs->rho33v8 +
                                                hCoeffs->rho33v8l *
                                                eulerlogxabs)))))));
        auxflm = v * v2 * hCoeffs->f33v3;
    break;
    case 2:
    rholm =
        1. + v * (hCoeffs->rho32v +
            v * (hCoeffs->rho32v2 +
            v * (hCoeffs->rho32v3 +
                v * (hCoeffs->rho32v4 +
                    v * (hCoeffs->rho32v5 +
                    v * (hCoeffs->rho32v6 +
                        hCoeffs->rho32v6l *
                        eulerlogxabs +
                        (hCoeffs->rho32v8 +
                        hCoeffs->rho32v8l *
                        eulerlogxabs) * v2))))));
    break;
    case 1:
    rholm =
        1. + v2 * (hCoeffs->rho31v2 +
            v * (hCoeffs->rho31v3 +
                v * (hCoeffs->rho31v4 +
                v * (hCoeffs->rho31v5 +
                    v * (hCoeffs->rho31v6 +
                    hCoeffs->rho31v6l * eulerlogxabs +
                    v * (hCoeffs->rho31v7 +
                        (hCoeffs->rho31v8 +
                        hCoeffs->rho31v8l *
                        eulerlogxabs) * v))))));
    auxflm = v * v2 * hCoeffs->f31v3;
    break;
    default:
        return CEV_FAILURE;
    break;
    }
    break;
    case 4:
    switch (modeM)
    {
    case 4:
    rholm = 1. + v2 * (hCoeffs->rho44v2
            + v * (hCoeffs->rho44v3 + v * (hCoeffs->rho44v4
                    +
                    v *
                    (hCoeffs->
                    rho44v5 +
                    (hCoeffs->
                        rho44v6 +
                        hCoeffs->
                        rho44v6l *
                        eulerlogxabs) *
                    v))));

    break;
    case 3:
    rholm =
        1. + v * (hCoeffs->rho43v +
            v * (hCoeffs->rho43v2 +
            v2 * (hCoeffs->rho43v4 +
                v * (hCoeffs->rho43v5 +
                    (hCoeffs->rho43v6 +
                    hCoeffs->rho43v6l * eulerlogxabs) *
                    v))));
    auxflm = v * hCoeffs->f43v;
    break;
    case 2:
    rholm = 1. + v2 * (hCoeffs->rho42v2
                + v * (hCoeffs->rho42v3 +
                    v * (hCoeffs->rho42v4 +
                    v * (hCoeffs->rho42v5 +
                        (hCoeffs->rho42v6 +
                        hCoeffs->rho42v6l *
                        eulerlogxabs) * v))));
    break;
    case 1:
    rholm =
        1. + v * (hCoeffs->rho41v +
            v * (hCoeffs->rho41v2 +
            v2 * (hCoeffs->rho41v4 +
                v * (hCoeffs->rho41v5 +
                    (hCoeffs->rho41v6 +
                    hCoeffs->rho41v6l * eulerlogxabs) *
                    v))));
    auxflm = v * hCoeffs->f41v;
    break;
    default:
        return CEV_FAILURE;
    break;
    }
    break;
    case 5:
    switch (modeM)
    {
    case 5:
    //RC: This terms are in Eq.A9 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
    rholm =
        1. + v2 * (hCoeffs->rho55v2 +
            v * (hCoeffs->rho55v3 +
            v * (hCoeffs->rho55v4 +
        v * (hCoeffs->rho55v5 +
            v * (hCoeffs->rho55v6 + hCoeffs->rho55v6l * eulerlogxabs +
            v2 * (hCoeffs->rho55v8 + hCoeffs->rho55v8l * eulerlogxabs +
        v2 * (hCoeffs->rho55v10 + hCoeffs->rho55v10l * eulerlogxabs )))))));
    //RC: This terms are in Eq.A12 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
    auxflm = v2 * v * (hCoeffs->f55v3 + v * (hCoeffs->f55v4));

    break;
    case 4:
    rholm = 1. + v2 * (hCoeffs->rho54v2 + v * (hCoeffs->rho54v3
                            + hCoeffs->rho54v4 * v));
    break;
    case 3:
    rholm = 1. + v2 * (hCoeffs->rho53v2
                + v * (hCoeffs->rho53v3 +
                    v * (hCoeffs->rho53v4 +
                    hCoeffs->rho53v5 * v)));
    break;
    case 2:
    rholm = 1. + v2 * (hCoeffs->rho52v2 + v * (hCoeffs->rho52v3
                            + hCoeffs->rho52v4 * v));
    break;
    case 1:
    rholm = 1. + v2 * (hCoeffs->rho51v2
                + v * (hCoeffs->rho51v3 +
                    v * (hCoeffs->rho51v4 +
                    hCoeffs->rho51v5 * v)));
    break;
    default:
        return CEV_FAILURE;
    break;
    }
    break;
    case 6:
    switch (modeM)
    {
    case 6:
    rholm = 1. + v2 * (hCoeffs->rho66v2 + v * (hCoeffs->rho66v3
                            + hCoeffs->rho66v4 * v));
    break;
    case 5:
    rholm = 1. + v2 * (hCoeffs->rho65v2 + hCoeffs->rho65v3 * v);
    break;
    case 4:
    rholm = 1. + v2 * (hCoeffs->rho64v2 + v * (hCoeffs->rho64v3
                            + hCoeffs->rho64v4 * v));
    break;
    case 3:
    rholm = 1. + v2 * (hCoeffs->rho63v2 + hCoeffs->rho63v3 * v);
    break;
    case 2:
    rholm = 1. + v2 * (hCoeffs->rho62v2 + v * (hCoeffs->rho62v3
                            + hCoeffs->rho62v4 * v));
    break;
    case 1:
    rholm = 1. + v2 * (hCoeffs->rho61v2 + hCoeffs->rho61v3 * v);
    break;
    default:
        return CEV_FAILURE;
    break;
    }
    break;
    case 7:
    switch (modeM)
    {
    case 7:
    rholm = 1. + v2 * (hCoeffs->rho77v2 + hCoeffs->rho77v3 * v);
    break;
    case 6:
    rholm = 1. + hCoeffs->rho76v2 * v2;
    break;
    case 5:
    rholm = 1. + v2 * (hCoeffs->rho75v2 + hCoeffs->rho75v3 * v);
    break;
    case 4:
    rholm = 1. + hCoeffs->rho74v2 * v2;
    break;
    case 3:
    rholm = 1. + v2 * (hCoeffs->rho73v2 + hCoeffs->rho73v3 * v);
    break;
    case 2:
    rholm = 1. + hCoeffs->rho72v2 * v2;
    break;
    case 1:
    rholm = 1. + v2 * (hCoeffs->rho71v2 + hCoeffs->rho71v3 * v);
    break;
    default:
        return CEV_FAILURE;
    break;
    }
    break;
    case 8:
    switch (modeM)
    {
    case 8:
    rholm = 1. + hCoeffs->rho88v2 * v2;
    break;
    case 7:
    rholm = 1. + hCoeffs->rho87v2 * v2;
    break;
    case 6:
    rholm = 1. + hCoeffs->rho86v2 * v2;
    break;
    case 5:
    rholm = 1. + hCoeffs->rho85v2 * v2;
    break;
    case 4:
    rholm = 1. + hCoeffs->rho84v2 * v2;
    break;
    case 3:
    rholm = 1. + hCoeffs->rho83v2 * v2;
    break;
    case 2:
    rholm = 1. + hCoeffs->rho82v2 * v2;
    break;
    case 1:
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

    /* Raise rholm to the lth power */
    rholmPwrl = 1.0;
    i = modeL;
    while (i--)
    {
        rholmPwrl *= rholm;
    }
    /* In the equal-mass odd m case, there is no contribution from nonspin terms,
    * and the only contribution comes from the auxflm term that is proportional to chiA (asymmetric spins).
    * In this case, we must ignore the nonspin terms directly, since the leading term defined by
    * CalculateThisMultipolePrefix in LALSimIMREOBNewtonianMultipole.c is not zero (see comments there).
    */
    if (eta == 0.25 && modeM % 2)
    {
        rholmPwrl = auxflm;
    }
    else
    {
        rholmPwrl += auxflm;
    }
    *rholmpwrl = rholmPwrl;
#endif
    if (CalculateFacWaveformAmpResV2(hCoeffs, eta, modeL, modeM, v, vPhi, vh3, eulerlogxabs, 1 / sqrt(r), dr, ncrv,
                                     &rholmPwrl) != CEV_SUCCESS)
        return CEV_FAILURE;
    *rholmpwrl = rholmPwrl;
    return CEV_SUCCESS;
}

INT XLALSimIMRSpinEOBSAGetAmplitudeResidualPrecV2(
    COMPLEX16 *rholmpwrl, REAL8Vector *values, /**< dyanmical variables: \f$(r,\phi,p_r,p_\phi)\f$ */
    const REAL8 v, const REAL8 dr, const REAL8 ncrv, const REAL8 Hreal, const INT modeL, const INT modeM,
    SpinEOBParams *params)
{
    INT4 i = 0;
    REAL8 eta;
    REAL8 eulerlogxabs;
    REAL8 rholm;
    REAL8 v2 = v * v, Omega, vh3, r;
    COMPLEX16 auxflm = 0.0;
    COMPLEX16 rholmPwrl;
    FacWaveformCoeffs *hCoeffs = params->hCoeffs;
    eulerlogxabs = CST_GAMMA + log(2.0 * (REAL8)modeM * v);
    REAL8 vPhi;

    if (abs(modeM) > (INT4)modeL)
    {
        return CEV_FAILURE;
    }
    if (modeM == 0)
    {
        return CEV_FAILURE;
    }
    eta = params->eta;

    /* Check our eta was sensible */
    if (eta > 0.25 && eta < 0.25 + 1e-4)
    {
        eta = 0.25;
    }
    if (eta > 0.25)
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Eta seems to be > 0.25 - this isn't allowed!");
        return CEV_FAILURE;
    }
    r = values->data[0];
    Omega = v2 * v;
    vh3 = Hreal * Omega;
    // vh = cbrt(vh3);

    vPhi = XLALSimIMRSpinAlignedEOBNonKeplerCoeff(values->data, params);
    /* Calculate the non-Keplerian velocity */

    vPhi = r * cbrt(vPhi);

    // if (debugPK)
    // 	XLAL_PRINT_INFO("In XLALSimIMRSpinEOBFluxCalculateNewtonianMultipole,
    // getting rW = %.12e\n", vPhi);
    vPhi *= Omega;

    if (CalculateFacWaveformAmpResV2(hCoeffs, eta, modeL, modeM, v, vPhi, vh3, eulerlogxabs, 1 / sqrt(r), dr, ncrv,
                                     &rholmPwrl) != CEV_SUCCESS)
        return CEV_FAILURE;
    *rholmpwrl = rholmPwrl;
    return CEV_SUCCESS;
}

REAL8
XLALInspiralPrecSpinFactorizedFluxV2(REAL8Vector *polvalues,          /**< \f$(r,\phi,p_r,p_\phi)\f$ */
                                     REAL8Vector *values,             /**< dynamical variables */
                                     EOBNonQCCoeffs *nqcCoeffs,       /**< pre-computed NQC coefficients */
                                     const REAL8 omega,               /**< orbital frequency */
                                     SpinEOBParams *ak,               /**< physical parameters */
                                     const REAL8 H,                   /**< real Hamiltonian */
                                     const INT4 lMax,                 /**< upper limit of the summation over l */
                                     const UINT SpinAlignedEOBversion /**< 1 for SEOBNRv1, 2 for SEOBNRv2 */
)
{
    // int	debugPK = 0;
    int i = 0;
    double radius =
        sqrt(values->data[0] * values->data[0] + values->data[1] * values->data[1] + values->data[2] * values->data[2]);
    if (radius < 1.)
    {
        return 0.;
    }
    if (1)
    {
        for (i = 0; i < 4; i++)
            if (isnan(polvalues->data[i]))
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "polvalues %3.10f %3.10f %3.10f %3.10f", polvalues->data[0],
                               polvalues->data[1], polvalues->data[2], polvalues->data[3]);
                PRINT_LOG_INFO(LOG_CRITICAL, "nan polvalues:  %3.10f %3.10f %3.10f %3.10f", polvalues->data[0],
                               polvalues->data[1], polvalues->data[2], polvalues->data[3]);
                return REAL8_FAIL_NAN;
            }

        for (i = 0; i < 12; i++)
            if (isnan(values->data[i]))
            {
                PRINT_LOG_INFO(LOG_CRITICAL,
                               "values %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f "
                               "%3.10f %3.10f %3.10f %3.10f %3.10f %3.10f",
                               values->data[0], values->data[1], values->data[2], values->data[3], values->data[4],
                               values->data[5], values->data[6], values->data[7], values->data[8], values->data[9],
                               values->data[10], values->data[11]);
                PRINT_LOG_INFO(LOG_CRITICAL,
                               "nan  in input values:  %3.10f %3.10f %3.10f "
                               "%3.10f %3.10f %3.10f "
                               "%3.10f %3.10f %3.10f %3.10f %3.10f %3.10f",
                               values->data[0], values->data[1], values->data[2], values->data[3], values->data[4],
                               values->data[5], values->data[6], values->data[7], values->data[8], values->data[9],
                               values->data[10], values->data[11]);
                return REAL8_FAIL_NAN;
            }
    }

    REAL8 flux = 0.0;
    REAL8 v;
    REAL8 omegaSq;
    COMPLEX16 hLM;
    INT4 l, m;

    // EOBNonQCCoeffs nqcCoeffs;
    //  if (lMax < 2)
    //  {
    //      XLAL_ERROR_REAL8(XLAL_EINVAL);
    //  }
    /* Omega is the derivative of phi */
    omegaSq = omega * omega;

    v = cbrt(omega);

#if 0
    /* Update the factorized multipole coefficients, w.r.t. new spins */
    if (0) {		/* {{{ */
        XLAL_PRINT_INFO("\nValues inside Flux:\n");
        for (i = 0; i < 11; i++)
            XLAL_PRINT_INFO("values[%d] = %.12e\n", i, values->data[i]);
        /*
            * Assume that initial conditions are available at this
            * point, to compute the chiS and chiA parameters. Calculate
            * the values of chiS and chiA, as given in Eq.16 of
            * Precessing EOB paper Pan et.al. arXiv:1307.6232 (or PRD 89, 084006 (2014)). Assuming \vec{L} to be pointing in
            * the direction of \vec{r}\times\vec{p}
            */
        REAL8		rcrossp  [3], rcrosspMag, s1dotL, s2dotL;
        REAL8		chiS    , chiA, tplspin;

        rcrossp[0] = values->data[1] * values->data[5] - values->data[2] * values->data[4];
        rcrossp[1] = values->data[2] * values->data[3] - values->data[0] * values->data[5];
        rcrossp[2] = values->data[0] * values->data[4] - values->data[1] * values->data[3];
        rcrosspMag = sqrt(rcrossp[0] * rcrossp[0] + rcrossp[1] * rcrossp[1] +
                    rcrossp[2] * rcrossp[2]);

        rcrossp[0] /= rcrosspMag;
        rcrossp[1] /= rcrosspMag;
        rcrossp[2] /= rcrosspMag;

        s1dotL = values->data[6] * rcrossp[0] + values->data[7] * rcrossp[1]
            + values->data[8] * rcrossp[2];
        s2dotL = values->data[9] * rcrossp[0] + values->data[10] * rcrossp[1]
            + values->data[11] * rcrossp[2];

        chiS = 0.5 * (s1dotL + s2dotL);
        chiA = 0.5 * (s1dotL - s2dotL);

        /*
            * Compute the test-particle limit spin of the deformed-Kerr
            * background
            */
        switch (SpinAlignedEOBversion) {
        case 1:
            tplspin = 0.0;
            break;
        case 2:
            tplspin = (1. - 2. * ak->eobParams->eta) * chiS + (ak->eobParams->m1
                                        - ak->eobParams->m2) / (ak->eobParams->m1 + ak->eobParams->m2) * chiA;
            break;
        default:
            XLALPrintError("XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
            break;
        }

        /* ************************************************* */
        /* Re-Populate the Waveform structures               */
        /* ************************************************* */

        /* Re-compute the spinning coefficients for hLM */
        //debugPK
            XLAL_PRINT_INFO("Re-calculating waveform coefficients in the Flux function with chiS, chiA = %e, %e!\n", chiS, chiA);
        chiS = 0.3039435650957116;
        chiA = -0.2959424290852973;
        XLAL_PRINT_INFO("Changed them to the correct values = %e, %e!\n", chiS, chiA);

        if (ak->alignedSpins==1) {
        if (XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients(ak->eobParams->hCoeffs,
            ak->eobParams->m1, ak->eobParams->m2, ak->eobParams->eta,
                                    tplspin, chiS, chiA, SpinAlignedEOBversion) == XLAL_FAILURE) {
            XLALDestroyREAL8Vector(values);
            XLAL_ERROR(XLAL_EFUNC);
        }
        }
        else {
            if (XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients(ak->eobParams->hCoeffs,
                                                                ak->eobParams->m1, ak->eobParams->m2, ak->eobParams->eta,
                                                                tplspin, chiS, chiA, 3) == XLAL_FAILURE) {
                XLALDestroyREAL8Vector(values);
                XLAL_ERROR(XLAL_EFUNC);
            }
        }
    }			/* }}} */
#endif
    // XLAL_PRINT_INFO("v = %.16e\n", v);
    COMPLEX16 hLMTab[lMax + 1][lMax + 1];
    if (XLALSimIMRSpinEOBFluxGetPrecSpinFactorizedWaveform(&hLMTab[0][0], polvalues, values, v, H, lMax, ak) ==
        CEV_FAILURE)
    {
        return REAL8_FAIL_NAN;
    }
    for (l = 2; l <= lMax; l++)
    {
        for (m = 1; m <= l; m++)
        {

            // if (debugPK)
            //     XLAL_PRINT_INFO("\nGetting (%d, %d) mode for flux!\n", l, m);
            // XLAL_PRINT_INFO("Stas, computing the waveform l = %d, m =%d\n", l,
            // m);
            hLM = hLMTab[l][m];
            // XLAL_PRINT_INFO("Stas: done\n");
            /*
             * For the 2,2 mode, we apply NQC correction to the
             * flux
             */
#if 0
            if (l == 2 && m == 2)
            {
                COMPLEX16	hNQC;
                /*
                    * switch ( SpinAlignedEOBversion ) { case 1:
                    * XLALSimIMRGetEOBCalibratedSpinNQC(
                    * &nqcCoeffs, l, m, ak->eobParams->eta,
                    * ak->a ); break; case 2: //
                    * XLALSimIMRGetEOBCalibratedSpinNQCv2(
                    * &nqcCoeffs, l, m, ak->eobParams->eta,
                    * ak->a );
                    * XLALSimIMRGetEOBCalibratedSpinNQC3D(
                    * &nqcCoeffs, l, m, ak->eobParams->eta,
                    * ak->a, (ak->chi1 - ak->chi2)/2. ); break;
                    * default: XLALPrintError( "XLAL Error - %s:
                    * Unknown SEOBNR version!\nAt present only
                    * v1 and v2 are available.\n", __func__);
                    * XLAL_ERROR( XLAL_EINVAL ); break; }
                    */
                // if (debugPK)
                //     XLAL_PRINT_INFO("\tl = %d, m = %d, NQC: a1 = %.16e, a2 = %.16e, a3 = %.16e, a3S = %.16e, a4 = %.16e, a5 = %.16e\n\tb1 = %.16e, b2 = %.16e, b3 = %.16e, b4 = %.16e\n",
                //             2, 2, nqcCoeffs->a1, nqcCoeffs->a2, nqcCoeffs->a3, nqcCoeffs->a3S, nqcCoeffs->a4, nqcCoeffs->a5,
                //             nqcCoeffs->b1, nqcCoeffs->b2, nqcCoeffs->b3, nqcCoeffs->b4);
                XLALSimIMREOBNonQCCorrection(&hNQC, polvalues, omega, nqcCoeffs);
                // if (debugPK)
                //     XLAL_PRINT_INFO("\tl = %d, m = %d, hNQC = %.16e + i%.16e, |hNQC| = %.16e\n", l, m,
                //             creal(hNQC), cimag(hNQC), sqrt(creal(hNQC) * creal(hNQC) + cimag(hLM) * cimag(hLM)));

                // if((m * m) * omegaSq * (creal(hLM) * creal(hLM) + cimag(hLM) * cimag(hLM)) > 5.)
                // {

                // XLAL_PRINT_INFO("\tl = %d, m = %d, mag(hLM) = %.17e, mag(hNQC) = %.17e, omega = %.16e\n",
                //     l, m, sqrt(creal(hLM) * creal(hLM) + cimag(hLM) * cimag(hLM)),
                //     sqrt(creal(hNQC) * creal(hNQC) + cimag(hNQC) * cimag(hNQC)), omega);

                // XLAL_PRINT_INFO("XLALInspiralPrecSpinFactorizedFlux (from input)::values %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f\n", values->data[0], values->data[1], values->data[2], values->data[3], values->data[4], values->data[5], values->data[6], values->data[7], values->data[8], values->data[9], values->data[10], values->data[11]);
                // }

                /* Eq. 16 */
                //FIXME
                hLM *= hNQC;
            }
#endif
            // if (debugPK)
            //     XLAL_PRINT_INFO("\tl = %d, m = %d, mag(hLM) = %.17e, omega =
            //     %.16e\n", l, m, sqrt(creal(hLM) * creal(hLM) + cimag(hLM) *
            //     cimag(hLM)), omega);

            /* Eq. 13 */
            flux += (REAL8)(m * m) * omegaSq * (creal(hLM) * creal(hLM) + cimag(hLM) * cimag(hLM));
        }
    }
    // if( (omegaSq > 1 || flux > 5) )
    // {
    //     if(debugPK) {
    //         XLAL_PRINT_INFO("In XLALInspiralPrecSpinFactorizedFlux: omegaSq =
    //         %3.12f, FLUX = %3.12f, r = %3.12f\n",
    //                 omegaSq, flux,radius);
    //     }
    //     flux = 0.;
    // }

    // if (debugPK)
    //     XLAL_PRINT_INFO("\tStas, FLUX = %.16e\n", flux * LAL_1_PI / 8.0);
    return flux * CST_1_PI / 8.0;
}

REAL8
XLALInspiralSASpinFactorizedFluxV2(REAL8Vector *polvalues,          /**< \f$(r,\phi,p_r,p_\phi)\f$ */
                                   EOBNonQCCoeffs *nqcCoeffs,       /**< pre-computed NQC coefficients */
                                   const REAL8 omega,               /**< orbital frequency */
                                   SpinEOBParams *ak,               /**< physical parameters */
                                   const REAL8 H,                   /**< real Hamiltonian */
                                   const INT4 lMax,                 /**< upper limit of the summation over l */
                                   const UINT SpinAlignedEOBversion /**< 1 for SEOBNRv1, 2 for SEOBNRv2 */
)
{
    // int	debugPK = 0;
    int i = 0;
    double radius = polvalues->data[0];
    if (radius < 1.)
    {
        return 0.;
    }
    for (i = 0; i < 4; i++)
        if (isnan(polvalues->data[i]))
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "polvalues %3.10f %3.10f %3.10f %3.10f", polvalues->data[0],
                           polvalues->data[1], polvalues->data[2], polvalues->data[3]);
            PRINT_LOG_INFO(LOG_CRITICAL, "nan polvalues:  %3.10f %3.10f %3.10f %3.10f", polvalues->data[0],
                           polvalues->data[1], polvalues->data[2], polvalues->data[3]);
            return REAL8_FAIL_NAN;
        }

    // for( i =0; i < 12; i++)
    //     if( isnan(values->data[i]) )
    //     {
    //         PRINT_LOG_INFO(LOG_CRITICAL, "values %3.10f %3.10f %3.10f %3.10f
    //         %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f",
    //             values->data[0], values->data[1], values->data[2],
    //             values->data[3], values->data[4], values->data[5],
    //             values->data[6], values->data[7], values->data[8],
    //             values->data[9], values->data[10], values->data[11]);
    //         PRINT_LOG_INFO(LOG_CRITICAL, "nan  in input values:  %3.10f %3.10f
    //         %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f
    //         %3.10f",
    //             values->data[0], values->data[1], values->data[2],
    //             values->data[3], values->data[4], values->data[5],
    //             values->data[6], values->data[7], values->data[8],
    //             values->data[9], values->data[10], values->data[11] );
    //         return REAL8_FAIL_NAN;
    //     }
    // }

    REAL8 flux = 0.0;
    REAL8 v;
    REAL8 omegaSq;
    COMPLEX16 hLM;
    INT4 l, m;

    // EOBNonQCCoeffs nqcCoeffs;
    //  if (lMax < 2)
    //  {
    //      XLAL_ERROR_REAL8(XLAL_EINVAL);
    //  }
    /* Omega is the derivative of phi */
    omegaSq = omega * omega;

    v = cbrt(omega);

    // XLAL_PRINT_INFO("v = %.16e\n", v);
    COMPLEX16 hLMTab[lMax + 1][lMax + 1];
    if (XLALSimIMRSpinEOBFluxGetSASpinFactorizedWaveform(&hLMTab[0][0], polvalues, v, H, lMax, ak) == CEV_FAILURE)
    {
        return REAL8_FAIL_NAN;
    }
    for (l = 2; l <= lMax; l++)
    {
        for (m = 1; m <= l; m++)
        {

            // if (debugPK)
            //     XLAL_PRINT_INFO("\nGetting (%d, %d) mode for flux!\n", l, m);
            // XLAL_PRINT_INFO("Stas, computing the waveform l = %d, m =%d\n", l,
            // m);
            hLM = hLMTab[l][m];
            // XLAL_PRINT_INFO("Stas: done\n");
            /*
             * For the 2,2 mode, we apply NQC correction to the
             * flux
             */
            // if (debugPK)
            //     XLAL_PRINT_INFO("\tl = %d, m = %d, mag(hLM) = %.17e, omega =
            //     %.16e\n", l, m, sqrt(creal(hLM) * creal(hLM) + cimag(hLM) *
            //     cimag(hLM)), omega);

            /* Eq. 13 */
            flux += (REAL8)(m * m) * omegaSq * (creal(hLM) * creal(hLM) + cimag(hLM) * cimag(hLM));
        }
    }
    // if( (omegaSq > 1 || flux > 5) )
    // {
    //     if(debugPK) {
    //         XLAL_PRINT_INFO("In XLALInspiralPrecSpinFactorizedFlux: omegaSq =
    //         %3.12f, FLUX = %3.12f, r = %3.12f\n",
    //                 omegaSq, flux,radius);
    //     }
    //     flux = 0.;
    // }

    // if (debugPK)
    //     XLAL_PRINT_INFO("\tStas, FLUX = %.16e\n", flux * LAL_1_PI / 8.0);
    return flux * CST_1_PI / 8.0;
}
