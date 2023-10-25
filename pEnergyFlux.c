/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "pEnergyFlux.h"
#include "pFactorizedWaveform.h"
#include "newFactorizedWaveform.h"
#include "pNQCcorrection.h"
#include "myLog.h"

REAL8
XLALInspiralSpinFactorizedFlux_SA (
                REAL8Vector *polvalues,  /**< polar dynamical variables */
				EOBNonQCCoeffs * nqcCoeffs,
							/**< pre-computed NQC coefficients */
				const REAL8 omega,	/**< orbital frequency */
                const REAL8 dr,
                const REAL8 ncrv,
				SpinEOBParams * ak,	/**< physical parameters */
				const REAL8 H,		/**< real Hamiltonian */
				const INT lMax	/**< upper limit of the summation over l */
)
{
    REAL8 flux = 0.0;
    REAL8 v;
    REAL8 omegaSq;
    COMPLEX16 hLM, hNQC;
    INT l, m;

    /* Omegs is the derivative of phi */
    omegaSq = omega * omega;

    v = cbrt (omega);

    //  printf( "v = %.16e\n", v );
    for (l = 2; l <= lMax; l++)
    {
        for (m = 1; m <= l; m++)
        {
            // if (CODE_VERSION == 3 && (l==2 || (l==3 && m==3) || (l==4 && m==4)))
            // {
            //     if (XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveformV2(
            //             &hLM, polvalues, values, v, dr, ncrv, H, l, m, ak) != CEV_SUCCESS) 
            //     {
            //         return REAL8_FAIL_NAN;
            //     }
            // }
            // else
            // {
                if (XLALSimIMRSpinEOBGetSpinFactorizedWaveform
                    (&hLM, polvalues, v, H, l, m, ak, 1) != CEV_SUCCESS)
                {
                    return REAL8_FAIL_NAN;
                }
            // }

            /* For the 2,2 mode, we apply NQC correction to the flux */
            // if (l == 2 && m == 2)
            // {
            //     XLALSimIMREOBNonQCCorrection (&hNQC, values, omega, nqcCoeffs);
            //     /* Eq. 16 */
            //     hLM *= hNQC;
            // }
            //printf( "l = %d, m = %d, mag(hLM) = %.17e, omega = %.16e\n", l, m, sqrt(creal(hLM)*creal(hLM)+cimag(hLM)*cimag(hLM)), omega );
            /* Eq. 13 */
            flux +=
                (REAL8) (m * m) * omegaSq * (creal (hLM) * creal (hLM) +
                            cimag (hLM) * cimag (hLM));
        }
    }
    return flux * CST_1_PI / 8.0;
}


REAL8
XLALInspiralSpinFactorizedFlux (
                REAL8Vector *polvalues,  /**< polar dynamical variables */
                REAL8Vector * values,	/**< cart dynamical variables */
				EOBNonQCCoeffs * nqcCoeffs,
							/**< pre-computed NQC coefficients */
				const REAL8 omega,	/**< orbital frequency */
                const REAL8 dr,
                const REAL8 ncrv,
				SpinEOBParams * ak,	/**< physical parameters */
				const REAL8 H,		/**< real Hamiltonian */
				const INT lMax	/**< upper limit of the summation over l */
)
{
    REAL8 flux = 0.0;
    REAL8 v;
    REAL8 omegaSq;
    COMPLEX16 hLM, hNQC;
    INT l, m;

    /* Omegs is the derivative of phi */
    omegaSq = omega * omega;

    v = cbrt (omega);

    //  printf( "v = %.16e\n", v );
    for (l = 2; l <= lMax; l++)
    {
        for (m = 1; m <= l; m++)
        {
            // if (CODE_VERSION == 3 && (l==2 || (l==3 && m==3) || (l==4 && m==4)))
            // {
            //     if (XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveformV2(
            //             &hLM, polvalues, values, v, dr, ncrv, H, l, m, ak) != CEV_SUCCESS) 
            //     {
            //         return REAL8_FAIL_NAN;
            //     }
            // }
            // else
            // {
                if (XLALSimIMRSpinEOBGetSpinFactorizedWaveform
                    (&hLM, polvalues, v, H, l, m, ak, 1) != CEV_SUCCESS)
                {
                    return REAL8_FAIL_NAN;
                }
            // }

            /* For the 2,2 mode, we apply NQC correction to the flux */
            // if (l == 2 && m == 2)
            // {
            //     XLALSimIMREOBNonQCCorrection (&hNQC, values, omega, nqcCoeffs);
            //     /* Eq. 16 */
            //     hLM *= hNQC;
            // }
            //printf( "l = %d, m = %d, mag(hLM) = %.17e, omega = %.16e\n", l, m, sqrt(creal(hLM)*creal(hLM)+cimag(hLM)*cimag(hLM)), omega );
            /* Eq. 13 */
            flux +=
                (REAL8) (m * m) * omegaSq * (creal (hLM) * creal (hLM) +
                            cimag (hLM) * cimag (hLM));
        }
    }
    return flux * CST_1_PI / 8.0;
}


REAL8
XLALInspiralPrecSpinFactorizedFlux(
				   REAL8Vector * polvalues,	/**< \f$(r,\phi,p_r,p_\phi)\f$ */
				   REAL8Vector * values,	/**< dynamical variables */
				   EOBNonQCCoeffs * nqcCoeffs,	/**< pre-computed NQC coefficients */
				   const REAL8 omega,	/**< orbital frequency */
                    const REAL8 dr,
                    const REAL8 ncrv,
				   SpinEOBParams * ak,	/**< physical parameters */
				   const REAL8 H,	/**< real Hamiltonian */
				   const INT4 lMax,	/**< upper limit of the summation over l */
				   const UINT SpinAlignedEOBversion	/**< 1 for SEOBNRv1, 2 for SEOBNRv2 */
)
{
    // print_log("omega, dr, ncrv, H = %.16e, %.16e, %.16e, %.16e\n",
    //     omega, dr, ncrv, H);
    // int	debugPK = 0;
    int i = 0;
    double radius = sqrt(values->data[0]*values->data[0] + values->data[1] *values->data[1]  + values->data[2] *values->data[2]  );
    // print_debug("radius = %.16e\n", radius);
    if (radius < 1.) 
    {
        return 0.;
    }
    if (1){
    for( i =0; i < 4; i++)
        if( isnan(polvalues->data[i]) ) 
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "polvalues %3.10f %3.10f %3.10f %3.10f", polvalues->data[0], polvalues->data[1], polvalues->data[2], polvalues->data[3]);
            PRINT_LOG_INFO(LOG_CRITICAL, "nan polvalues:  %3.10f %3.10f %3.10f %3.10f",polvalues->data[0], polvalues->data[1], polvalues->data[2], polvalues->data[3] );
            return REAL8_FAIL_NAN;
        }
    for( i =0; i < 12; i++)
        if( isnan(values->data[i]) ) 
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "values %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f", 
                values->data[0], values->data[1], values->data[2], values->data[3], values->data[4], values->data[5], values->data[6], values->data[7], values->data[8], values->data[9], values->data[10], values->data[11]);
            PRINT_LOG_INFO(LOG_CRITICAL, "nan  in input values:  %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f",  
                values->data[0], values->data[1], values->data[2], values->data[3], values->data[4], values->data[5], values->data[6], values->data[7], values->data[8], values->data[9], values->data[10], values->data[11] );
            return REAL8_FAIL_NAN;
        }
    }

    REAL8		flux = 0.0;
    REAL8		v;
    REAL8		omegaSq;
    COMPLEX16	hLM;
    INT4		l        , m;

    //EOBNonQCCoeffs nqcCoeffs;
    // if (lMax < 2) 
    // {
    //     XLAL_ERROR_REAL8(XLAL_EINVAL);
    // }
    /* Omega is the derivative of phi */
    omegaSq = omega * omega;

    v = cbrt(omega);
    // print_log("flux::v = %.16e\n", v);
    // return 0.0;
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
    //XLAL_PRINT_INFO("v = %.16e\n", v);
    COMPLEX16 hLMTab[lMax+1][lMax+1];
    if (XLALSimIMRSpinEOBFluxGetPrecSpinFactorizedWaveform(&hLMTab[0][0], polvalues, values, v, H,
                lMax, ak) == CEV_FAILURE) {
        return REAL8_FAIL_NAN;
    }
    for (l = 2; l <= lMax; l++) 
    {
        for (m = 1; m <= l; m++) 
        {

            // if (debugPK)
            //     XLAL_PRINT_INFO("\nGetting (%d, %d) mode for flux!\n", l, m);
            //XLAL_PRINT_INFO("Stas, computing the waveform l = %d, m =%d\n", l, m);
            // if (CODE_VERSION == 3 && (l==2 || (l==3 && m==3) || (l==4 && m==4)))
            // {
            //     if (XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveformV2(
            //             &hLM, polvalues, values, v, dr, ncrv, H, l, m, ak) != CEV_SUCCESS) 
            //     {
            //         return REAL8_FAIL_NAN;
            //     }
            // }
            // else
                hLM = hLMTab[l][m];
            //XLAL_PRINT_INFO("Stas: done\n");
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
            //     XLAL_PRINT_INFO("\tl = %d, m = %d, mag(hLM) = %.17e, omega = %.16e\n", l, m, sqrt(creal(hLM) * creal(hLM) + cimag(hLM) * cimag(hLM)), omega);

            /* Eq. 13 */
            flux += (REAL8) (m * m) * omegaSq * (creal(hLM) * creal(hLM) + cimag(hLM) * cimag(hLM));
        }
    }
    // if( (omegaSq > 1 || flux > 5) ) 
    // {
    //     if(debugPK) {
    //         XLAL_PRINT_INFO("In XLALInspiralPrecSpinFactorizedFlux: omegaSq = %3.12f, FLUX = %3.12f, r = %3.12f\n",
    //                 omegaSq, flux,radius);
    //     }
    //     flux = 0.;
    // }

    // if (debugPK)
    //     XLAL_PRINT_INFO("\tStas, FLUX = %.16e\n", flux * LAL_1_PI / 8.0);
    return flux * CST_1_PI / 8.0;
}

INT CalculateRRForceSpinCoeffs(RRForceCoeffs *coeffsFf,
                               RRForceCoeffs *coeffsFr,
                               REAL8 m1, REAL8 m2,
                               REAL8 chi1, REAL8 chi2)
{
    REAL8 Mtot = m1 + m2;
    REAL8 dm = (m1 - m2) / Mtot;
    REAL8 X1 = m1 / Mtot;
    REAL8 X1sq = X1*X1;
    REAL8 X2 = m2 / Mtot;
    REAL8 X2sq = X2*X2;
    REAL8 eta = m1*m2 / Mtot / Mtot;
    REAL8 eta2 = eta*eta;
    REAL8 chi1sq = chi1*chi2, chi2sq = chi2*chi2;
    REAL8 C1ESsq = 1., C2ESsq = 1.;
    coeffsFf->SOPre = 2.*eta2/15.;
    coeffsFf->SOf200 = chi1 * (22.*eta-9.*dm-9.) + chi2 * (22.*eta+9.*dm-9.);
    coeffsFf->SOf020 = chi1*5.*(3.+3.*dm+16.*eta) + chi2*5.*(3.-3.*dm+16.*eta);
    coeffsFf->SOf011 = chi1*(66.*eta-23.*dm-23.) + chi2*(66.*eta+23.*dm-23.);
    coeffsFf->SOf101 = chi1*(179.+179.*dm-146.*eta) + chi2*(179.-179.*dm-146.*eta);
    coeffsFf->SOf110 = chi1*6.*(9.+9.*dm-22.*eta) + chi2*6.*(9.-9.*dm-22.*eta);

    coeffsFr->SOPre = 2.*eta2/15.;
    coeffsFr->SOf001 = chi1*(179.*dm-146.*eta+179.) + chi2*(-179.*dm-146.*eta+179.);
    coeffsFr->SOf100 = chi1*(22.*eta-9.*dm-9.) + chi2*(22.*eta+9.*dm-9.);
    coeffsFr->SOf010 = -chi1*5.*(3.+3.*dm+16.*eta) -5.*chi2*(3.-3.*dm+16.*eta);

    coeffsFf->SSPre = eta2/30.;
    coeffsFf->SSf100 = chi1sq*(24.*X1sq*C1ESsq*15. + 400.*eta*dm - 290.*dm - 400.*eta2 + 818.*eta - 209.) + 
        2.*eta*chi1*chi2*(378.-400.*eta) + chi2sq*(24.*X2sq*C2ESsq*15. - 400.*eta*dm + 290.*dm - 400.*eta2 + 818.*eta - 209.);
    coeffsFf->SSf010 = chi1sq*(24.*X1sq*C1ESsq*(-90.) + 361.*dm-632.*eta*dm+632.*eta2-1354.*eta+361.) + 
        2.*eta*chi1*chi2*(362.*eta-2250.) + chi2sq*(24.*X2sq*C2ESsq*(-90.) - 361.*dm+632.*eta*dm+632.*eta2-1354.*eta+361.);
    coeffsFf->SSf001 = chi1sq*(24.*X1sq*C1ESsq*(-49.) + 355.*dm-704.*eta*dm+704.*eta2-1414.*eta+355.) + 
        2.*eta*chi1*chi2*(704.*eta-1182.) + chi2sq*(24.*X2sq*C2ESsq*(-49.) - 355.*dm+704.*eta*dm+704.*eta2-1414.*eta+355.);
    return CEV_SUCCESS;
}

INT CalculateRRForceCoeffs(RRForceCoeffs *coeffsFf,
                           RRForceCoeffs *coeffsFr,
                           SpinEOBParams *params)
{
    if (!params || !coeffsFr || !coeffsFf)
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Core params is NULL");
        return CEV_FAILURE;
    }
    memset (coeffsFf, 0, sizeof (RRForceCoeffs));
    memset (coeffsFr, 0, sizeof (RRForceCoeffs));
    REAL8 eta = params->eta;
    REAL8 eta2 = eta*eta;
    // Ff
    coeffsFf->LOPre = 8 * eta2 / 15;
    coeffsFf->LOf100 = 10;
    coeffsFf->LOf010 = -39;
    coeffsFf->LOf001 = -22;

    coeffsFf->P1Pre = eta2/105.;
    coeffsFf->P1f110 = 18.*(4.*eta+93.);
    coeffsFf->P1f200 = -6.*(4.*eta+93.);
    coeffsFf->P1f020 = 180.*(7.*eta-5.);
    coeffsFf->P1f011 = (9780.*eta + 19198.);
    coeffsFf->P1f101 = -(484.*eta+3833.);
    coeffsFf->P1f002 = 1684.*eta + 6213.;

    coeffsFf->TPre = CST_PI * eta2;
    coeffsFf->Tf201 = 334./15.;
    coeffsFf->Tf102 = -718./15.;
    coeffsFf->Tf111 = -334./5.;
    coeffsFf->Tf012 = -308./15.;
    coeffsFf->Tf030 = 49./225.;

    coeffsFf->P2Pre = eta2;
    coeffsFf->P2f300 = 5276.*eta2/315. + 109609.*eta/630. - 9175./378.;
    coeffsFf->P2f201 = 1519./54. - 9964.*eta2/315. - 100847.*eta/630.;
    coeffsFf->P2f102 = 3512.*eta2/315. - 4234.*eta/45. + 3355./21.;
    coeffsFf->P2f210 = 9175./126. - 5276.*eta2/105. - 109609.*eta/210.;
    coeffsFf->P2f111 = 104296./945. - 15121.*eta2/105. - 254732.*eta/315.;
    coeffsFf->P2f003 = 52./3. - 152.*eta2/9. - 310.*eta/9.;
    coeffsFf->P2f120 = -185.*eta2/63. + 2171.*eta/63. - 269./63.;
    coeffsFf->P2f021 = 4112./35. - 278.*eta2/45. + 344.*eta/35.;
    coeffsFf->P2f012 = -1604.*eta2/315. - 2054.*eta/7. - 8803./21.;
    coeffsFf->P2f003 = 8.*eta2/15. + 9728.*eta/315. - 190244./2835.;
    
    // Fr
    coeffsFr->LOPre = -16.*eta2/15.;
    coeffsFr->LOf100 = -5.;
    coeffsFr->LOf010 = 12.;
    coeffsFr->LOf001 = 11.;

    coeffsFr->P1Pre = eta2/105;
    coeffsFr->P1f020 = 180.*(7.*eta - 5.);
    coeffsFr->P1f020 = -6.*(4.*eta + 93.);
    coeffsFr->P1f011 = 4.*(691.*eta + 3958.);
    coeffsFr->P1f110 = -198.*(6.*eta - 13.);
    coeffsFr->P1f101 = -(484.*eta + 3833.);
    coeffsFr->P1f002 = 1684.*eta + 6213.;

    coeffsFr->TPre = CST_PI * eta2;
    coeffsFr->Tf101 = 334./15.;
    coeffsFr->Tf020 = -1./18.;
    coeffsFr->Tf011 = -7034./45.;
    coeffsFr->Tf002 = -718./15.;

    coeffsFr->P2Pre = eta2;
    coeffsFr->P2f300 = 5276.*eta2/315. + 109609.*eta/630. - 9175./378.;
    coeffsFr->P2f201 = 1519./54. - 9964.*eta2/315. - 100847.*eta/630.;
    coeffsFr->P2f102 = 3512.*eta2/315. - 4234.*eta/45. + 3355./21;
    coeffsFr->P2f210 = 9713./126. - 2129.*eta2/45. - 350537*eta/630.;
    coeffsFr->P2f111 = 26459./945. - 1745.*eta2/63. - 449227.*eta/315.;
    coeffsFr->P2f030 = 52./3. - 152.*eta2/9. - 310.*eta/9.;
    coeffsFr->P2f120 = 293.*eta2/21. + 1447.*eta/21. - 1361./63.;
    coeffsFr->P2f021 = 41612./315. - 1894.*eta2/45. + 28648.*eta/315.;
    coeffsFr->P2f012 = 34528.*eta2/945. - 481624.*eta/945. - 18761./45.;
    coeffsFr->P2f003 = 8.*eta2/15. + 9728.*eta/315. - 190244./2835.;

    CalculateRRForceSpinCoeffs(coeffsFf, coeffsFr, params->m1, params->m2, params->chi1, params->chi2);
    return CEV_SUCCESS;
}


static INT CalculateRRComponents(RRForceCoeffs *coeffs, 
                          REAL8 r, REAL8 pr2, REAL8 p2,
                          REAL8 *tLO, REAL8 *tPN1, REAL8 *tPN2, REAL8 *tTail,
                          REAL8 *tSO, REAL8 *tSS)
{
    REAL8 t100, t010, t001;
    REAL8 t200, t020, t002, t110, t101, t011;
    REAL8 t300, t030, t003, t210, t201, t120, t021, t102, t012, t111;
    t100 = p2;
    t010 = pr2;
    t001 = 1./r;
    t200 = t100 * t100;
    t020 = t010 * t010;
    t002 = t001 * t001;
    t110 = t100 * t010;
    t011 = t010 * t001;
    t101 = t100 * t001;
    t300 = t200 * t100;
    t030 = t020 * t010;
    t003 = t002 * t001;
    t210 = t200 * t010;
    t201 = t200 * t001;
    t120 = t100 * t020;
    t021 = t020 * t001;
    t102 = t100 * t002;
    t012 = t010 * t002;
    t111 = t110 * t001;
    REAL8 rr_LO, rr_P1, rr_Tail, rr_P2, rr_SO, rr_SS;
    rr_LO = coeffs->LOPre * (coeffs->LOf001 * t001 + coeffs->LOf010 * t010 + coeffs->LOf001 * t001);
    rr_P1 = coeffs->P1Pre * (coeffs->P1f002 * t002 + coeffs->P1f020 * t020 + coeffs->P1f002 * t002 + 
                coeffs->P1f110 * t110 + coeffs->P1f011 * t011 + coeffs->P1f101 * t101);
    rr_Tail = coeffs->TPre * (coeffs->Tf201 * t201 + coeffs->Tf102 * t102 + coeffs->Tf111 * t111 + 
                coeffs->Tf012 * t012 + coeffs->Tf030 * t030 + coeffs->Tf101 * t101 + coeffs->Tf020 * t020 + 
                coeffs->Tf011 * t011 + coeffs->Tf002 * t002);
    rr_P2 = coeffs->P2Pre * (coeffs->P2f300 * t300 + coeffs->P2f030 * t030 + coeffs->P2f003 * t003 + 
                coeffs->P2f210 * t210 + coeffs->P2f201 * t201 + coeffs->P2f120 * t120 + coeffs->P2f021 * t021 + 
                coeffs->P2f102 * t102 + coeffs->P2f012 * t012 + coeffs->P2f111 * t111);
    rr_SS = coeffs->SSPre * (coeffs->SSf001 * t001 + coeffs->SSf010 * t010 + coeffs->SSf100 * t100);
    rr_SO = coeffs->SOPre * (coeffs->SOf001 * t001 + coeffs->SOf010 * t010 + coeffs->SOf100 * t100 + 
        coeffs->SOf020 * t020 + coeffs->SOf200 * t200 + coeffs->SOf011 * t011 + 
        coeffs->SOf101 * t101 + coeffs->SOf110 * t110);
    *tLO = rr_LO;
    *tTail = rr_Tail;
    *tPN1 = rr_P1;
    *tPN2 = rr_P2;
    *tSO = rr_SO;
    *tSS = rr_SS;
    return CEV_SUCCESS;
}

INT CalculateRRForce(RRForceCoeffs *coeffsFf,
                       RRForceCoeffs *coeffsFr,
                       REAL8 *rrFf, REAL8 *rrFr,
                       REAL8 r, REAL8 pr, REAL8 p2, REAL8 pf)
{
    REAL8 r2 = r*r;
    REAL8 r3 = r2*r;
    REAL8 pr2 = pr*pr;
    REAL8 Ff_LO, Ff_P1, Ff_Tail, Ff_P2, Ff_SO, Ff_SS;
    CalculateRRComponents(coeffsFf, r, pr2, p2, &Ff_LO, &Ff_P1, &Ff_P2, &Ff_Tail, &Ff_SO, &Ff_SS);
    Ff_LO *= pf / r3;
    Ff_P1 *= pf / r3;
    Ff_Tail *= 1./r2;
    Ff_P2 *= pf / r3;
    Ff_SO *= 1./r3/r2;
    Ff_SS *= pf / r3 / r2;
    REAL8 Fr_LO, Fr_P1, Fr_Tail, Fr_P2, Fr_SO, Fr_SS;
    CalculateRRComponents(coeffsFr, r, pr2, p2, &Fr_LO, &Fr_P1, &Fr_P2, &Fr_Tail, &Fr_SO, &Fr_SS);
    Fr_LO *= pr / r3;
    Fr_P1 *= pr / r3;
    Fr_Tail *= pf*pr / r3 / r;
    Fr_P2 *= pr / r3;
    Fr_SO *= pr * pf /r3 / r2;
    Fr_SS *= pr / r3 / r2;
    *rrFf = Ff_LO + Ff_P1 + Ff_Tail + Ff_P2 + Ff_SO + Ff_SS;
    *rrFr = Fr_LO + Fr_P1 + Fr_Tail + Fr_P2 + Fr_SO + Fr_SS;
    return CEV_SUCCESS;
}


INT CalculateEccCorrectionToFlux(REAL8 eta, REAL8 chi1, REAL8 chi2,
                                     REAL8 r, REAL8 vr, REAL8 prDot,
                                     REAL8 *cFr, REAL8 *cFf)
{
    REAL8 h0, g0, g0_2, g0_3, g0_4;
    REAL8 h1, f1, f1_2, f1_3;
    REAL8 h2, f2, f2_2, f2_3;
    REAL8 sqh0, sqh, h0ph2sq;
    REAL8 chi1_2, chi2_2, eta_2;
    h0 = 1./r;
    h1 = vr*vr;
    h2 = r*prDot;
    g0 = sqrt(h0);
    g0_2 = h0;
    g0_3 = g0_2*g0;
    g0_4 = g0_3*g0;
    f1 = h1 / h0;
    f1_2 = f1*f1;
    f1_3 = f1_2*f1;
    f2 = h2 / h0;
    f2_2 = f2*f2;
    f2_3 = f2_2*f2;
    chi1_2 = chi1*chi1;
    chi2_2 = chi2*chi2;
    eta_2 = eta*eta;
    sqh = sqrt((h0+h2)/h0);
    h0ph2sq = pow(h0+h2, 2.);
    REAL8 corrFrPN0, corrFfPN0;
    REAL8 BaseFr = 1., BaseFf = 1.;
    corrFrPN0 = ((6*h0 + 7*h1 - 5*h2)*(h0 + h2))/(6.*h0*h0);
    corrFfPN0 = ((12*h0 + 29*h1 - 10*h2)*(h0 + h2))/(12.*h0*h0);
    BaseFf = (12. + 29* f1 - 10* f2) *(1 + f2);
    BaseFr = (6. + 7 *f1 - 5* f2)* (1 + f2);
    REAL8 dm;
    dm = sqrt(1.-4.*eta);
    REAL8    sqf2, sqf2_2, sqf2_4;
    sqf2 = sqrt(GET_MAX(1.+f2, 0.0));
    sqf2_2 = sqf2*sqf2;
    sqf2_4 = sqf2_2 * sqf2_2;
    REAL8 corrFfPN1 = 0., corrFfPN32 = 0., corrFfPN2 = 0.;
    REAL8 corrFrPN1 = 0., corrFrPN32 = 0., corrFrPN2 = 0.;
    REAL8 chi12_2;
    chi12_2 = (chi1+chi2) * (chi1+chi2);
    // NLO terms
    // PN1
        corrFfPN1 = (24*f1_2*(149 - 460*f2 + eta*(-124 + 79*f2)) + 4*f2*(-62 + f2*(2203 + 417*f2) - 12*eta*(105 + f2*(-17 + 32*f2))) + 
            f1*(14713 + f2*(5317 + 11436*f2) + 12*eta*(-4219 + f2*(-4063 + 660*f2))))/336./BaseFf;
        corrFrPN1 = (6*f1_2*(-181 - 475*f2 + 2*eta*(43 + 92*f2)) + f1*(-8009 + f2*(-11708 + 249*f2) + 18*eta*(-66 + 5*f2*(13 + 29*f2))) + 
            f2*(-62 + f2*(2203 + 417*f2) - 12*eta*(105 + f2*(-17 + 32*f2))))/168./BaseFr;
    // PN1.5
        corrFfPN32 = (-480.*(22. + 29*f1)*(1 + f2)*CST_PI + 30.*sqf2_4*CST_PI*(160 - 167*sqf2) + (10770*(1 + f2) + f1*(20400 + 10020*f1 - 49*f1_2 + 5010*f2))*CST_PI*sqf2 + 
            5*chi1*((73*(1 + dm) - 56*eta)*(22 + 29*f1)*(1 + f2) - 6*((1 + dm)*(155 - 21*f2 + 4*f1*(106 + 15*f1 + 9*f2)) - 2*eta*(67 + 23*f2 + f1*(142 + 15*f1 + 44*f2)))*
            sqf2 + 2*sqf2_4*(eta*(280 - 66*sqf2) + (1 + dm)*(-365 + 27*sqf2))) + 
            5*chi2*((73 - 73*dm - 56*eta)*(22 + 29*f1)*(1 + f2) + 6*((-1 + dm)*(155 - 21*f2 + 4*f1*(106 + 15*f1 + 9*f2)) + 2*eta*(67 + 23*f2 + f1*(142 + 15*f1 + 44*f2)))*
            sqf2 - 2*sqf2_4*((-1 + dm)*(-365 + 27*sqf2) + eta*(-280 + 66*sqf2))))/120./BaseFf;
        corrFrPN32 = (-384*(11 + 7*f1)*(1 + f2)*CST_PI + 12*sqf2_4*CST_PI*(160 - 167*sqf2) + (4308 + f1*(12064 + 5*f1))*(1 + f2)*CST_PI*sqf2 + 
            4*chi1*((73*(1 + dm) - 56*eta)*(11 + 7*f1)*(1 + f2) + 3*((1 + dm)*(-155 - 88*f1 + 3*(7 + 8*f1)*f2) + 2*eta*(67 + 23*f2 + f1*(57 + 29*f2)))*sqf2 + 
            sqf2_4*(eta*(280 - 66*sqf2) + (1 + dm)*(-365 + 27*sqf2))) - 
            4*chi2*((73*(-1 + dm) + 56*eta)*(11 + 7*f1)*(1 + f2) - 3*((1 - dm)*(-155 - 88*f1 + 3*(7 + 8*f1)*f2) + 2*eta*(67 + 23*f2 + f1*(57 + 29*f2)))*sqf2 + 
            sqf2_4*((-1 + dm)*(-365 + 27*sqf2) + eta*(-280 + 66*sqf2)))) / 96./BaseFr;
    // PN2
        /*-2688.*(-100731. + 2*eta*(15399 + 707*eta))*h0_4 + */
        corrFfPN2 = (2016*f1_3*(-55655 - 43604*f2 + 3*eta*(108227 + 102218*f2 + 3*eta*(4246 + 4607*f2))) + 
            72*f1_2*(-9410885 + 2*(-6875383 - 102312*chi2_2*(-1 + dm) + 102312*chi1_2*(1 + dm) - 1125495*f2)*f2 + 
            84*eta_2*(62602 + f2*(63458 + 4872*chi12_2 + 15475*f2)) + 
            3*eta*(11335354 + f2*(18428251 + 136416*chi2_2*(-2 + dm) - 136416*chi1_2*(2 + dm) + 4403546*f2))) + 
            f1*(-676564289 - 1091884877*f2 - 127008*chi1*chi2*eta*(-509 + 187*f2 + 8*eta*(9 + 2*f2*(-51 + 10*f2))) - 
            31752*chi2_2*(16*eta_2*(9 + 2*f2*(-51 + 10*f2)) - (-1 + dm)*(-235 + f2*(-427 + 160*f2)) + 2*eta*(163 + 72*dm + f2*(1243 - 816*dm + 160*(-2 + dm)*f2))) + 
            31752*chi1_2*(-16*eta_2*(9 + 2*f2*(-51 + 10*f2)) - (1 + dm)*(-235 + f2*(-427 + 160*f2)) + 2*eta*(-163 + 72*dm + f2*(-1243 - 816*dm + 160*(2 + dm)*f2))) + 
            36*(3*f2_2*(-2684653 + 119616*f2) + 84*eta_2*(92743 + f2*(224507 + 16*f2*(5218 + 195*f2))) + 2*eta*(7460071 + f2*(40269487 + 6*f2*(4360423 + 44856*f2)))))
            + 4*f2*(-4371470 + 41899357*f2 - 127008*chi1*chi2*eta*(49 - 83*f2 + 12*eta*(-1 + 10*f2)) - 
            31752*chi2_2*(-23*(-1 + dm)*(1 + f2) + 24*eta_2*(-1 + 10*f2) + eta*(-22 - 24*dm - 286*f2 + 240*dm*f2)) + 
            31752*chi1_2*(eta_2*(24 - 240*f2) - 23*(1 + dm)*(1 + f2) + eta*(22 - 24*dm + 286*f2 + 240*dm*f2)) + 
            9*(f2_2*(6574215 + 1177652*f2) - 336*eta_2*(399 + f2*(2569 + f2*(5669 + 2707*f2))) - 4*eta*(287308 + f2*(5572558 + 3*f2*(2322597 + 728119*f2))))))/1016064./BaseFf;

        corrFrPN2 = (63504*chi1*chi2*eta*(-352 + f1*(187 + 19*f2) + 2*f2*(-225 + 83*f2) + 8*eta*(44 + 3*f1 + (47 + 2*f1*(18 + 7*f1))*f2 - 10*(3 + f1)*f2_2)) + 
            504*f1_3*(-58295 - 38819*f2 + 3*eta*(107017 + 99973*f2 + 2*eta*(5609 + 6308*f2))) + 
            54*f1_2*(-4422473 - f2*(4678771 + 32928*chi2_2*(-1 + dm) + 750890*f2) + 28*eta_2*(40630 + f2*(53345 + 2352*chi2_2 + 17755*f2)) + 
            2*eta*(9068956 + f2*(11658445 + 32928*chi2_2*(-2 + dm) + 2212693*f2))) + 
            f2*(-4371470 + f2*(41899357 + 9*f2*(6574215 + 1177652*f2)) - 
            31752*chi2_2*(-23*(-1 + dm)*(1 + f2) + 24*eta_2*(-1 + 10*f2) + eta*(-22 - 24*dm - 286*f2 + 240*dm*f2)) - 
            3024*eta_2*(399 + f2*(2569 + f2*(5669 + 2707*f2))) - 36*eta*(287308 + f2*(5572558 + 3*f2*(2322597 + 728119*f2)))) + 
            f1*(-378748211 - f2*(511759760 + 9*f2*(15676013 + 1023792*f2)) + 504*eta_2*(36806 + f2*(120878 + 3*f2*(27789 + 5045*f2))) + 
            18*eta*(58491518 + f2*(113256329 + 3*f2*(16610749 + 359632*f2))) - 
            15876*chi2_2*(16*eta_2*(-3 + 2*f2*(-18 + 5*f2)) - (-1 + dm)*(-125 + f2*(-221 + 80*f2)) + 2*eta*(149 - 24*dm + f2*(509 - 288*dm + 80*(-2 + dm)*f2)))) - 
            15876*chi1_2*(112*(-1 - dm + 2*(2 + dm - eta)*eta)*f1_2*f2 + 2*f2*(23*(1 + dm)*(1 + f2) + 24*eta_2*(-1 + 10*f2) - 2*eta*(11 - 12*dm + 143*f2 + 120*dm*f2)) + 
            f1*(16*eta_2*(-3 + 2*f2*(-18 + 5*f2)) + (1 + dm)*(-125 + f2*(-221 + 80*f2)) + 2*eta*(149 + 24*dm + f2*(509 + 288*dm - 80*(2 + dm)*f2))))) / 508032. / BaseFr;
    // *cFf = corrFfPN0 * (1. + corrFfPN1 * g0_2 -+corrFfPN32 * g0_3 + (corrFfPN2)*g0_4);
    *cFf = corrFfPN0 / (1. - corrFfPN1 * g0_2 - corrFfPN32 * g0_3 + (corrFfPN1*corrFfPN1 - corrFfPN2)*g0_4);
    *cFr = corrFrPN0 / (1. - corrFrPN1 * g0_2 - corrFrPN32 * g0_3 + (corrFrPN1*corrFrPN1 - corrFrPN2)*g0_4);
    // print_debug("%e\t%e\t%e\t%e\n", corrFfPN0, corrFfPN1, 1+f2, corrFfPN32);

    return CEV_SUCCESS;
}

INT CalculateEccCorrectionToFluxV2(REAL8 eta, REAL8 chi1, REAL8 chi2,
                                     REAL8 r, REAL8 prt, REAL8 prDot,
                                     REAL8 *cFr, REAL8 *cFf)
{
    REAL8 h0, g0, g0_2, g0_3, g0_4;
    REAL8 h1, f1, f1_2, f1_3;
    REAL8 h2, f2, f2_2, f2_3, f2_4;
    REAL8 chi1_2, chi2_2, eta_2;
    h0 = 1./r;
    h1 = prt*prt;
    h2 = r*prDot;
    g0 = sqrt(h0);
    g0_2 = h0;
    g0_3 = g0_2*g0;
    g0_4 = g0_3*g0;
    f1 = h1 / h0;
    f1_2 = f1*f1;
    f1_3 = f1_2*f1;
    f2 = h2 / h0;
    f2_2 = f2*f2;
    f2_3 = f2_2*f2;
    f2_4 = f2_3*f2;
    chi1_2 = chi1*chi1;
    chi2_2 = chi2*chi2;
    eta_2 = eta*eta;
    REAL8 corrFrPN0, corrFfPN0;
    REAL8 BaseFr = 1., BaseFf = 1.;
    corrFrPN0 = (6 + 7*f1 - 5*f2)/6.;
    corrFfPN0 = 1 + (29*f1)/12. - (5*f2)/6.;
    REAL8 dm;
    dm = sqrt(1.-4.*eta);
    REAL8    sqf2, sqf2_2, sqf2_3, sqf2_4, sqf2_5, sqf2_6;
    sqf2 = sqrt(GET_MAX(1.+f2, 0.0001));
    sqf2_2 = sqf2*sqf2;
    sqf2_3 = sqf2_2 * sqf2;
    sqf2_4 = sqf2_2 * sqf2_2;
    sqf2_5 = sqf2_4 * sqf2;
    sqf2_6 = sqf2_5 * sqf2;
    REAL8 corrFfPN1 = 0., corrFfPN32 = 0., corrFfPN2 = 0.;
    REAL8 corrFrPN1 = 0., corrFrPN32 = 0., corrFrPN2 = 0.;
    REAL8 chi12_2;
    chi12_2 = (chi1+chi2) * (chi1+chi2);

        corrFrPN1 = -((6*f1_2*(377 + 671*f2 + 2*eta*(55 + 6*f2)) + f2*(62 - f2*(2203 + 417*f2) + 12*eta*(105 + f2*(-17 + 32*f2))) + 
            f1*(11537 + f2*(16412 + 927*f2) - 6*eta*(-2 + f2*(195 + 239*f2)))))/(168.*(6 + 7*f1 - 5*f2));
        corrFfPN1 = -((24*f1_2*(257 + 866*f2 + eta*(530 + 327*f2)) + 4*f2*(62 - f2*(2203 + 417*f2) + 12*eta*(105 + f2*(-17 + 32*f2))) + 
            f1*(14519 + (33659 - 1692*f2)*f2 + 12*eta*(3407 + f2*(4063 + 152*f2)))))/(336.*(12 + 29*f1 - 10*f2));

        corrFrPN32 = ((-384*(11 + 7*f1)*(1 + f2)*CST_PI + (4308 + f1*(12064 + 5*f1))*CST_PI*sqf2_3 + 12*CST_PI*(160 - 167*sqf2)*sqf2_4 - 
            4*chi1*((-73*(1 + dm) + 56*eta)*(11 + 7*f1)*(1 + f2) + 
            3*((-1 - dm)*(-155 - 88*f1 + 3*(7 + 8*f1)*f2) - 2*eta*(67 + 23*f2 + f1*(57 + 29*f2)))*sqf2 + 
            5*(73*(1 + dm) - 56*eta)*sqf2_4 - 3*(9 + 9*dm - 22*eta)*sqf2_5) + 
            4*chi2*((73 - 73*dm - 56*eta)*(11 + 7*f1)*(1 + f2) + 5*(73*(-1 + dm) + 56*eta)*sqf2_4 - 
            3*((-1 + dm)*(-155 - 88*f1 + 3*(7 + 8*f1)*f2) - 2*eta*(67 + 23*f2 + f1*(57 + 29*f2)))*sqf2 - 3*(-9 + 9*dm + 22*eta)*sqf2_5))
            )/(96.*(6 + 7*f1 - 5*f2));

        corrFfPN32 = (sqf2*(-480*(22 + 29*f1)*CST_PI*sqf2 + 30*CST_PI*(160 - 167*sqf2)*sqf2_3 + 
            CST_PI*(f1*(20400 + 10020*f1 - 49*f1_2 + 5010*f2) + 10770*sqf2_2) + 
            5*chi1*(-6*(1 + dm)*(155 - 21*f2 + 4*f1*(106 + 15*f1 + 9*f2)) + 12*eta*(67 + 23*f2 + f1*(142 + 15*f1 + 44*f2)) + 
            (73 + 73*dm - 56*eta)*(22 + 29*f1)*sqf2 + 2*sqf2_3*(eta*(280 - 66*sqf2) + (1 + dm)*(-365 + 27*sqf2))) + 
            5*chi2*(6*((-1 + dm)*(155 - 21*f2 + 4*f1*(106 + 15*f1 + 9*f2)) + 2*eta*(67 + 23*f2 + f1*(142 + 15*f1 + 44*f2))) + 
            (73 - 73*dm - 56*eta)*(22 + 29*f1)*sqf2 - 2*sqf2_3*((-1 + dm)*(-365 + 27*sqf2) + eta*(-280 + 66*sqf2)))))/
            (120.*(12 + 29*f1 - 10*f2));

        corrFrPN2 = ((504*f1_3*(-38207 + 2437*f2 + 6*eta_2*(5753 + 5276*f2) + eta*(336711 + 329691*f2)) + 
            54*f1_2*(-3609073 - 3*(992273 + 10976*chi2_2*(-1 + dm))*f2 - 478562*f2_2 + 
            28*eta_2*(40366 + (52181 + 2352*chi2_2)*f2 + 11647*f2_2) + 
            2*eta*(9162504 + 3*(3869143 + 10976*chi2_2*(-2 + dm))*f2 + 2213953*f2_2)) - 
            f2*(4371470 - 41899357*f2 - 59167935*f2_2 - 10598868*f2_3 + 3024*eta_2*(399 + 2569*f2 + 5669*f2_2 + 2707*f2_3) + 
            36*eta*(287308 + 5572558*f2 + 6967791*f2_2 + 2184357*f2_3) + 
            31752*chi2_2*(23*(1 + f2) - 23*dm*(1 + f2) + 24*dm*eta*(-1 + 10*f2) + 24*eta_2*(-1 + 10*f2) - 22*eta*(1 + 13*f2))) + 
            63504*chi1*chi2*eta*(-352 - 450*f2 + 166*f2_2 + f1*(187 + 19*f2) + 
            8*eta*(44 + 47*f2 + 14*f1_2*f2 - 30*f2_2 + f1*(31 + 64*f2 - 10*f2_2))) - 
            15876*chi1_2*(-112*(1 + dm - 4*eta - 2*dm*eta + 2*eta_2)*f1_2*f2 + 
            2*f2*(dm*eta*(24 - 240*f2) + 23*(1 + f2) + 23*dm*(1 + f2) + 24*eta_2*(-1 + 10*f2) - 22*eta*(1 + 13*f2)) + 
            f1*(-349 - 445*f2 + 80*f2_2 + eta*(1194 + 1914*f2 - 320*f2_2) + 16*eta_2*(-31 - 64*f2 + 10*f2_2) + 
            dm*(-349 - 445*f2 + 80*f2_2 - 16*eta*(-31 - 64*f2 + 10*f2_2)))) + 
            f1*(-306090563 - 382087616*f2 - 101841669*f2_2 - 3108672*f2_3 + 504*eta_2*(36734 + 131498*f2 + 84951*f2_2 + 3003*f2_3) + 
            18*eta*(58189286 + 112646993*f2 + 50618319*f2_2 + 930216*f2_3) - 
            15876*chi2_2*(-349 - 445*f2 + 80*f2_2 + eta*(1194 + 1914*f2 - 320*f2_2) + 16*eta_2*(-31 - 64*f2 + 10*f2_2) + 
            dm*(349 + 445*f2 - 80*f2_2 + 16*eta*(-31 - 64*f2 + 10*f2_2))))))/(508032.*(6 + 7*f1 - 5*f2));

        corrFfPN2 = (((-810040625 + 36927576*chi1_2*(1 + dm) + 72*(-512883*chi2_2*(-1 + dm) + 
            (15377953 - 882*(-1018*chi1*chi2 + chi2_2*(2019 - 856*dm) + chi1_2*(2019 + 856*dm)))*eta + 
            42*(51859 + 17976*chi12_2)*eta_2))*f1 + 72*(-10929983 + 12*eta*(3030017 + 526372*eta))*f1_2 - 17485880*f2 + 
            2016*(-1449*(-(chi2_2*(-1 + dm)) + chi1_2*(1 + dm)) - 
            2*(10261 + 63*(98*chi1*chi2 + chi1_2*(-11 + 12*dm) - chi2_2*(11 + 12*dm)))*eta + 126*(-19 + 12*chi12_2)*eta_2)*f2
            + (5*(-232959913 + 8604792*chi1_2*(1 + dm)) + 
            72*(-597555*chi2_2*(-1 + dm) + (49340479 - 882*(374*chi1*chi2 + chi2_2*(3099 - 1744*dm) + chi1_2*(3099 + 1744*dm)))*eta + 
            42*(221507 + 36624*chi12_2)*eta_2))*f1*f2 + 
            4*(41899357 - 730296*chi1_2*(1 + dm) + 72*(10143*chi2_2*(-1 + dm) + 
            (-2786279 + 882*(166*chi1*chi2 + chi2_2*(143 - 120*dm) + chi1_2*(143 + 120*dm)))*eta - 
            294*(367 + 360*chi12_2)*eta_2))*f2_2 + 
            36*(56*(-51767 + 341151*eta + 61758*eta_2)*f1_3 + 
            2*(16*(-731063 - 12789*chi2_2*(-1 + dm) + 12789*chi1_2*(1 + dm)) - 
            87*(-634993 - 4704*chi2_2*(-2 + dm) + 4704*chi1_2*(2 + dm))*eta + 4704*(1612 + 87*chi12_2)*eta_2)*f1_2*f2 - 
            3*(3286121 - 47040*chi2_2*(-1 + dm) + 47040*chi1_2*(1 + dm) - 
            96*(195341 - 980*chi2_2*(-2 + dm) + 980*chi1_2*(2 + dm))*eta + 560*(-6521 + 168*chi12_2)*eta_2)*f1*f2_2 - 
            3*(-2191405 + 9290388*eta + 634928*eta_2)*f2_3) + 
            1008*f2*(2*(4132 + 352356*eta + 50391*eta_2)*f1_3 + 3*(-47823 + eta*(308599 + 33854*eta))*f1_2*f2 + 
            36*(369 - 316*eta + 6*eta_2)*f1*f2_2 + (42059 - 3*eta*(104017 + 10828*eta))*f2_3)))/
            (1.016064e6*(12 + 29*f1 - 10*f2));

    *cFf = corrFfPN0 * sqf2_6 / (sqf2_4 - sqf2_2 * corrFfPN1 * g0_2 - sqf2_2 * corrFfPN32 * g0_3 + (corrFfPN1*corrFfPN1 - corrFfPN2 * sqf2_2)*g0_4);
    *cFr = corrFrPN0 * sqf2_6 / (sqf2_4 - sqf2_2 * corrFrPN1 * g0_2 - sqf2_2 * corrFrPN32 * g0_3 + (corrFrPN1*corrFrPN1 - corrFrPN2 * sqf2_2)*g0_4);
    return CEV_SUCCESS;
}

void CalculateEccCorrectionCoeffs(REAL8 eta, REAL8 chi1, REAL8 chi2, EccCorrectionCoeffs *coeffs)
{
    REAL8 pi = CST_PI;
    REAL8 dm;
    dm = sqrt(1.-4.*eta);
    REAL8 eta2 = eta*eta;
    REAL8 chi1_2 = chi1*chi1;
    REAL8 chi2_2 = chi2*chi2;
    REAL8 chi12_2 = (chi1+chi2) * (chi1+chi2);

    coeffs->FrPN1_f21 = (-181. + 86.*eta)/28.;
    coeffs->FrPN1_f20 = (-181. + 86.*eta)/28.;
    coeffs->FrPN1_f12 = -6341./336. + (226.*eta)/7.;
    coeffs->FrPN1_f11 = (-22247. + 8024.*eta)/336.;
    coeffs->FrPN1_f10 = -8009./168. - (99.*eta)/14.;
    coeffs->FrPN1_f03 = (3919. - 4188.*eta)/336.;
    coeffs->FrPN1_f02 = (2029. - 268.*eta)/336.;
    coeffs->FrPN1_f01 = -857./168. + (67.*eta)/7.;

    coeffs->FfPN1_f21 = 149./14. - (62.*eta)/7.;
    coeffs->FfPN1_f20 = 149./14. - (62.*eta)/7.;
    coeffs->FfPN1_f12 = (-11497. + 29436.*eta)/336.;
    coeffs->FfPN1_f11 = (460. - 2881.*eta)/42.;
    coeffs->FfPN1_f10 = (14713. - 50628.*eta)/336.;
    coeffs->FfPN1_f03 = 3919./168. - (349.*eta)/14.;
    coeffs->FfPN1_f02 = 2029./168. - (67.*eta)/42.;
    coeffs->FfPN1_f01 = -857./84. + (134.*eta)/7.;

    coeffs->FrPN32_f21 = (-5.*pi)/96.;
    coeffs->FrPN32_f20 = (-5.*pi)/96.;
    coeffs->FrPN32_f11 = (chi1*(-331. - 331.*dm + 78.*eta) + chi2*(-331. + 331.*dm + 78.*eta) - 3016.*pi)/24.;
    coeffs->FrPN32_f10 = (chi1*(-247. - 247.*dm + 50.*eta) + chi2*(-247. + 247.*dm + 50.*eta) - 2344.*pi)/24.;
    coeffs->FrPN32_f02 = (2.*chi1*(79. + 79.*dm - 57.*eta) - 2.*chi2*(-79. + 79.*dm + 57.*eta) + 501.*pi)/24.;
    coeffs->FrPN32_f01 = (chi2*(26. - 26.*dm - 70.*eta) + chi1*(26. + 26.*dm - 70.*eta) - 555.*pi)/24.;

    coeffs->FfPN32_f30 = (49.*pi)/120.;
    coeffs->FfPN32_f20 = (15.*chi1*(2. + 2.*dm - eta) - 15.*chi2*(-2. + 2.*dm + eta) - 167.*pi)/2.;
    coeffs->FfPN32_f11 = (chi1*(-857. - 857.*dm + 516.*eta) + chi2*(-857. + 857.*dm + 516.*eta) - 1002.*pi)/24.;
    coeffs->FfPN32_f10 = (chi2*(427. - 427.*dm - 80.*eta) + chi1*(427. + 427.*dm - 80.*eta) - 1296.*pi)/24.;
    coeffs->FfPN32_f02 = (2.*chi1*(79. + 79.*dm - 57.*eta) - 2.*chi2*(-79. + 79.*dm + 57.*eta) + 501.*pi)/12.;
    coeffs->FfPN32_f01 = (chi2*(26. - 26.*dm - 70.*eta) + chi1*(26. + 26.*dm - 70.*eta) - 555.*pi)/12.;

    coeffs->FrPN2Topfm2_f20 = (3. - 13.*eta + 4.*eta2)/6.;
    coeffs->FrPN2Topfm2_f10 = (143. - 571.*eta - 4.*eta2)/126.;
    coeffs->FrPN2Topfm2_f00 = (-11.*(-22. + 71.*eta + 68.*eta2))/441.;
    coeffs->FrPN2Topfm1_f20 = (93. - 368.*eta - 16.*eta2)/196.;
    coeffs->FrPN2Topfm1_f10 = (36983. - 125676.*eta - 89024.*eta2)/21168.;
    coeffs->FrPN2Topfm1_f00 = (-6803. + 6456.*eta + 83024.*eta2)/21168.;

    coeffs->FrPN2Topf0_f30 = (-58295. + 321051.*eta + 33654.*eta2)/1008.;
    coeffs->FrPN2Topf0_f21 = (-719923. + 4400322.*eta + 555796.*eta2)/9408.;
    coeffs->FrPN2Topf0_f20 = (-4431641. + 18175960.*eta + 1132136.*eta2)/9408.;
    coeffs->FrPN2Topf0_f12 = (1403989. - 4613140.*eta + 7636304.*eta2)/112896.;
    coeffs->FrPN2Topf0_f11 = (-17323609. + 193297692.*eta + 7028624.*eta2)/112896.;
    coeffs->FrPN2Topf0_f10 = (chi1*chi2*eta*(187. + 24.*eta))/8. + 
        (chi1_2*(125. + dm*(125. - 48.*eta) - 298.*eta + 48.*eta2))/32. + 
        (chi2_2*(125. - 125.*dm - 298.*eta + 48.*dm*eta + 48.*eta2))/32. + 
        (-380212379. + 1058165820.*eta + 20702928.*eta2)/508032.;
    coeffs->FrPN2Topf0_f03 = (936881. - 13825864.*eta - 3945616.*eta2)/112896.;
    coeffs->FrPN2Topf0_f02 = (7157819. - 40117320.*eta + 91184.*eta2)/112896.;
    coeffs->FrPN2Topf0_f01 = 16168237./508032. - (21115.*eta)/588. + 
        (23.*chi1*chi2*eta)/4. - (16889.*eta2)/2646. - (15.*chi12_2*eta2)/2. - 
        (chi2_2*(23. - 23.*dm - 166.*eta + 120.*dm*eta))/16. + (chi1_2*(-23. + 166.*eta + dm*(-23. + 120.*eta)))/16.;
    coeffs->FrPN2Topf0_f00 = -4813./21168. + (431./294. - 44.*chi1*chi2)*eta 
        + (-2945./1323. + (15.*chi1_2)/2. - (15.*chi12_2)/2. + 59.*chi1*chi2 + 
        (15.*chi2_2)/2.)*eta2;

    coeffs->FfPN2Topfm2_f23 = (-29.*(3. - 13.*eta + 4.*eta2))/42.;
    coeffs->FfPN2Topfm2_f14 = (1331. - 5062.*eta - 1048.*eta2)/441.;
    coeffs->FfPN2Topfm2_f05 = (22.*(-22. + 71.*eta + 68.*eta2))/441.;

    coeffs->FfPN2Topfm1_f22 = (627. - 2602.*eta + 376.*eta2)/98.;
    coeffs->FfPN2Topfm1_f13 = (-107093. + 327840.*eta + 402128.*eta2)/10584.;
    coeffs->FfPN2Topfm1_f04 = (51277. - 180984.*eta - 96496.*eta2)/10584.;

    coeffs->FfPN2Topf0_f30 = (-55655. + 324681.*eta + 38214.*eta2)/504.;
    coeffs->FfPN2Topf0_f21 = (-788870. + 4644559.*eta + 432720.*eta2)/4704.;
    coeffs->FfPN2Topf0_f20 = (-9410885. + 34006062.*eta + 5258568.*eta2)/14112.;
    coeffs->FfPN2Topf0_f12 = (16230547. - 51071064.*eta + 35968592.*eta2)/338688.;
    coeffs->FfPN2Topf0_f11 = (-8917985. + 100392366.*eta - 1753660.*eta2)/42336.;
    coeffs->FfPN2Topf0_f10 = -676564289./1016064. + (7460071.*eta)/14112. + 
        (509.*chi1*chi2*eta)/8. + (13249.*eta2)/48. - (9.*chi12_2*eta2)/2. - 
        (chi2_2*(-235. + 235.*dm + 326.*eta + 144.*dm*eta))/32. + (chi1_2*(235. - 
        326.*eta + dm*(235. + 144.*eta)))/32.;
    coeffs->FfPN2Topf0_f03 = (2176067. - 39181656.*eta - 10867376.*eta2)/169344.;
    coeffs->FfPN2Topf0_f02 = (21922177. - 122048088.*eta - 121456.*eta2)/169344.;
    coeffs->FfPN2Topf0_f01 = 15773941./254016. - (57635.*eta)/882. + 
        (23.*chi1*chi2*eta)/2. - (871.*eta2)/63. - 15.*chi12_2*eta2 - 
        (chi2_2*(23. - 23.*dm - 166.*eta + 120.*dm*eta))/8. + (chi1_2*(-23. + 
        166.*eta + dm*(-23. + 120.*eta)))/8.;
    return;
}

INT CalculateEccCorrectionToFluxV3X(REAL8 r, REAL8 vr, REAL8 prDot,
                                     REAL8 *cFr, REAL8 *cFf, REAL8 e0,
                                     EccCorrectionCoeffs *coeffs)
{
    if (fabs(e0) < 1e-3)
    {
        *cFr = 1.0;
        *cFf = 1.0;
        return CEV_SUCCESS;
    }
    REAL8 h0, g0, g0_2, g0_3, g0_4;
    REAL8 h1, f1, f1_2, f1_3;
    REAL8 h2, f2, f2_2, f2_3, f2_4;
    // REAL8 chi1_2, chi2_2, eta_2;
    h0 = 1./r;
    h1 = vr*vr;
    h2 = r*prDot;
    g0 = sqrt(h0);
    g0_2 = h0;
    g0_3 = g0_2*g0;
    g0_4 = g0_3*g0;
    f1 = h1 / h0;
    f1_2 = f1*f1;
    f1_3 = f1_2*f1;
    f2 = h2 / h0;
    f2_2 = f2*f2;
    f2_3 = f2_2*f2;
    // f2_4 = f2_3*f2;
    // chi1_2 = chi1*chi1;
    // chi2_2 = chi2*chi2;
    // eta_2 = eta*eta;
    REAL8 corrFrPN0, corrFfPN0;
    // REAL8 BaseFr = 1., BaseFf = 1.;
    // REAL8 dm;
    // dm = sqrt(1.-4.*eta);
    REAL8    sqf2, sqf2_2;
    sqf2 = sqrt(GET_MAX(1.+f2, 0.0001));
    sqf2_2 = sqf2*sqf2;
    // sqf2_3 = sqf2_2 * sqf2;
    // sqf2_4 = sqf2_2 * sqf2_2;
    // sqf2_5 = sqf2_4 * sqf2;
    // sqf2_6 = sqf2_5 * sqf2;
    REAL8 invsqf2_2 = 1./sqf2_2;
    corrFrPN0 = (6 + 7*f1 - 5*f2)/6.;
    corrFfPN0 = 1 + (29*f1)/12. - (5*f2)/6.;
    REAL8 corrFfPN1 = 0., corrFfPN32 = 0., corrFfPN2 = 0.;
    REAL8 corrFrPN1 = 0., corrFrPN32 = 0., corrFrPN2 = 0.;
    REAL8 FrPN2T0 = 0., FrPN2Topfm1 = 0., FrPN2Topfm2 = 0.;
    REAL8 FfPN2T0 = 0., FfPN2Topfm1 = 0., FfPN2Topfm2 = 0.;

    corrFrPN1 = coeffs->FrPN1_f01 * f2 + coeffs->FrPN1_f02 * f2_2 + coeffs->FrPN1_f03 * f2_3 +
        (coeffs->FrPN1_f10 + coeffs->FrPN1_f11 * f2 + coeffs->FrPN1_f12 * f2_2)*f1 +
        (coeffs->FrPN1_f20 + coeffs->FrPN1_f21 * f2) * f1_2;
    corrFrPN1 /=  (sqf2_2*(6.+7.*f1-5.*f2));
    corrFfPN1 = coeffs->FfPN1_f01 * f2 + coeffs->FfPN1_f02 * f2_2 + coeffs->FfPN1_f03 * f2_3 +
        (coeffs->FfPN1_f10 + coeffs->FfPN1_f11 * f2 + coeffs->FfPN1_f12 * f2_2)*f1 +
        (coeffs->FfPN1_f20 + coeffs->FfPN1_f21 * f2) * f1_2;
    corrFfPN1 /= (sqf2_2*(12.+29.*f1-10.*f2));
// REAL8 tmp = coeffs->FfPN1_f01 * f2 + coeffs->FfPN1_f02 * f2_2 + coeffs->FfPN1_f03 * f2_3 +
//         (coeffs->FfPN1_f10 + coeffs->FfPN1_f11 * f2 + coeffs->FfPN1_f12 * f2_2)*f1 +
//         (coeffs->FfPN1_f20 + coeffs->FfPN1_f21 * f2) * f1_2;
// print_debug("Num = %.16e\n", tmp);
// print_debug("FfPN1_f01 = %.16e\n", coeffs->FfPN1_f01);
// print_debug("FfPN1_f02 = %.16e\n", coeffs->FfPN1_f02);
// print_debug("FfPN1_f03 = %.16e\n", coeffs->FfPN1_f03);

// print_debug("FfPN1_f10 = %.16e\n", coeffs->FfPN1_f10);
// print_debug("FfPN1_f11 = %.16e\n", coeffs->FfPN1_f11);
// print_debug("FrPN1_f12 = %.16e\n", coeffs->FrPN1_f12);

// print_debug("FfPN1_f20 = %.16e\n", coeffs->FfPN1_f20);
// print_debug("FfPN1_f21 = %.16e\n", coeffs->FfPN1_f21);


    corrFrPN32 = coeffs->FrPN32_f01 * f2 + coeffs->FrPN32_f02 * f2_2 +
        (coeffs->FrPN32_f10 + coeffs->FrPN32_f11 * f2) * f1 +
        (coeffs->FrPN32_f20 + coeffs->FrPN32_f21 * f2) * f1_2;
    corrFrPN32 /= (sqf2*(-6.-7.*f1+5.*f2));
    corrFfPN32 = coeffs->FfPN32_f01 * f2 + coeffs->FfPN32_f02 * f2_2 +
        (coeffs->FfPN32_f10 + coeffs->FfPN32_f11 * f2) * f1 +
        (coeffs->FfPN32_f20 + coeffs->FfPN32_f30 * f1) * f1_2;
    corrFfPN32 /= (sqf2*(-12.-29.*f1+10.*f2));

    FrPN2T0 = coeffs->FrPN2Topf0_f00 +
        coeffs->FrPN2Topf0_f01*f2 + coeffs->FrPN2Topf0_f02*f2_2 + coeffs->FrPN2Topf0_f03*f2_3 +
        (coeffs->FrPN2Topf0_f10 + coeffs->FrPN2Topf0_f11*f2 + coeffs->FrPN2Topf0_f12*f2_2)*f1 + 
        (coeffs->FrPN2Topf0_f20 + coeffs->FrPN2Topf0_f21*f2) * f1_2 +
        coeffs->FrPN2Topf0_f30 * f1_3;
    FrPN2Topfm1 = coeffs->FrPN2Topfm1_f00 + coeffs->FrPN2Topfm1_f10 * f1 + coeffs->FrPN2Topfm1_f20 * f1_2;
    FrPN2Topfm2 = coeffs->FrPN2Topfm2_f00 + coeffs->FrPN2Topfm2_f10 * f1 + coeffs->FrPN2Topfm2_f20 * f1_2;
    corrFrPN2 = (FrPN2T0 + (FrPN2Topfm1 + FrPN2Topfm2*invsqf2_2)*invsqf2_2) / (6. + 7.*f1 - 5.*f2);

    FfPN2T0 = coeffs->FfPN2Topf0_f01*f2 + coeffs->FfPN2Topf0_f02*f2_2 + coeffs->FfPN2Topf0_f03*f2_3 +
        (coeffs->FfPN2Topf0_f10 + coeffs->FfPN2Topf0_f11*f2 + coeffs->FfPN2Topf0_f12*f2_2)*f1 + 
        (coeffs->FfPN2Topf0_f20 + coeffs->FfPN2Topf0_f21*f2) * f1_2 +
        coeffs->FfPN2Topf0_f30 * f1_3;
    FfPN2Topfm1 = (coeffs->FfPN2Topfm1_f04 * f2_2 + coeffs->FfPN2Topfm1_f13 * f2 * f1 + coeffs->FfPN2Topfm1_f22 * f1_2) * f2_2;
    FfPN2Topfm2 = (coeffs->FfPN2Topfm2_f05 * f2_2 + coeffs->FfPN2Topfm2_f14 * f2 * f1 + coeffs->FfPN2Topfm2_f23 * f1_2) * f2_3;
    corrFfPN2 = (FfPN2T0 + (FfPN2Topfm1 + FfPN2Topfm2*invsqf2_2)*invsqf2_2) / (12. + 29.*f1 - 10.*f2);

// print_err("\n");
// print_debug("FrPN2Coeffs:\n");
// print_err("T0:\n\tf01,f02,f03 = %.16e, %.16e, %.16e\n\tf10,f11,f12 = %.16e, %.16e, %.16e\n\tf20,f21 = %.16e, %.16e\n\tf30 = %.16e\n",
//     coeffs->FrPN2Topf0_f01, coeffs->FrPN2Topf0_f02, coeffs->FrPN2Topf0_f03,
//     coeffs->FrPN2Topf0_f10, coeffs->FrPN2Topf0_f11, coeffs->FrPN2Topf0_f12,
//     coeffs->FrPN2Topf0_f20, coeffs->FrPN2Topf0_f21,
//     coeffs->FfPN2Topf0_f30);
// print_err("fm1:\n\tf00,f10,f20 = %.16e, %.16e, %.16e\n",
//     coeffs->FrPN2Topfm1_f00, coeffs->FrPN2Topfm1_f10, coeffs->FrPN2Topfm1_f20);
// print_err("fm2:\n\tf00,f10,f20 = %.16e, %.16e, %.16e\n",
//     coeffs->FrPN2Topfm2_f00, coeffs->FrPN2Topfm2_f10, coeffs->FrPN2Topfm2_f20);

// print_debug("FfPN2Coeffs:\n");
// print_err("T0:\n\tf01,f02,f03 = %.16e, %.16e, %.16e\n\tf10,f11,f12 = %.16e, %.16e, %.16e\n\tf20,f21 = %.16e, %.16e\b\tf30 = %.16e\n",
//     coeffs->FfPN2Topf0_f01, coeffs->FfPN2Topf0_f02, coeffs->FfPN2Topf0_f03,
//     coeffs->FfPN2Topf0_f10, coeffs->FfPN2Topf0_f11, coeffs->FfPN2Topf0_f12,
//     coeffs->FfPN2Topf0_f20, coeffs->FfPN2Topf0_f21,
//     coeffs->FfPN2Topf0_f30);
// print_err("fm1:\n\tf04,f13,f22 = %.16e, %.16e, %.16e\n",
//     coeffs->FfPN2Topfm1_f04, coeffs->FfPN2Topfm1_f13, coeffs->FfPN2Topfm1_f22);
// print_err("fm2:\n\tf05,f14,f23 = %.16e, %.16e, %.16e\n",
//     coeffs->FfPN2Topfm2_f05, coeffs->FfPN2Topfm2_f14, coeffs->FfPN2Topfm2_f23);

// print_debug("r, vr, prDot = %.16e, %.16e, %.16e\n", r, vr, prDot);
// print_debug("FrPN0 = %.16e, FrPN1 = %.16e, FrPN32 = %.16e, FrPN2 = %.16e\n",
//         corrFrPN0, corrFrPN1, corrFrPN32, corrFrPN2);
// print_debug("FfPN0 = %.16e, FfPN1 = %.16e, FfPN32 = %.16e, FfPN2 = %.16e\n",
//         corrFfPN0, corrFfPN1, corrFfPN32, corrFfPN2);
// print_err("\n");

    *cFf = corrFfPN0 / (1 - corrFfPN1 * g0_2 - corrFfPN32 * g0_3 + (corrFfPN1*corrFfPN1 - corrFfPN2)*g0_4);
    *cFr = corrFrPN0 / (1 - corrFrPN1 * g0_2 - corrFrPN32 * g0_3 + (corrFrPN1*corrFrPN1 - corrFrPN2)*g0_4);
    return CEV_SUCCESS;
}


INT CalculateEccCorrectionToFluxV3(REAL8 eta, REAL8 chi1, REAL8 chi2,
                                     REAL8 r, REAL8 vr, REAL8 prDot,
                                     REAL8 *cFr, REAL8 *cFf, REAL8 e0)
{
    if (fabs(e0) < 1e-3)
    {
        *cFr = 1.0;
        *cFf = 1.0;
        return CEV_SUCCESS;
    }
// print_debug("r, vr, prDot = %.16e, %.16e, %.16e\n", r, vr, prDot);
    REAL8 h0, g0, g0_2, g0_3, g0_4;
    REAL8 h1, f1, f1_2, f1_3;
    REAL8 h2, f2, f2_2, f2_3, f2_4;
    REAL8 chi1_2, chi2_2, eta_2;
    h0 = 1./r;
    h1 = vr*vr;
    h2 = r*prDot;
    g0 = sqrt(h0);
    g0_2 = h0;
    g0_3 = g0_2*g0;
    g0_4 = g0_3*g0;
    f1 = h1 / h0;
    f1_2 = f1*f1;
    f1_3 = f1_2*f1;
    f2 = h2 / h0;
    f2_2 = f2*f2;
    f2_3 = f2_2*f2;
    f2_4 = f2_3*f2;
    chi1_2 = chi1*chi1;
    chi2_2 = chi2*chi2;
    eta_2 = eta*eta;
    REAL8 corrFrPN0, corrFfPN0;
    REAL8 BaseFr = 1., BaseFf = 1.;
    REAL8 dm;
    dm = sqrt(1.-4.*eta);
    REAL8    sqf2, sqf2_2, sqf2_3, sqf2_4, sqf2_5, sqf2_6;
    sqf2 = sqrt(GET_MAX(1.+f2, 0.0001));
    sqf2_2 = sqf2*sqf2;
    sqf2_3 = sqf2_2 * sqf2;
    sqf2_4 = sqf2_2 * sqf2_2;
    sqf2_5 = sqf2_4 * sqf2;
    sqf2_6 = sqf2_5 * sqf2;
    corrFrPN0 = (6 + 7*f1 - 5*f2)/6.;
    corrFfPN0 = 1 + (29*f1)/12. - (5*f2)/6.;
    REAL8 corrFfPN1 = 0., corrFfPN32 = 0., corrFfPN2 = 0.;
    REAL8 corrFrPN1 = 0., corrFrPN32 = 0., corrFrPN2 = 0.;
    REAL8 chi12_2;
    chi12_2 = (chi1+chi2) * (chi1+chi2);

        corrFrPN1 = ((12*(-181 + 86*eta)*f1_2*(1 + f2) + f2*(-1714 + f2*(2029 + 3919*f2) - 4*eta*(-804 + f2*(67 + 1047*f2))) + 
            f1*(-16018 - f2*(22247 + 6341*f2) + 8*eta*(-297 + f2*(1003 + 1356*f2)))))/(336.*(6 + 7*f1 - 5*f2)*(1 + f2));
        corrFfPN1 = ((-24*(-149 + 124*eta)*f1_2*(1 + f2) + 2*f2*(-1714 + f2*(2029 + 3919*f2) - 4*eta*(-804 + f2*(67 + 1047*f2))) + 
            f1*(14713 + (3680 - 11497*f2)*f2 + 4*eta*(-12657 + f2*(-5762 + 7359*f2)))))/(336.*(12 + 29*f1 - 10*f2)*(1 + f2));
// REAL8 tmpNum, tmpDen;
// tmpNum = ((-24*(-149 + 124*eta)*f1_2*(1 + f2) + 2*f2*(-1714 + f2*(2029 + 3919*f2) - 4*eta*(-804 + f2*(67 + 1047*f2))) + 
//             f1*(14713 + (3680 - 11497*f2)*f2 + 4*eta*(-12657 + f2*(-5762 + 7359*f2)))));
// tmpDen = (336.*(12 + 29*f1 - 10*f2)*(1 + f2));
// print_debug("num = %.16e, den = %.16e\n", tmpNum, tmpDen);

        corrFrPN32 = -((8*chi2*f2*(eta*(35 + 57*f2) + (-1 + dm)*(13 + 79*f2)) - 4*chi2*f1*(eta*(50 + 78*f2) + (-1 + dm)*(247 + 331*f2)) + 
            4*chi1*(-2*eta*f1*(25 + 39*f2) + 2*eta*f2*(35 + 57*f2) - 2*(1 + dm)*f2*(13 + 79*f2) + (1 + dm)*f1*(247 + 331*f2)) + 
            (12*(185 - 167*f2)*f2 + 5*f1_2*(1 + f2) + 32*f1*(293 + 377*f2))*CST_PI))/(96.*(-6 - 7*f1 + 5*f2)*sqf2);

        corrFfPN32 = ((-5*chi2*(180*(-2 + 2*dm + eta)*f1_2 + eta*f1*(80 - 516*f2) + 4*eta*f2*(35 + 57*f2) + 4*(-1 + dm)*f2*(13 + 79*f2) - 
            (-1 + dm)*f1*(-427 + 857*f2)) + 5*chi1*(180*(2 + 2*dm - eta)*f1_2 - 4*eta*f2*(35 + 57*f2) + 4*(1 + dm)*f2*(13 + 79*f2) + 
            eta*f1*(-80 + 516*f2) - (1 + dm)*f1*(-427 + 857*f2)) + 
            (-10020*f1_2 + 49*f1_3 + 30*f2*(-185 + 167*f2) - 30*f1*(216 + 167*f2))*CST_PI))/(120.*(-12 - 29*f1 + 10*f2)*sqf2);

        corrFrPN2 = ((8*(2927618 - 24830856*eta + 5069352*eta_2 - 73984547*f1 + 41891769*eta*f1 + 5859372*eta_2*f1 - 50108193*f1_2 + 
            185971113*eta*f1_2 + 7780590*eta_2*f1_2 - 7345170*f1_3 + 40452426*eta*f1_3 + 4240404*eta_2*f1_3 + 
            15876*chi1*chi2*eta*(-398 + 187*f1 + 8*eta*(59 + 3*f1)) + 
            3969*chi2_2*(46 + 125*f1 + 48*eta_2*(5 + f1) - 2*eta*(166 + 149*f1) + dm*(-46 - 125*f1 + 48*eta*(5 + f1))) - 
            3969*chi1_2*(-46 - 125*f1 - 48*eta_2*(5 + f1) + eta*(332 + 298*f1) + dm*(-46 - 125*f1 + 48*eta*(5 + f1)))) + 
            (1152*(-1 + 4*eta)*(11 + 7*f1)*(-44 - 34*eta + 21*(-3 + eta)*f1))/sqf2_4 - 
            (48*(-1 + 4*eta)*(-6803 - 20756*eta + f1*(36983 + 22256*eta + 108*(93 + 4*eta)*f1)))/(1 + f2) + 
            (-71208481 - 1460592*chi1_2*(1 + dm) + 24*(60858*chi2_2*(-1 + dm) + 
            9*(1445957 + 588*(46*chi1*chi2 + chi2_2*(83 - 60*dm) + chi1_2*(83 + 60*dm)))*eta - 
            70*(68249 + 4536*chi12_2)*eta_2) - 63*(2875941 + 4*eta*(-7232999 + 294428*eta))*f1 + 
            108*(-719923 + 4400322*eta + 555796*eta_2)*f1_2)*(1 + f2) + 
            9*(4347176 + 16*eta*(85017 + 745502*eta) + (1403989 - 4613140*eta + 7636304*eta_2)*f1)*sqf2_4 - 
            9*(-936881 + 8*eta*(1728233 + 493202*eta))*sqf2_6))/(1.016064e6*(6 + 7*f1 - 5*f2));

        corrFfPN2 = ((7461720*chi1_2*dm*f1 - 7461720*chi2_2*dm*f1 + 4572288*chi1_2*dm*eta*f1 - 4572288*chi2_2*dm*eta*f1 + 
            (-676564289 + 7461720*chi1_2 + 72*(103635*chi2_2 + (7460071 - 882*(163*chi1_2 - 1018*chi1*chi2 + 163*chi2_2))*eta - 
            294*(-13249 + 216*chi12_2)*eta_2))*f1 + 63095764*f2 - 2921184*chi1_2*dm*f2 + 2921184*chi2_2*dm*f2 + 
            15240960*chi1_2*dm*eta*f2 - 15240960*chi2_2*dm*eta*f2 + 
            288*(-10143*(chi1_2 + chi2_2) + (-230540 + 882*(83*chi1_2 + 46*chi1*chi2 + 83*chi2_2))*eta - 
            56*(871 + 945*chi12_2)*eta_2)*f2 + 
            (1152*(-1 + 4*eta)*f2_3*(-29*f1 + 22*f2)*(21*(-3 + eta)*f1 + 2*(22 + 17*eta)*f2))/sqf2_4 - 
            (96*(-1 + 4*eta)*f2_2*(-108*(-627 + 94*eta)*f1_2 - (107093 + 100532*eta)*f1*f2 + (51277 + 24124*eta)*f2_2))/(1 + f2) + 
            6*(12*(-9410885 + 6*eta*(5667677 + 876428*eta))*f1_2 - 4*(8917985 + 2*eta*(-50196183 + 876830*eta))*f1*f2 + 
            (21922177 - 8*eta*(15256011 + 15182*eta))*f2_2) + 
            3*(672*(-55655 + 324681*eta + 38214*eta_2)*f1_3 + 72*(-788870 + eta*(4644559 + 432720*eta))*f1_2*f2 + 
            (16230547 + 88*eta*(-580353 + 408734*eta))*f1*f2_2 - 2*(-2176067 + 8*eta*(4897707 + 1358422*eta))*f2_3)))/
            (1.016064e6*(12 + 29*f1 - 10*f2));
// print_err("\n");
// print_debug("eta = %.16e, chi1 = %.16e, chi2 = %.16e\n", eta, chi1, chi2);
// print_debug("r, vr, prDot = %.16e, %.16e, %.16e\n", r, vr, prDot);
// print_debug("FrPN0 = %.16e, FrPN1 = %.16e, FrPN32 = %.16e, FrPN2 = %.16e\n",
//         corrFrPN0, corrFrPN1, corrFrPN32, corrFrPN2);
// print_debug("FfPN0 = %.16e, FfPN1 = %.16e, FfPN32 = %.16e, FfPN2 = %.16e\n",
//         corrFfPN0, corrFfPN1, corrFfPN32, corrFfPN2);
// print_err("\n");
    *cFf = corrFfPN0 / (1 - corrFfPN1 * g0_2 - corrFfPN32 * g0_3 + (corrFfPN1*corrFfPN1 - corrFfPN2)*g0_4);
    *cFr = corrFrPN0 / (1 - corrFrPN1 * g0_2 - corrFrPN32 * g0_3 + (corrFrPN1*corrFrPN1 - corrFrPN2)*g0_4);
    return CEV_SUCCESS;
}

INT CalculateEccCorrectionToFluxV4(REAL8 eta, REAL8 chi1, REAL8 chi2,
                                     REAL8 r, REAL8 prt, REAL8 prDot,
                                     REAL8 *cFr, REAL8 *cFf)
{
    REAL8 h0, g0, g0_2, g0_3, g0_4;
    REAL8 h1, f1, f1_2, f1_3;
    REAL8 h2, f2, f2_2, f2_3, f2_4;
    REAL8 chi1_2, chi2_2, eta_2;
    h0 = 1./r;
    h1 = prt*prt;
    h2 = r*prDot;
    g0 = sqrt(h0);
    g0_2 = h0;
    g0_3 = g0_2*g0;
    g0_4 = g0_3*g0;
    f1 = h1 / h0;
    f1_2 = f1*f1;
    f1_3 = f1_2*f1;
    f2 = h2 / h0;
    f2_2 = f2*f2;
    f2_3 = f2_2*f2;
    f2_4 = f2_3*f2;
    chi1_2 = chi1*chi1;
    chi2_2 = chi2*chi2;
    eta_2 = eta*eta;
    REAL8 corrFrPN0, corrFfPN0;
    REAL8 BaseFr = 1., BaseFf = 1.;
    REAL8 dm;
    dm = sqrt(1.-4.*eta);
    REAL8    sqf2, sqf2_2, sqf2_3, sqf2_4, sqf2_5, sqf2_6;
    sqf2 = sqrt(GET_MAX(1.+f2, 0.0001));
    sqf2_2 = sqf2*sqf2;
    sqf2_3 = sqf2_2 * sqf2;
    sqf2_4 = sqf2_2 * sqf2_2;
    sqf2_5 = sqf2_4 * sqf2;
    sqf2_6 = sqf2_5 * sqf2;
    REAL8 corrFfPN1 = 0., corrFfPN32 = 0., corrFfPN2 = 0.;
    REAL8 corrFrPN1 = 0., corrFrPN32 = 0., corrFrPN2 = 0.;
    REAL8 chi12_2;
    corrFrPN0 = (6 + 7*f1 - 5*f2)/6.;
    corrFfPN0 = 1 + (29*f1)/12. - (5*f2)/6.;
    chi12_2 = (chi1+chi2) * (chi1+chi2);

        corrFrPN1 = -((12*(377 + 110*eta)*f1_2*(1 + f2) + f1*(23074 + f2*(31655 + 8693*f2) - 8*eta*(-3 + 59*f2*(17 + 18*f2))) + 
            f2*(1714 - f2*(2029 + 3919*f2) + 4*eta*(-804 + f2*(67 + 1047*f2)))))/(336.*(6 + 7*f1 - 5*f2)*sqf2_2);
        corrFfPN1 = -((24*(257 + 530*eta)*f1_2*(1 + f2) + f1*(14519 + f2*(35296 + 21241*f2) + 4*eta*(10221 + (5762 - 4923*f2)*f2)) + 
            2*f2*(1714 - f2*(2029 + 3919*f2) + 4*eta*(-804 + f2*(67 + 1047*f2)))))/(336.*(12 + 29*f1 - 10*f2)*sqf2_2);

        corrFrPN32 = -((8*chi2*f2*(eta*(35 + 57*f2) + (-1 + dm)*(13 + 79*f2)) - 4*chi2*f1*(eta*(50 + 78*f2) + (-1 + dm)*(247 + 331*f2)) + 
            4*chi1*(-2*eta*f1*(25 + 39*f2) + 2*eta*f2*(35 + 57*f2) - 2*(1 + dm)*f2*(13 + 79*f2) + (1 + dm)*f1*(247 + 331*f2)) + 
            (12*(185 - 167*f2)*f2 + 5*f1_2*(1 + f2) + 32*f1*(293 + 377*f2))*CST_PI))/(96.*(-6 - 7*f1 + 5*f2)*sqf2);

        corrFfPN32 = ((-5*chi2*(180*(-2 + 2*dm + eta)*f1_2 + eta*f1*(80 - 516*f2) + 4*eta*f2*(35 + 57*f2) + 4*(-1 + dm)*f2*(13 + 79*f2) - 
            (-1 + dm)*f1*(-427 + 857*f2)) + 5*chi1*(180*(2 + 2*dm - eta)*f1_2 - 4*eta*f2*(35 + 57*f2) + 4*(1 + dm)*f2*(13 + 79*f2) + 
            eta*f1*(-80 + 516*f2) - (1 + dm)*f1*(-427 + 857*f2)) + 
            (-10020*f1_2 + 49*f1_3 + 30*f2*(-185 + 167*f2) - 30*f1*(216 + 167*f2))*CST_PI))/(120.*(-12 - 29*f1 + 10*f2)*sqf2);

        corrFrPN2 = ((44706816*chi1*chi2*(-1 + eta)*eta + 11081448*chi1_2*dm*f1 - 11081448*chi2_2*dm*f1 - 15748992*chi1_2*dm*eta*f1 + 
            15748992*chi2_2*dm*eta*f1 - 1460592*chi1_2*dm*f2 + 1460592*chi2_2*dm*f2 + 7620480*chi1_2*dm*eta*f2 - 
            7620480*chi2_2*dm*eta*f2 + (1152*(-1 + 4*eta)*f2_3*(-7*f1 + 11*f2)*(21*(-3 + eta)*f1 + 2*(22 + 17*eta)*f2))/sqf2_4 + 
            2*((-306090563 + 5540724*chi1_2 + 36*(153909*chi2_2 + (29094643 - 882*(597*chi1_2 - 374*chi1*chi2 + 597*chi2_2))*eta + 
            28*(18367 + 7812*chi12_2)*eta_2))*f1 + 15773941*f2 + 
            72*(-10143*(chi1_2 + chi2_2) + (-230540 + 882*(83*chi1_2 + 46*chi1*chi2 + 83*chi2_2))*eta - 
            56*(871 + 945*chi12_2)*eta_2)*f2) - 
            (48*(-1 + 4*eta)*f2_2*(-36*(-1357 + 86*eta)*f1_2 - (147191 + 8816*eta)*f1*f2 + (51277 + 24124*eta)*f2_2))/sqf2_2 + 
            3*(36*(-3609073 + 8*eta*(2290626 + 141281*eta))*f1_2 + (-14792459 + 148*eta*(3846909 + 215228*eta))*f1*f2 + 
            (21922177 - 8*eta*(15256011 + 15182*eta))*f2_2) + 
            3*(336*(-38207 + 336711*eta + 34518*eta_2)*f1_3 + 36*(-338967 + 4504590*eta + 293812*eta_2)*f1_2*f2 + 
            (13759775 + 4*eta*(-6383523 + 3151516*eta))*f1*f2_2 + (2176067 - 8*eta*(4897707 + 1358422*eta))*f2_3)))/
            (1.016064e6*(6 + 7*f1 - 5*f2));

        corrFfPN2 = ((36927576*chi1_2*dm*f1 - 36927576*chi2_2*dm*f1 - 54359424*chi1_2*dm*eta*f1 + 54359424*chi2_2*dm*eta*f1 + 
            (-810040625 + 36927576*chi1_2 + 72*(512883*chi2_2 + (15377953 - 882*(2019*chi1_2 - 1018*chi1*chi2 + 2019*chi2_2))*eta + 
            42*(51859 + 17976*chi12_2)*eta_2))*f1 + 63095764*f2 - 2921184*chi1_2*dm*f2 + 2921184*chi2_2*dm*f2 + 
            15240960*chi1_2*dm*eta*f2 - 15240960*chi2_2*dm*eta*f2 + 
            288*(-10143*(chi1_2 + chi2_2) + (-230540 + 882*(83*chi1_2 + 46*chi1*chi2 + 83*chi2_2))*eta - 
            56*(871 + 945*chi12_2)*eta_2)*f2 + 
            (1152*(-1 + 4*eta)*f2_3*(-29*f1 + 22*f2)*(21*(-3 + eta)*f1 + 2*(22 + 17*eta)*f2))/sqf2_4 - 
            (96*(-1 + 4*eta)*f2_2*(36*(2287 + 124*eta)*f1_2 - 25*(5453 + 2852*eta)*f1*f2 + (51277 + 24124*eta)*f2_2))/sqf2_2 + 
            6*(12*(-10929983 + 12*eta*(3030017 + 526372*eta))*f1_2 + 4*(-4759733 + 94943874*eta + 6258932*eta_2)*f1*f2 + 
            (21922177 - 8*eta*(15256011 + 15182*eta))*f2_2) + 
            3*(672*(-51767 + 341151*eta + 61758*eta_2)*f1_3 + 72*(-529912 + eta*(4600725 + 334552*eta))*f1_2*f2 + 
            (33665923 + 8*eta*(-8921901 + 1868722*eta))*f1*f2_2 - 2*(-2176067 + 8*eta*(4897707 + 1358422*eta))*f2_3)))/
            (1.016064e6*(12 + 29*f1 - 10*f2));

    *cFf = corrFfPN0 / (1 - corrFfPN1 * g0_2 - corrFfPN32 * g0_3 + (corrFfPN1*corrFfPN1 - corrFfPN2)*g0_4);
    *cFr = corrFrPN0 / (1 - corrFrPN1 * g0_2 - corrFrPN32 * g0_3 + (corrFrPN1*corrFrPN1 - corrFrPN2)*g0_4);
    return CEV_SUCCESS;
}