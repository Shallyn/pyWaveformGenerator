/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "newFactorizedWaveformPrec.h"
#include "myLog.h"
#include "pHamiltonian.h"
#include "pFactorizedWaveform.h"
#include <gsl/gsl_sf_gamma.h>

/**
 * This function calculates coefficients for hlm mode factorized-resummed waveform.
 * The coefficients are pre-computed and stored in the SpinEOBParams structure.
 * Eq. 17 and the entire Appendix of PRD 86, 024011 (2012) + changes
 * described in ￼the section "Factorized waveforms" of https://dcc.ligo.org/T1400476
 */

int
EccPrec_CalcSpinPrecFacWaveformCoefficients(
                        FacWaveformCoeffs * const coeffs,	/**< OUTPUT, pre-computed waveform coefficients */
                        const REAL8 m1,	/**< mass 1 */
                        const REAL8 m2,	/**< mass 2 */
                        const REAL8 eta,	/**< symmetric mass ratio */
                        const REAL8 tmpa,	/**< Kerr spin parameter for test-particle terms */
                        const REAL8 chiS,
                        const REAL8 chiA,
                        const REAL8 csx,	/**< (chi1+chi2)/2 */
                        const REAL8 csy,
                        const REAL8 csz,
                        const REAL8 cax,	/**< (chi1-chi2)/2 */
                        const REAL8 cay,
                        const REAL8 caz,
                        UINT SpinAlignedEOBversion	/**< 1 for SEOBNRv1; 2 for SEOBNRv2; 4 for the coefficients
                        in the flux of v4P and v4Pwave for the coefficients in the waveform of v4P */
                        )
{
	// int		debugPK = 0;
	// if (debugPK) {
	// 	XLAL_PRINT_INFO("In XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients: Renewing hlm coefficients.\n");
	// 		XLAL_PRINT_INFO("PK:: chiS = %.12e, chiA = %.12e, a = %.12e (UNUSED), EOBVERSION = %d\n",
	// 			   chiS, chiA, tmpa, SpinAlignedEOBversion);
	// }
	REAL8	a = tmpa;
	INT waveform = 0;
	REAL8	eta2 = eta * eta;
	REAL8	eta3 = eta2 * eta;

	if (SpinAlignedEOBversion == 451) {
		SpinAlignedEOBversion = 4;
		waveform = 1;
	}
	 //RC: Here I am taking care of the fact that the flux used to calibrate SEOBNRv4 is missing some terms
	 //in the 21 mode that I have added later in SEOBNRv4HM. For this reason when computing the Coefficients
	 //for the 21 mode one has to specify if they are needed for the flux or for the waveform.

	REAL8		dM	  , dM2, chiA2, chiS2;
	//dM3;
	REAL8		aDelta  , a2, a3;

	/* Combination which appears a lot */
	REAL8		m1Plus3eta, m1Plus3eta2;
	REAL8 m1Plus3eta3;

	dM2 = 1. - 4. * eta;
    // REAL8 chiA = caz;
    // REAL8 chiS = csz;
	chiA2 = chiA*chiA;
	chiS2 = chiS*chiS;
	REAL8 chiS3 = chiS2*chiS;
	REAL8 chiA3 = chiA2*chiA;


	//XLAL_PRINT_INFO("****************************** a = %e *********************************\n", a);

	/* Check that deltaM has a reasonable value */
	if (dM2 < 0) 
    {
		PRINT_LOG_INFO(LOG_CRITICAL, "eta seems to be > 0.25 - this isn't allowed!");
		return CEV_FAILURE;
	}
	dM = sqrt(dM2);
	if (m1 < m2) 
    {
		dM = -dM;
	}
	//dM3 = dM2 * dM;

	aDelta = 0.;
	//a value in delta_lm is 0 in both SEOBNRv1 and SEOBNRv2
	a2 = a * a;
	a3 = a2 * a;

	m1Plus3eta = -1. + 3. * eta;
	m1Plus3eta2 = m1Plus3eta * m1Plus3eta;
	m1Plus3eta3 = m1Plus3eta * m1Plus3eta2;

	/* Initialize all coefficients to zero */
	/* This is important, as we will not set some if dM is zero */
	memset(coeffs, 0, sizeof(FacWaveformCoeffs));


	/*
	 * l = 2, Eqs. A8a and A8b for rho, Eq. A15a for f, Eqs. 20 and 21 of
	 * DIN and Eqs. 27a and 27b of PBFRT for delta as well as eqns 28-29 of PBFRT
	 */

	coeffs->delta22vh3 = 7. / 3.;
	coeffs->delta22vh6 = (-4. * aDelta) / 3. + (428. * CST_PI) / 105.;
	if (SpinAlignedEOBversion == 4)
    {
			coeffs->delta22vh6 =
	-4. / 3. * (dM * chiA + chiS * (1 - 2 * eta)) +
	(428. * CST_PI) / 105.;
		}
	coeffs->delta22v8 = (20. * aDelta) / 63.;
	coeffs->delta22vh9 = -2203. / 81. + (1712. * CST_PI * CST_PI) / 315.;
	coeffs->delta22v5 = -24. * eta;
	coeffs->delta22v6 = 0.0;
	if (SpinAlignedEOBversion == 2 && chiS + chiA * dM / (1. - 2. * eta) > 0.) {
		coeffs->delta22v6 = -540. * eta * (chiS + chiA * dM / (1. - 2. * eta));
//		double chi = (chiS + chiA * dM / (1. - 2. * eta));
//		coeffs->delta22v6 = eta*(1./4.*(1. - 1080.*chi - sqrt((1. - 1080.*chi)*(1. - 1080*chi) + 8.*(13.5 +270.*chi +13.5*chi*chi))));
		}
   if (SpinAlignedEOBversion == 3 /*&& chiS + chiA * dM / (1. - 2. * eta) > 0.*/) {
//				coeffs->delta22v6 = -540. * eta * (chiS + chiA * dM / (1. - 2. * eta));
	   double chi = (chiS + chiA * dM / (1. - 2. * eta));
	   coeffs->delta22v6 = eta*(1./4.*(1. - 1080.*chi - sqrt((1. - 1080.*chi)*(1. - 1080*chi) + 8.*(13.5 +270.*chi +13.5*chi*chi))));
   }
	coeffs->rho22v2 = -43. / 42. + (55. * eta) / 84.;
	coeffs->rho22v3 = (-2. * (chiS + chiA * dM - chiS * eta)) / 3.;
	switch (SpinAlignedEOBversion) {
	case 1:
		coeffs->rho22v4 = -20555. / 10584. + 0.5 * a2
			- (33025. * eta) / 21168. + (19583. * eta2) / 42336.;
		break;
	case 2:
		coeffs->rho22v4 = -20555. / 10584. + 0.5 * (chiS + chiA * dM) * (chiS + chiA * dM)
			- (33025. * eta) / 21168. + (19583. * eta2) / 42336.;
		break;
	case 3:
		coeffs->rho22v4 = -20555. / 10584. + 0.5 * (chiS + chiA * dM) * (chiS + chiA * dM)
		- (33025. * eta) / 21168. + (19583. * eta2) / 42336.;
		break;
	case 4:
		coeffs->rho22v4 =
			-20555. / 10584. + 0.5 * (chiS + chiA * dM) * (chiS + chiA * dM) -
			(33025. * eta) / 21168. + (19583. * eta2) / 42336.;
		break;
	default:
		PRINT_LOG_INFO(LOG_CRITICAL, "wrong SpinAlignedEOBversion value, must be 1, 2 or 3!");
		return CEV_FAILURE;
		break;
	}
	coeffs->rho22v5 = (-34. * a) / 21.;
if (SpinAlignedEOBversion == 4)
			{
				coeffs->rho22v5 =
		(-34. / 21. + 49. * eta / 18. + 209. * eta2 / 126.) * chiS +
		(-34. / 21. - 19. * eta / 42.) * dM * chiA;
			}
	coeffs->rho22v6 = 1556919113. / 122245200. + (89. * a2) / 252. - (48993925. * eta) / 9779616.
		- (6292061. * eta2) / 3259872. + (10620745. * eta3) / 39118464.
		+ (41. * eta * CST_PI * CST_PI) / 192.;
	coeffs->rho22v6l = -428. / 105.;
	coeffs->rho22v7 = (18733. * a) / 15876. + a * a2 / 3.;
	if (SpinAlignedEOBversion == 4)
			{
				coeffs->rho22v7 =
		a3 / 3. + chiA * dM * (18733. / 15876. + (50140. * eta) / 3969. +
							 (97865. * eta2) / 63504.) +
		chiS * (18733. / 15876. + (74749. * eta) / 5292. -
			(245717. * eta2) / 63504. + (50803. * eta3) / 63504.);
			}
	switch (SpinAlignedEOBversion) {
	case 1:
		coeffs->rho22v8 = eta * (-5.6 - 117.6 * eta) - 387216563023. / 160190110080. + (18353. * a2) / 21168. - a2 * a2 / 8.;
		break;
	case 2:
		coeffs->rho22v8 = -387216563023. / 160190110080. + (18353. * a2) / 21168. - a2 * a2 / 8.;
		break;
	case 3:
		coeffs->rho22v8 = -387216563023. / 160190110080. + (18353. * a2) / 21168. - a2 * a2 / 8.;
		break;
	case 4:
		coeffs->rho22v8 =
			-387216563023. / 160190110080. + (18353. * a2) / 21168. -
			a2 * a2 / 8.;
		break;
	default:
		PRINT_LOG_INFO(LOG_CRITICAL, "wrong SpinAlignedEOBversion value, must be 1, 2 or 3!");
		return CEV_FAILURE;
		break;
	}
	coeffs->rho22v8l = 9202. / 2205.;
	coeffs->rho22v10 = -16094530514677. / 533967033600.;
	coeffs->rho22v10l = 439877. / 55566.;

	// if (debugPK) {
	// 	XLAL_PRINT_INFO("\nPK:: dM, eta, chiS, chiA while renewing hlm coeffs: %e, %e, %e, %e\n",
	// 		   dM, eta, chiS, chiA);
	// 	XLAL_PRINT_INFO("PK:: Renewed rho-lm coeffs: v2 = %.16e, v3 = %.16e, v4 = %.16e, v5 = %.16e\nv6 = %.16e, v6l = %.16e v7 = %.16e v8 = %.16e, v8l = %.16e v10 = %.16e v10l = %.16e\n",
	// 		   coeffs->rho22v2, coeffs->rho22v3, coeffs->rho22v4, coeffs->rho22v5,
	// 		   coeffs->rho22v6, coeffs->rho22v6l, coeffs->rho22v7, coeffs->rho22v8,
	// 		 coeffs->rho22v8l, coeffs->rho22v10, coeffs->rho22v10l);
	// }
	if(waveform == 1){
		a = 0;
		a2 = 0;
		a3 = 0;
	}
	coeffs->delta21vh3 = 2. / 3.;
	coeffs->delta21vh6 = (-17. * aDelta) / 35. + (107. * CST_PI) / 105.;
	coeffs->delta21vh7 = (3. * aDelta * aDelta) / 140.;
	coeffs->delta21vh9 = -272. / 81. + (214. * CST_PI * CST_PI) / 315.;
	coeffs->delta21v5 = -493. * eta / 42.;
	coeffs->delta21v7 = 0.0;
	if (dM2) {

		//coeffs->rho21v1 = (-3. * (chiS + chiA / dM)) / (4.);
		coeffs->rho21v1 = 0.0;
		//coeffs->rho21v2 = -59. / 56 - (9. * chiAPlusChiSdM * chiAPlusChiSdM) / (32. * dM2) + (23. * eta) / 84.;
		switch (SpinAlignedEOBversion) {
			case 1:
				coeffs->rho21v2 = -59. / 56 + (23. * eta) / 84. - 9. / 32. * a2;
				coeffs->rho21v3 = 1177. / 672. * a - 27. / 128. * a3;
				break;
			case 2:
				coeffs->rho21v2 = -59. / 56 + (23. * eta) / 84.;
				coeffs->rho21v3 = 0.0;
				break;
			case 3:
				coeffs->rho21v2 = -59. / 56 + (23. * eta) / 84.;
				coeffs->rho21v3 = 0.0;
				break;
			case 4:
				coeffs->rho21v2 = -59. / 56 + (23. * eta) / 84.;
				coeffs->rho21v3 = 0.0;
				break;
			default:
                PRINT_LOG_INFO(LOG_CRITICAL, "wrong SpinAlignedEOBversion value, must be 1, 2 or 3!");
                return CEV_FAILURE;
				break;
		}
		/*
		 * coeffs->rho21v3   = (-567.*chiA*chiA*chiA -
		 * 1701.*chiA*chiA*chiS*dM + chiA*(-4708. + 1701.*chiS*chiS -
		 * 2648.*eta)*(-1. + 4.*eta) + chiS* dM3 *(4708. -
		 * 567.*chiS*chiS + 1816.*eta))/(2688.*dM3);
		 */
		coeffs->rho21v4 = -47009. / 56448. - (865. * a2) / 1792. - (405. * a2 * a2) / 2048. - (10993. * eta) / 14112.
			+ (617. * eta2) / 4704.;
		coeffs->rho21v5 = (-98635. * a) / 75264. + (2031. * a * a2) / 7168. - (1701. * a2 * a3) / 8192.;
		coeffs->rho21v6 = 7613184941. / 2607897600. + (9032393. * a2) / 1806336. + (3897. * a2 * a2) / 16384.
			- (15309. * a3 * a3) / 65536.;
		coeffs->rho21v6l = -107. / 105.;
		coeffs->rho21v7 = (-3859374457. * a) / 1159065600. - (55169. * a3) / 16384.
			+ (18603. * a2 * a3) / 65536. - (72171. * a2 * a2 * a3) / 262144.;
		coeffs->rho21v7l = 107. * a / 140.;
		coeffs->rho21v8 = -1168617463883. / 911303737344.;
		coeffs->rho21v8l = 6313. / 5880.;
		coeffs->rho21v10 = -63735873771463. / 16569158860800.;
		coeffs->rho21v10l = 5029963. / 5927040.;

		coeffs->f21v1 = (-3. * (chiS + chiA / dM)) / (2.);
		switch (SpinAlignedEOBversion) {
		case 1:
			coeffs->f21v3 = 0.0;
			break;
		case 2:
			coeffs->f21v3 = (chiS * dM * (427. + 79. * eta) + chiA * (147. + 280. * dM * dM + 1251. * eta)) / 84. / dM;
			break;
		case 3:
			coeffs->f21v3 = (chiS * dM * (427. + 79. * eta) + chiA * (147. + 280. * dM * dM + 1251. * eta)) / 84. / dM;
			break;
		case 4:
			coeffs->f21v3 =
			(chiS * dM * (427. + 79. * eta) +
				chiA * (147. + 280. * dM * dM + 1251. * eta)) / 84. / dM;
			/*RC: New terms for SEOBNRv4HM, they are put to zero if use_hm == 0 */
			if (waveform == 0)
			{
				coeffs->f21v4 = 0.0;
				coeffs->f21v5 = 0.0;
				coeffs->f21v6 = 0.0;
				coeffs->f21v7c = 0;
			}
			else{
				coeffs->f21v4 = (-3.-2.*eta)*chiA2 + (-3.+ eta/2.)*chiS2 + (-6.+21.*eta/2.)*chiS*chiA/dM;
				coeffs->f21v5 = (3./4.-3.*eta)*chiA3/dM + (-81./16. +1709.*eta/1008. + 613.*eta2/1008.+(9./4.-3*eta)*chiA2)*chiS + 3./4.*chiS3
				+ (-81./16. - 703.*eta2/112. + 8797.*eta/1008.+(9./4. - 6.*eta)*chiS2)*chiA/dM;
				coeffs->f21v6 = (4163./252.-9287.*eta/1008. - 85.*eta2/112.)*chiA2 + (4163./252. - 2633.*eta/1008. + 461.*eta2/1008.)*chiS2 + (4163./126.-1636.*eta/21. + 1088.*eta2/63.)*chiS*chiA/dM;
			}
			/* End new terms for SEOBNRv4HM */
			break;
		default:
            PRINT_LOG_INFO(LOG_CRITICAL, "wrong SpinAlignedEOBversion value, must be 1, 2 or 3!");
            return CEV_FAILURE;
			break;
		}
	} else {
		coeffs->f21v1 = -3. * chiA / 2.; // Odd modes in dm->0 limit, note that there is a typo in Taracchini et.al. discussion after A7 (even->odd)
		switch (SpinAlignedEOBversion) {
		case 1:
			coeffs->f21v3 = 0.0;
			break;
		case 2:
			coeffs->f21v3 = (chiS * dM * (427. + 79. * eta) + chiA * (147. + 280. * dM * dM + 1251. * eta)) / 84.;
			break;
		case 3:
			coeffs->f21v3 = (chiS * dM * (427. + 79. * eta) + chiA * (147. + 280. * dM * dM + 1251. * eta)) / 84.;
			break;
		case 4:
			coeffs->f21v3 =
			(chiS * dM * (427. + 79. * eta) +
				chiA * (147. + 280. * dM * dM + 1251. * eta)) / 84.;
			/* New terms for SEOBNRv4HM, they are put to zero if use_hm == 0 */
			if (waveform == 0)
			{
				coeffs->f21v4 = 0.0;
				coeffs->f21v5 = 0.0;
				coeffs->f21v6 = 0.0;
			}
			else{
				coeffs->f21v4 = (-6+21*eta/2.)*chiS*chiA;
				coeffs->f21v5 = (3./4.-3.*eta)*chiA3  + (-81./16. - 703.*eta2/112. + 8797.*eta/1008. + (9./4. - 6.*eta)*chiS2)*chiA;
				coeffs->f21v6 = (4163./126.-1636.*eta/21. + 1088.*eta2/63.)*chiS*chiA;
				// printf("f21v1 = %.16f\n f21v3 = %.16f\n f21v4 = %.16f\n f21v5 = %.16f\n f21v6 = %.16f\n", coeffs->f21v1, coeffs->f21v3, coeffs->f21v4, coeffs->f21v5, coeffs->f21v6);
			}
			/* End new terms for SEOBNRv4HM */
			break;
		default:
            PRINT_LOG_INFO(LOG_CRITICAL, "wrong SpinAlignedEOBversion value, must be 1, 2 or 3!");
            return CEV_FAILURE;
			break;
		}
	}

	/*
	 * l = 3, Eqs. A9a - A9c for rho, Eqs. A15b and A15c for f, Eqs. 22 -
	 * 24 of DIN and Eqs. 27c - 27e of PBFRT for delta
	 */
	coeffs->delta33vh3 = 13. / 10.;
	coeffs->delta33vh6 = (-81. * aDelta) / 20. + (39. * CST_PI) / 7.;
	coeffs->delta33vh9 = -227827. / 3000. + (78. * CST_PI * CST_PI) / 7.;
	coeffs->delta33v5 = -80897. * eta / 2430.; 
	if (dM2) {
		coeffs->rho33v2 = -7. / 6. + (2. * eta) / 3.;
		//coeffs->rho33v3 = (chiS * dM * (-4. + 5. * eta) + chiA * (-4. + 19. * eta)) / (6. * dM);
		coeffs->rho33v3 = 0.0;
		coeffs->rho33v4 = -6719. / 3960. + a2 / 2. - (1861. * eta) / 990. + (149. * eta2) / 330.;
		coeffs->rho33v5 = (-4. * a) / 3.;
		coeffs->rho33v6 = 3203101567. / 227026800. + (5. * a2) / 36.;
		coeffs->rho33v6l = -26. / 7.;
		coeffs->rho33v7 = (5297. * a) / 2970. + a * a2 / 3.;
		coeffs->rho33v8 = -57566572157. / 8562153600.;
		coeffs->rho33v8l = 13. / 3.;
		if(waveform == 1){
			//RC: This terms are in Eq.A6 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
			coeffs->rho33v6 = 3203101567. / 227026800. + (5. * a2) / 36. + (-129509./25740. + 41./192. * CST_PI*CST_PI)*eta - 274621./154440.*eta2+ 12011./46332.*eta3;
			coeffs->rho33v10 = -903823148417327./30566888352000.;
			coeffs->rho33v10l = 87347./13860.;
		}

		coeffs->f33v3 = (chiS * dM * (-4. + 5. * eta) + chiA * (-4. + 19. * eta)) / (2. * dM);
		coeffs->f33v4 = 0;
		coeffs->f33v5 = 0;
		coeffs->f33v6 = 0;
		coeffs->f33vh6 = 0;
		if(waveform == 1){
			//RC: This terms are in Eq.A10 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
			coeffs->f33v4 = (3./2. * chiS2 * dM + (3. - 12 * eta) * chiA * chiS + dM * (3./2. -6. * eta) * chiA2)/(dM);
			coeffs->f33v5 = (dM * (241./30. * eta2 + 11./20. * eta + 2./3.) * chiS + (407./30. * eta2 - 593./60. * eta + 2./3.)* chiA)/(dM);
			coeffs->f33v6 = (dM * (6. * eta2 -27. / 2. * eta - 7./ 4.) * chiS2 + (44. * eta2 - 1. * eta - 7./2.) * chiA * chiS + dM * (-12 * eta2 + 11./2. * eta - 7./4.) * chiA2)/dM;
			coeffs->f33vh6 = (dM * (593. / 108. * eta - 81./20.) * chiS + (7339./540. * eta - 81./20.) * chiA)/(dM);
		}
	}
	else {
		coeffs->f33v3 = chiA * 3. / 8.;
		if(waveform == 1){
			//RC: This terms are in Eq.A10 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
			coeffs->f33v4 = ((3. - 12 * eta) * chiA * chiS);
			coeffs->f33v5 = ((407./30. * eta2 - 593./60. * eta + 2./3.)* chiA);
			coeffs->f33v6 = ((44. * eta2 - 1. * eta - 7./2.) * chiA * chiS);
			coeffs->f33vh6 = ((7339./540. * eta - 81./20.) * chiA);
		}
	}

	coeffs->delta32vh3 = (10. + 33. * eta) / (-15. * m1Plus3eta);
	coeffs->delta32vh4 = 4. * aDelta;
	coeffs->delta32vh6 = (-136. * aDelta) / 45. + (52. * CST_PI) / 21.;
	coeffs->delta32vh9 = -9112. / 405. + (208. * CST_PI * CST_PI) / 63.;

	coeffs->rho32v = (4. * chiS * eta) / (-3. * m1Plus3eta);
	/** TODO The term proportional to eta^2 a^2 in coeffs->rho32v2 is wrong, but it was used in the calibration of SEOBNRv2 */
	coeffs->rho32v2 = (-4. * a2 * eta2) / (9. * m1Plus3eta2) + (328. - 1115. * eta
						 + 320. * eta2) / (270. * m1Plus3eta);
	if (SpinAlignedEOBversion == 4) {
		coeffs->rho32v2 = (328. - 1115. * eta +
			320. * eta2) / (270. * m1Plus3eta);
	}
	//coeffs->rho32v3 = (2. * (45. * a * m1Plus3eta3
	//	 - a * eta * (328. - 2099. * eta + 5. * (733. + 20. * a2) * eta2
	//		  - 960. * eta3))) / (405. * m1Plus3eta3);
	coeffs->rho32v3 = 2. / 9. * a;
	coeffs->rho32v4 = a2 / 3. + (-1444528.
			   + 8050045. * eta - 4725605. * eta2 - 20338960. * eta3
			   + 3085640. * eta2 * eta2) / (1603800. * m1Plus3eta2);
	coeffs->rho32v5 = (-2788. * a) / 1215.;
	coeffs->rho32v6 = 5849948554. / 940355325. + (488. * a2) / 405.;
	coeffs->rho32v6l = -104. / 63.;
	coeffs->rho32v8 = -10607269449358. / 3072140846775.;
	coeffs->rho32v8l = 17056. / 8505.;

	if (dM2) {
		coeffs->delta31vh3 = 13. / 30.;
		coeffs->delta31vh6 = (61. * aDelta) / 20. + (13. * CST_PI) / 21.;
		coeffs->delta31vh7 = (-24. * aDelta * aDelta) / 5.;
		coeffs->delta31vh9 = -227827. / 81000. + (26. * CST_PI * CST_PI) / 63.;
		coeffs->delta31v5 = -17. * eta / 10.;

		coeffs->rho31v2 = -13. / 18. - (2. * eta) / 9.;
		//coeffs->rho31v3 = (chiA * (-4. + 11. * eta) + chiS * dM * (-4. + 13. * eta)) / (6. * dM);
		coeffs->rho31v3 = 0.0;
		coeffs->rho31v4 = 101. / 7128.
			- (5. * a2) / 6. - (1685. * eta) / 1782. - (829. * eta2) / 1782.;
		coeffs->rho31v5 = (4. * a) / 9.;
		coeffs->rho31v6 = 11706720301. / 6129723600. - (49. * a2) / 108.;
		coeffs->rho31v6l = -26. / 63.;
		coeffs->rho31v7 = (-2579. * a) / 5346. + a * a2 / 9.;
		coeffs->rho31v8 = 2606097992581. / 4854741091200.;
		coeffs->rho31v8l = 169. / 567.;

		coeffs->f31v3 = (chiA * (-4. + 11. * eta) + chiS * dM * (-4. + 13. * eta)) / (2. * dM);
	} else {
		coeffs->f31v3 = -chiA * 5. / 8.;
	}

	/*
	 * l = 4, Eqs. A10a - A10d for rho, Eq. A15d for f Eqs. 25 - 28 of
	 * DIN and Eqs. 27f - 27i of PBFRT for delta
	 */

	coeffs->delta44vh3 = (112. + 219. * eta) / (-120. * m1Plus3eta);
	coeffs->delta44vh6 = (-464. * aDelta) / 75. + (25136. * CST_PI) / 3465.;
	coeffs->delta44vh9 = 0.;
  	if(waveform == 1){
		//RC: This terms are in Eq.A15 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
  		coeffs->delta44vh9 = -55144./375. + 201088.*CST_PI*CST_PI/10395.;
  	}

	coeffs->rho44v2 = (1614. - 5870. * eta + 2625. * eta2) / (1320. * m1Plus3eta);
	coeffs->rho44v3 = (chiA * (10. - 39. * eta) * dM + chiS * (10. - 41. * eta
					+ 42. * eta2)) / (15. * m1Plus3eta);
	coeffs->rho44v4 = a2 / 2. + (-511573572.
		+ 2338945704. * eta - 313857376. * eta2 - 6733146000. * eta3
		  + 1252563795. * eta2 * eta2) / (317116800. * m1Plus3eta2);
	coeffs->rho44v5 = (-69. * a) / 55.;
	coeffs->rho44v8 = 0.;
 	coeffs->rho44v8l = 0.;
 	coeffs->rho44v10 = 0.;
 	coeffs->rho44v10l = 0;
  	if(waveform == 1){
	  //RC: This terms are in Eq.A8 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
	  coeffs->rho44v4 =
		(-511573572. + 2338945704. * eta - 313857376. * eta2 -
		 6733146000. * eta3 +
		 1252563795. * eta2 * eta2) / (317116800. * m1Plus3eta2)
		+ chiS2/2. + dM*chiS*chiA + dM2*chiA2/2.;
	  coeffs->rho44v5 = chiA*dM*(-8280. + 42716.*eta - 57990.*eta2 + 8955*eta3)/(6600.*m1Plus3eta2)
		+ chiS*(-8280. + 66284.*eta-176418.*eta2+128085.*eta3 + 88650*eta2*eta2)/(6600.*m1Plus3eta2);
	  coeffs->rho44v8 = -172066910136202271./19426955708160000.;
	  coeffs->rho44v8l = 845198./190575.;
	  coeffs->rho44v10 = - 17154485653213713419357./568432724020761600000.;
	  coeffs->rho44v10l = 22324502267./3815311500;
	}
	coeffs->rho44v6 = 16600939332793. / 1098809712000. + (217. * a2) / 3960.;
	coeffs->rho44v6l = -12568. / 3465.;

	if (dM2) {
		coeffs->delta43vh3 = (486. + 4961. * eta) / (810. * (1. - 2. * eta));
		coeffs->delta43vh4 = (11. * aDelta) / 4.;
		coeffs->delta43vh6 = 1571. * CST_PI / 385.;

		//coeffs->rho43v = (5. * (chiA - chiS * dM) * eta) / (8. * dM * (-1. + 2. * eta));
		coeffs->rho43v = 0.0;
		coeffs->rho43v2 = (222. - 547. * eta + 160. * eta2) / (176. * (-1. + 2. * eta));
		coeffs->rho43v4 = -6894273. / 7047040. + (3. * a2) / 8.;
		coeffs->rho43v5 = (-12113. * a) / 6160.;
		coeffs->rho43v6 = 1664224207351. / 195343948800.;
		coeffs->rho43v6l = -1571. / 770.;

		coeffs->f43v = (5. * (chiA - chiS * dM) * eta) / (2. * dM * (-1. + 2. * eta));
	} else {
		coeffs->f43v = -5. * chiA / 4.;
	}

	coeffs->delta42vh3 = (7. * (1. + 6. * eta)) / (-15. * m1Plus3eta);
	coeffs->delta42vh6 = (212. * aDelta) / 75. + (6284. * CST_PI) / 3465.;

	coeffs->rho42v2 = (1146. - 3530. * eta + 285. * eta2) / (1320. * m1Plus3eta);
	coeffs->rho42v3 = (chiA * (10. - 21. * eta) * dM + chiS * (10. - 59. * eta
					+ 78. * eta2)) / (15. * m1Plus3eta);
	coeffs->rho42v4 = a2 / 2. + (-114859044. + 295834536. * eta + 1204388696. * eta2 - 3047981160. * eta3
		   - 379526805. * eta2 * eta2) / (317116800. * m1Plus3eta2);
	coeffs->rho42v5 = (-7. * a) / 110.;
	coeffs->rho42v6 = 848238724511. / 219761942400. + (2323. * a2) / 3960.;
	coeffs->rho42v6l = -3142. / 3465.;

	if (dM2) {
		coeffs->delta41vh3 = (2. + 507. * eta) / (10. * (1. - 2. * eta));
		coeffs->delta41vh4 = (11. * aDelta) / 12.;
		coeffs->delta41vh6 = 1571. * CST_PI / 3465.;

		//coeffs->rho41v = (5. * (chiA - chiS * dM) * eta) / (8. * dM * (-1. + 2. * eta));
		coeffs->rho41v = 0.0;
		coeffs->rho41v2 = (602. - 1385. * eta + 288. * eta2) / (528. * (-1. + 2. * eta));
		coeffs->rho41v4 = -7775491. / 21141120. + (3. * a2) / 8.;
		coeffs->rho41v5 = (-20033. * a) / 55440. - (5 * a * a2) / 6.;
		coeffs->rho41v6 = 1227423222031. / 1758095539200.;
		coeffs->rho41v6l = -1571. / 6930.;

		coeffs->f41v = (5. * (chiA - chiS * dM) * eta) / (2. * dM * (-1. + 2. * eta));
	} else {
		coeffs->f41v = -5. * chiA / 4.;
	}

	/*
	 * l = 5, Eqs. A11a - A11e for rho, Eq. 29 of DIN and Eqs. E1a and
	 * E1b of PBFRT for delta
	 */
	coeffs->delta55vh3 =
		(96875. + 857528. * eta) / (131250. * (1 - 2 * eta));
	coeffs->delta55vh6 = 0;
	coeffs->delta55vh9 = 0;
 	if(waveform == 1){
		//RC: This terms are in Eq.A16 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
		coeffs->delta55vh6 = 3865./429.*CST_PI;
		coeffs->delta55vh9 = (-7686949127. + 954500400.*CST_PI*CST_PI)/31783752.;
 	}
	if (dM2) {

		coeffs->rho55v2 = (487. - 1298. * eta + 512. * eta2) / (390. * (-1. + 2. * eta));
		coeffs->rho55v3 = (-2. * a) / 3.;
		coeffs->rho55v4 = -3353747. / 2129400. + a2 / 2.;
		coeffs->rho55v5 = -241. * a / 195.;
		coeffs->rho55v6 = 0.;
		coeffs->rho55v6l = 0.;
		coeffs->rho55v8 = 0.;
		coeffs->rho55v8l = 0.;
		coeffs->rho55v10 = 0.;
		coeffs->rho55v10l = 0.;
		coeffs->f55v3 = 0.;
		coeffs->f55v4 = 0.;
		coeffs->f55v5c = 0;
		if(waveform == 1){
			//RC: This terms are in Eq.A9 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
			coeffs->rho55v6 = 190606537999247./11957879934000.;
			coeffs->rho55v6l = - 1546./429.;
			coeffs->rho55v8 = - 1213641959949291437./118143853747920000.;
			coeffs->rho55v8l = 376451./83655.;
			coeffs->rho55v10 = -150082616449726042201261./4837990810977324000000.;
			coeffs->rho55v10l = 2592446431./456756300.;

			coeffs->f55v3 = chiA/dM *(10./(3.*(-1.+2.*eta)) - 70.*eta/(3.*(-1.+2.*eta)) + 110.*eta2/(3.*(-1.+2.*eta)) ) +
				chiS*(10./(3.*(-1.+2.*eta)) -10.*eta/(-1.+2.*eta) + 10*eta2/(-1.+2.*eta));
			coeffs->f55v4 = chiS2*(-5./(2.*(-1.+2.*eta)) + 5.*eta/(-1.+2.*eta)) +
				chiA*chiS/dM *(-5./(-1.+2.*eta) + 30.*eta/(-1.+2.*eta) - 40.*eta2/(-1.+2.*eta)) +
				chiA2*(-5./(2.*(-1.+2.*eta)) + 15.*eta/(-1.+2.*eta) - 20.*eta2/(-1.+2.*eta));
			coeffs->f55v5c = 0; //RC: this is the calibration parameter which is initially set to 0.
	  }
	}
	else{
		coeffs->f55v3 = 0;
		coeffs->f55v4 = 0;
		coeffs->f55v5c = 0;
		if(waveform == 1){
			//RC: This terms are in Eq.A12 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
			coeffs->f55v3 = chiA *(10./(3.*(-1.+2.*eta)) - 70.*eta/(3.*(-1.+2.*eta)) + 110.*eta2/(3.*(-1.+2.*eta)) );
			coeffs->f55v4 = chiA*chiS *(-5./(-1.+2.*eta) + 30.*eta/(-1.+2.*eta) - 40.*eta2/(-1.+2.*eta));
			coeffs->f55v5c = 0;
		}
	}
	coeffs->delta54vh3 = 8. / 15.;
	coeffs->delta54vh4 = 12. * aDelta / 5.;

	coeffs->rho54v2 = (-17448. + 96019. * eta - 127610. * eta2
		  + 33320. * eta3) / (13650. * (1. - 5. * eta + 5. * eta2));
	coeffs->rho54v3 = (-2. * a) / 15.;
	coeffs->rho54v4 = -16213384. / 15526875. + (2. * a2) / 5.;

	if (dM2) {
		coeffs->delta53vh3 = 31. / 70.;

		coeffs->rho53v2 = (375. - 850. * eta + 176. * eta2) / (390. * (-1. + 2. * eta));
		coeffs->rho53v3 = (-2. * a) / 3.;
		coeffs->rho53v4 = -410833. / 709800. + a2 / 2.;
		coeffs->rho53v5 = -103. * a / 325.;
	}
	coeffs->delta52vh3 = 4. / 15.;
	coeffs->delta52vh4 = 6. * aDelta / 5.;

	coeffs->rho52v2 = (-15828. + 84679. * eta - 104930. * eta2
		  + 21980. * eta3) / (13650. * (1. - 5. * eta + 5. * eta2));
	coeffs->rho52v3 = (-2. * a) / 15.;
	coeffs->rho52v4 = -7187914. / 15526875. + (2. * a2) / 5.;

	if (dM2) {
		coeffs->delta51vh3 = 31. / 210.;

		coeffs->rho51v2 = (319. - 626. * eta + 8. * eta2) / (390. * (-1. + 2. * eta));
		coeffs->rho51v3 = (-2. * a) / 3.;
		coeffs->rho51v4 = -31877. / 304200. + a2 / 2.;
		coeffs->rho51v5 = 139. * a / 975.;
	}
	/*
	 * l = 6, Eqs. A12a - A12f for rho, Eqs. E1c and E1d of PBFRT for
	 * delta
	 */

	coeffs->delta66vh3 = 43. / 70.;

	coeffs->rho66v2 = (-106. + 602. * eta - 861. * eta2
			   + 273. * eta3) / (84. * (1. - 5. * eta + 5. * eta2));
	coeffs->rho66v3 = (-2. * a) / 3.;
	coeffs->rho66v4 = -1025435. / 659736. + a2 / 2.;

	if (dM2) {
		coeffs->delta65vh3 = 10. / 21.;

		coeffs->rho65v2 = (-185. + 838. * eta - 910. * eta2
				+ 220. * eta3) / (144. * (dM2 + 3. * eta2));
		coeffs->rho65v3 = -2. * a / 9.;
	}
	coeffs->delta64vh3 = 43. / 105.;

	coeffs->rho64v2 = (-86. + 462. * eta - 581. * eta2
			   + 133. * eta3) / (84. * (1. - 5. * eta + 5. * eta2));
	coeffs->rho64v3 = (-2. * a) / 3.;
	coeffs->rho64v4 = -476887. / 659736. + a2 / 2.;

	if (dM2) {
		coeffs->delta63vh3 = 2. / 7.;

		coeffs->rho63v2 = (-169. + 742. * eta - 750. * eta2
				+ 156. * eta3) / (144. * (dM2 + 3. * eta2));
		coeffs->rho63v3 = -2. * a / 9.;
	}
	coeffs->delta62vh3 = 43. / 210.;

	coeffs->rho62v2 = (-74. + 378. * eta - 413. * eta2
			+ 49. * eta3) / (84. * (1. - 5. * eta + 5. * eta2));
	coeffs->rho62v3 = (-2. * a) / 3.;
	coeffs->rho62v4 = -817991. / 3298680. + a2 / 2.;

	if (dM2) {
		coeffs->delta61vh3 = 2. / 21.;

		coeffs->rho61v2 = (-161. + 694. * eta - 670. * eta2
				+ 124. * eta3) / (144. * (dM2 + 3. * eta2));
		coeffs->rho61v3 = -2. * a / 9.;
	}
	/*
	 * l = 7, Eqs. A13a - A13g for rho, Eqs. E1e and E1f of PBFRT for
	 * delta
	 */
	if (dM2) {
		coeffs->delta77vh3 = 19. / 36.;

		coeffs->rho77v2 = (-906. + 4246. * eta - 4963. * eta2
				   + 1380. * eta3) / (714. * (dM2 + 3. * eta2));
		coeffs->rho77v3 = -2. * a / 3.;
	}
	coeffs->rho76v2 = (2144. - 16185. * eta + 37828. * eta2 - 29351. * eta3
		 + 6104. * eta2 * eta2) / (1666. * (-1 + 7 * eta - 14 * eta2
							+ 7 * eta3));

	if (dM2) {
		coeffs->delta75vh3 = 95. / 252.;

		coeffs->rho75v2 = (-762. + 3382. * eta - 3523. * eta2
				+ 804. * eta3) / (714. * (dM2 + 3. * eta2));
		coeffs->rho75v3 = -2. * a / 3.;
	}
	coeffs->rho74v2 = (17756. - 131805. * eta + 298872. * eta2 - 217959. * eta3
		+ 41076. * eta2 * eta2) / (14994. * (-1. + 7. * eta - 14. * eta2
						 + 7. * eta3));

	if (dM2) {
		coeffs->delta73vh3 = 19. / 84.;

		coeffs->rho73v2 = (-666. + 2806. * eta - 2563. * eta2
				+ 420. * eta3) / (714. * (dM2 + 3. * eta2));
		coeffs->rho73v3 = -2. * a / 3.;
	}
	coeffs->rho72v2 = (16832. - 123489. * eta + 273924. * eta2 - 190239. * eta3
		+ 32760. * eta2 * eta2) / (14994. * (-1. + 7. * eta - 14. * eta2
						 + 7. * eta3));

	if (dM2) {
		coeffs->delta71vh3 = 19. / 252.;

		coeffs->rho71v2 = (-618. + 2518. * eta - 2083. * eta2
				+ 228. * eta3) / (714. * (dM2 + 3. * eta2));
		coeffs->rho71v3 = -2. * a / 3.;
	}
	/* l = 8, Eqs. A14a - A14h */

	coeffs->rho88v2 = (3482. - 26778. * eta + 64659. * eta2 - 53445. * eta3
		 + 12243. * eta2 * eta2) / (2736. * (-1. + 7. * eta - 14. * eta2
						 + 7. * eta3));

	if (dM2) {
		coeffs->rho87v2 = (23478. - 154099. * eta + 309498. * eta2 - 207550. * eta3
		+ 38920 * eta2 * eta2) / (18240. * (-1 + 6 * eta - 10 * eta2
							+ 4 * eta3));
	}
	coeffs->rho86v2 = (1002. - 7498. * eta + 17269. * eta2 - 13055. * eta3
		   + 2653. * eta2 * eta2) / (912. * (-1. + 7. * eta - 14. * eta2
						 + 7. * eta3));

	if (dM2) {
		coeffs->rho85v2 = (4350. - 28055. * eta + 54642. * eta2 - 34598. * eta3
				   + 6056. * eta2 * eta2) / (3648. * (-1. + 6. * eta - 10. * eta2
								  + 4. * eta3));
	}
	coeffs->rho84v2 = (2666. - 19434. * eta + 42627. * eta2 - 28965. * eta3
		  + 4899. * eta2 * eta2) / (2736. * (-1. + 7. * eta - 14. * eta2
						 + 7. * eta3));

	if (dM2) {
		coeffs->rho83v2 = (20598. - 131059. * eta + 249018. * eta2 - 149950. * eta3
				   + 24520. * eta2 * eta2) / (18240. * (-1. + 6. * eta - 10. * eta2
								  + 4. * eta3));
	}
	coeffs->rho82v2 = (2462. - 17598. * eta + 37119. * eta2 - 22845. * eta3
		  + 3063. * eta2 * eta2) / (2736. * (-1. + 7. * eta - 14. * eta2
						 + 7. * eta3));

	if (dM2) {
		coeffs->rho81v2 = (20022. - 126451. * eta + 236922. * eta2 - 138430. * eta3
				   + 21640. * eta2 * eta2) / (18240. * (-1. + 6. * eta - 10. * eta2
								  + 4. * eta3));
	}

    /* New Factorized waveform coeffs */
    coeffs->h22T0ff00 = 1./2.;
    coeffs->h22T0ff02 = 1./2.;
    coeffs->h22T0ff11 = 1;
    coeffs->h22T0ff20 = -1./2.;

    coeffs->h22T2ff40 = -1./14.-15.*eta/28.;
    coeffs->h22T2ff31 = (1./7. + 15.*eta/14.);
    coeffs->h22T2ff20 = -9./7. + 31.*eta/28.;
    coeffs->h22T2ff13 = (1./7. + 15.*eta/14.);
    coeffs->h22T2ff11 = (25./42. - 16.*eta/7.);
    coeffs->h22T2ff04 = 1./14. + 15.*eta/28.;
    coeffs->h22T2ff02 = -13./21. + 3.*eta/28.;
    coeffs->h22T2ff00 = -3./2. + eta/2. + 0.5*I*(cax+cay + dM*(csx+csy));

    coeffs->h22T3ff30 = -4.*dM*caz/3. - (4./3. - eta) * csz;
    coeffs->h22T3ff21 = -8.*dM*caz/3. - (8./3. - 2.*eta) * csz;
    coeffs->h22T3ff12 = 4.*dM*caz/3. + (4./3. - eta) * csz;
    coeffs->h22T3ff10 = (1.*dM/3.)*caz + (1./3.+2.*eta/3.)*csz;
    // coeffs->h22T3ff03 = (4.*dM/3.)*caz + (4./3.-eta)*csz;
    coeffs->h22T3ff01 = -dM*caz + (-1. + 7.*eta/6.)*csz;

    coeffs->h22T4ff60 = -1./42. - 179.*eta/336. - 25.*eta2/48.;
    coeffs->h22T4ff51 = (1./21. + 179.*eta/168. + 25.*eta2/24.);
    coeffs->h22T4ff42 = -1./42. - 179.*eta/336. -25.*eta2/48.;
    coeffs->h22T4ff40 = -11./21. - 923.*eta/336. + 335.*eta2/336.;
    coeffs->h22T4ff33 = (2./21. + 179.*eta/84. + 25.*eta2/12.);
    coeffs->h22T4ff31 = (71./63. + 935.*eta/252. - 131.*eta2/63.);
    coeffs->h22T4ff24 = 1./42. + 179.*eta/336. + 25.*eta2/48.;
    coeffs->h22T4ff22 = 23./84. - 45.*eta/56. - 109.*eta2/168.;
    coeffs->h22T4ff20 = -425./504. + 1355.*eta/252. - 20.*eta2/63.+ cax*I*(-3./2.-19.*eta/28.) + cay*(-3./2.-37.*eta/42.) + csx*dM*I*(-3./2.-9.*eta/28.) + 
        csy*dM*(-3./2.-11.*eta/14.) + 0.5*(pow(dM*cay+(1.-2.*eta)*csy, 2.) + pow(dM*caz+(1.-2.*eta)*csy, 2.));
    coeffs->h22T4ff15 = (1./21. + 179.*eta/168. + 25.*eta2/24.);
    coeffs->h22T4ff13 = (-1./168. + 29.*eta/12. - 58.*eta2/21.);
    coeffs->h22T4ff11 = I*(-359./378. - 592.*eta/189. - 307.*eta2/378.) + cax*(-3./2.-51.*eta/28.) + cay*I*(5./2.-5.*eta/12.) + csx*dM*(-3./2.+125.*eta/84.) + 
        csy*I*dM*(5./2.+eta/12.) + (dM*cax + (1.-2.*eta)*csx)*(dM*cay + (1.-2.*eta)*csy);
    coeffs->h22T4ff06 = 1./42. + 179.*eta/336. + 25.*eta2/48.;
    coeffs->h22T4ff04 = 317./1008. - 1045.*eta/1008. - 101.*eta2/144.;
    coeffs->h22T4ff02 = -10379./3024. - 635.*eta/189. + 26.*eta2/189. + cax*I*(9./4.-89.*eta/14.) + cay*(3./4.+71.*eta/84.) + csx*I*dM*(9./4.-55.*eta/42.) + 
        csy*dM*(3./4.+13.*eta/84.) - I*(dM*cax + (1.-2.*eta)*csx) * (dM*cay + (1.-2.*eta)*csy) + 0.5*(pow(dM*cay+(1.-2.*eta)*csy, 2.) + pow(dM*caz+(1.-2.*eta)*csy, 2.));
    coeffs->h22T4ff00 = 65./252. - 31.*eta/36. + 205.*eta2/252. + cax*I*(-1.+19.*eta/14.) + cay*(-1./3.-19.*eta/21.) + csx*I*dM*(-1.+9.*eta/14.) + 
        csy*dM*(-1./3.+5.*eta/21.) -2.*I*(dM*cax+(1.+eta)*csx)*(dM*cay+(1.+eta)*csy) + 6.*I*csx*csy + 
        0.5*(pow(dM*cay+(1.+eta)*csy, 2.) + pow(dM*caz+(1.+eta)*csy, 2.) - 3.*eta2 * (csy*csy + csz*csz));

    // h21
    if (dM2)
    {
        coeffs->h21T0ff01 = 1.;
        coeffs->h21T1ff00 = -3.*caz / 2. / dM - 3.*csz / 2.;

        coeffs->h21T2ff21 = -20./7. - 11.*eta/14.;
        coeffs->h21T2ff12 = 83./14. - 6.*eta/7.;
        coeffs->h21T2ff03 = 5./28. - 5.*eta/14.;
        coeffs->h21T2ff01 = -16./7. + 11.*eta/7. -6.*I*cax - I*(6.+eta) * csx / dM + 4.*eta * csy/dM;
        coeffs->h21T2ff10 = -3.*cay + (-3. + 8.*eta) * csy/dM;
        coeffs->h21T3ff20 = (9./2.+95.*eta/28.)*caz/dM + (9./2.+87.*eta/28.)*csz;
        coeffs->h21T3ff11 = (-55./6.+44.*eta/21.)*caz/dM + (-55./6.-3.*eta/7.)*csz;
        coeffs->h21T3ff02 = (-9./4.+327.*eta/28.)*caz/dM + (-9./4.+107.*eta/28.)*csz;
        coeffs->h21T3ff00 = (6.-95.*eta/14.)*caz/dM + (6.-59.*eta/14.)*csz + I*(9.-36.*eta)*cax*caz/dM + 9.*I*(caz*csx+cax*csz) + 9.*I*csx*csz/dM;
    }
    else
    {
        coeffs->h21T1ff00 = -3.*caz/2.;
        coeffs->h21T2ff10 = (-3. + 8.*eta) * csy;
        coeffs->h21T2ff01 = -I*(6.+eta)*csx + 4.*eta*csy;
        coeffs->h21T3ff20 = (9./2.+95.*eta/28.)*caz;
        coeffs->h21T3ff11 = (-55./6.+44.*eta/21.)*caz;
        coeffs->h21T3ff02 = (-9./4.+327.*eta/28.)*caz;
        coeffs->h21T3ff00 = (6.-95.*eta/14.)*caz + I*(9.-36.*eta)*cax*caz + 9.*I*csx*csz;
    }

    // h33
    if (dM2)
    {
        coeffs->h33T0ff30 = -2./9.;
        coeffs->h33T0ff21 = -2./3.;
        coeffs->h33T0ff12 = 2./3.;
        coeffs->h33T0ff10 = 4./9.;
        coeffs->h33T0ff03 = 2./9.;
        coeffs->h33T0ff01 = 7./9.;

        coeffs->h33T2ff50 = -2./27.-8.*eta/27.;
        coeffs->h33T2ff41 = -2./9. - 8.*eta/9.;
        coeffs->h33T2ff32 = 4./27. + 16.*eta/27.;
        coeffs->h33T2ff30 = -28./27. + 20.*eta/27.;
        coeffs->h33T2ff23 = -4./27. - 16.*eta/27.;
        coeffs->h33T2ff21 = -20./9. + 19.*eta/9.;
        coeffs->h33T2ff14 = 2./9. + 8.*eta/9.;
        coeffs->h33T2ff12 = 2./3.-2.*eta;
        coeffs->h33T2ff10 = -37./81.-4.*eta/81. + (-I)*(4.*eta/9.)*(csx-csy)/dM;
        coeffs->h33T2ff05 = 2./27. + 8.*eta/27.;
        coeffs->h33T2ff03 = -11./18.+2.*eta/3.;
        coeffs->h33T2ff01 = -80./27 + 22.*eta/27. - (16.*eta*I/9.)*(csx+csy)/dM;

        coeffs->h33T3ff40 = (8./9.-32.*eta/9.)*caz/dM + (8./9.-2.*eta/3.)*csz;
        coeffs->h33T3ff31 = (-8./3.+32.*eta/3.)*caz/dM + (-8./3.+2.*eta)*csz;
        coeffs->h33T3ff22 = (-8./3.+32.*eta/3.)*caz/dM + (-8./3.+2.*eta)*csz;
        coeffs->h33T3ff13 = (8./9.-32.*eta/9.)*caz/dM + (8./9.-2.*eta/3.)*csz;
        coeffs->h33T3ff20 = (-10./9. + 13.*eta/3.) * caz / dM + (-10./9. + eta/3.)*csz;
        coeffs->h33T3ff11 = (10./9.-35.*eta/9.)*caz/dM + (-10./9. + 10.*eta/9.)*csz;
        coeffs->h33T3ff02 = (-4./3. + 119.*eta/18.)*caz / dM + (-4./3. + 35.*eta/18.)*csz;
        coeffs->h33T3ff00 = (-2./9. + 10.*eta/9.)*caz/dM + (-2./9.+eta/3.)*csz;
    }
    else
    {
        coeffs->h33T2ff10 = (-I)*(4.*eta/9.)*(csx-csy);
        coeffs->h33T2ff01 = - (16.*eta*I/9.)*(csx+csy);

        coeffs->h33T3ff40 = (8./9.-32.*eta/9.)*caz;
        coeffs->h33T3ff31 = (-8./3.+32.*eta/3.)*caz;
        coeffs->h33T3ff22 = (-8./3.+32.*eta/3.)*caz;
        coeffs->h33T3ff13 = (8./9.-32.*eta/9.)*caz;
        coeffs->h33T3ff20 = (-10./9. + 13.*eta/3.) * caz;
        coeffs->h33T3ff11 = (10./9.-35.*eta/9.)*caz;
        coeffs->h33T3ff02 = (-4./3. + 119.*eta/18.)*caz;
        coeffs->h33T3ff00 = (-2./9. + 10.*eta/9.)*caz;
    }

    // h32
    coeffs->h32T0ff11 = 1./4.;
    coeffs->h32T0ff02 = 1.;
    coeffs->h32T1ff10 = -eta*chiS / m1Plus3eta;
    coeffs->h32T1ff01 = -4.*eta*chiS / m1Plus3eta;
    coeffs->h32T2ff31 = (16. - eta*(47. + 10.*eta)) / (24.) / m1Plus3eta;
    coeffs->h32T2ff22 = (28. - 5.*eta*(16.+5.*eta)) / 10. / m1Plus3eta;
    coeffs->h32T2ff13 = (-31.+eta*(107.-26.*eta)) / 8. / m1Plus3eta;
    coeffs->h32T2ff11 = (-73. + eta*(239.-26.*eta)) / 36. / m1Plus3eta;
    coeffs->h32T2ff04 = (43. + 5.*eta*(-37.+eta)) / 60. / m1Plus3eta;
    coeffs->h32T2ff02 = (527. + 5.*eta*(-347.+161.*eta)) / 180. / m1Plus3eta;

    // h31
    if (dM2)
    {
        coeffs->h31T0ff30 = -6.;
        coeffs->h31T0ff21 = -6.;
        coeffs->h31T0ff12 = -6.;
        coeffs->h31T0ff10 = 12.;
        coeffs->h31T0ff03 = -6.;
        coeffs->h31T0ff01 = 7.;

        coeffs->h31T2ff50 = -2.-8.*eta;
        coeffs->h31T2ff41 = -2.-8.*eta;
        coeffs->h31T2ff32 = -4.-16.*eta;
        coeffs->h31T2ff30 = -28.+20.*eta;
        coeffs->h31T2ff23 = -4.-16.*eta;
        coeffs->h31T2ff21 = -20.+19.*eta;
        coeffs->h31T2ff14 = -2.-8.*eta;
        coeffs->h31T2ff12 = -6.+26.*eta;
        coeffs->h31T2ff10 = -37./3. - 4.*eta/3.;
        coeffs->h31T2ff05 = -2.-8.*eta;
        coeffs->h31T2ff03 = 53./2.+6.*eta;
        coeffs->h31T2ff01 = -80./3. + 22.*eta/3.;

        coeffs->h31T3ff20 = (6.-25.*eta)*chiA / dM + (6.-9.*eta)*chiS;
        coeffs->h31T3ff11 = (-6.+31.*eta)*chiA / dM - (6.+11.*eta)*chiS;
        coeffs->h31T3ff02 = (-12.+87.*eta/2.)*chiA/dM + (-12.+19.*eta/2)*chiS;
        coeffs->h31T3ff00 = (-2.+10.*eta)*chiA/dM + (-2.+3.*eta)*chiS;
    }
    else
    {
        coeffs->h31T3ff20 = (6.-25.*eta)*chiA;
        coeffs->h31T3ff11 = (-6.+31.*eta)*chiA;
        coeffs->h31T3ff02 = (-12.+87.*eta/2.)*chiA;
        coeffs->h31T3ff00 = (-2.+10.*eta)*chiA;
    }

    // h44
    coeffs->h44T0ff40 = (3./32.);
    coeffs->h44T0ff31 = (-3./8.);
    coeffs->h44T0ff22 = (-9./16.);
    coeffs->h44T0ff20 = (-9./32.);
    coeffs->h44T0ff13 = (3./8.);
    coeffs->h44T0ff11 = (27./32.);
    coeffs->h44T0ff04 = (3./32.);
    coeffs->h44T0ff02 = (51./64.); 
    coeffs->h44T0ff00 = 7./64.;    
    
    REAL8 m2Plus35eta = -2. + 35.*eta;
    coeffs->h44T2ff60 = 9.*(-4. + eta*m2Plus35eta)/704./m1Plus3eta;
    coeffs->h44T2ff51 =  -9.*(-4. + eta*m2Plus35eta)/176./m1Plus3eta;
    coeffs->h44T2ff42 = 45. * (4. - eta*m2Plus35eta)/704./m1Plus3eta;
    coeffs->h44T2ff40 = -9.*(46. + eta*(-186 + 109*eta))/704./m1Plus3eta;
    coeffs->h44T2ff31 = 3. * (2348. + 25*eta*(-386+243*eta))/3520./m1Plus3eta;
    coeffs->h44T2ff24 = 45. * ( 4. - eta * m2Plus35eta) / 704./m1Plus3eta;
    coeffs->h44T2ff22 = 3.*(5426. + 5.*eta*(-4546. + 3237.*eta))/7040./m1Plus3eta;
    coeffs->h44T2ff20 = (579. + eta*(-2290.+783.*eta) + 132*eta*(I*cax+cay-I*dM*(csx-I*csy)) )/1408./m1Plus3eta;
    coeffs->h44T2ff15 = 9.*(-4.+eta*m2Plus35eta)/176./m1Plus3eta;
    coeffs->h44T2ff13 = -3.*(577.+15.*eta*(-196+295*eta))/3520./m1Plus3eta;
    coeffs->h44T2ff11 = (1416. - 5.*eta*(536.+465.*eta) - 825*eta*(I*cax+cay-dM*I*(csx-I*csy)) )/1760./m1Plus3eta;
    coeffs->h44T2ff06 = 9.*(-4.+eta*m2Plus35eta)/704./m1Plus3eta;
    coeffs->h44T2ff04 = 9. *(859.+10.*eta*(-512.+457.*eta))/14080./m1Plus3eta;
    coeffs->h44T2ff02 = (53913. + 30.*eta*(-5746. + 1291.*eta) - 16180*eta*(I*cax+cay-dM*(I*csx+csy)))/14080./m1Plus3eta;
    coeffs->h44T2ff00 = (397.+eta*(-1358.+621.*eta) - 132*eta*(I*cax+cay-dM*(I*csx+csy)) )/704./m1Plus3eta;

    // h43
    if (dM2)
    {
        coeffs->h43T0ff21 = -2./27.;
        coeffs->h43T0ff12 = 10./27.;
        coeffs->h43T0ff03 = 23. / 27;
        coeffs->h43T0ff01 = 4./27.;
        coeffs->h43T1ff20 = (chiA/dM - chiS) * 5. * eta / 27. / (1. - 2. * eta);
        coeffs->h43T1ff11 = (-chiA/dM + chiS) * 25. * eta / 27. / (1. - 2. * eta);
        coeffs->h43T1ff02 = (-chiA/dM + chiS) * 115. * eta / 54. / (1. - 2. * eta);
        coeffs->h43T1ff00 = (-chiA/dM + chiS) * 10. * eta / 27. / (1. - 2. * eta);
    }
    else
    {
        coeffs->h43T1ff20 = 5.*chiA/54.;
        coeffs->h43T1ff11 = -25.*chiA/54.;
        coeffs->h43T1ff02 = -115.*chiA/108.;
        coeffs->h43T1ff00 = -5.*chiA/27.;
    }


    // h42
    coeffs->h42T0ff40 = 3./2.;
    coeffs->h42T0ff31 = -3.;
    coeffs->h42T0ff20 = -9./2.;
    coeffs->h42T0ff13 = -3.;
    coeffs->h42T0ff11 = 27./4.;
    coeffs->h42T0ff04 = -3./2.;
    coeffs->h42T0ff02 = 3./4.;
    coeffs->h42T0ff00 = 7./4.;

    coeffs->h42T2ff60 = 9.*(-4. + eta*m2Plus35eta) / 44. / m1Plus3eta;
    coeffs->h42T2ff51 = 9.*(-4. + eta*m2Plus35eta) / (-22.) / m1Plus3eta;
    coeffs->h42T2ff42 = coeffs->h42T2ff60;
    coeffs->h42T2ff40 = -9.*(46.+eta*(-186+109*eta)) / 44. / m1Plus3eta;
    coeffs->h42T2ff33 = 9.*(-4. + eta*m2Plus35eta) / (-11.) / m1Plus3eta;
    coeffs->h42T2ff31 = 3.*(2348.+25.*eta*(-386.+243.*eta))/440./m1Plus3eta;
    coeffs->h42T2ff24 = -coeffs->h42T2ff60;
    coeffs->h42T2ff22 = (822.-15.*eta*(54.+199.*eta)) / 440. / m1Plus3eta;
    coeffs->h42T2ff20 = (579.+eta*(-2290.+783.*eta)) / 88. / m1Plus3eta;
    coeffs->h42T2ff15 = coeffs->h42T2ff51;
    coeffs->h42T2ff13 = 3.*(523+5.*eta*(-952.+1535*eta)) / 440. / m1Plus3eta;
    coeffs->h42T2ff11 = (1416. - 5.*eta*(536.+465.*eta)) / 220. / m1Plus3eta;
    coeffs->h42T2ff06 = -coeffs->h42T2ff60;
    coeffs->h42T2ff04 = (-8361.+30.*eta*(838.-35.*eta)) / 880. / m1Plus3eta;
    coeffs->h42T2ff02 = 3.*(919.+10.*eta*(-344.+113.*eta)) / 880. / m1Plus3eta;
    coeffs->h42T2ff00 = (397.+eta*(-1358.+621.*eta)) / 44. / m1Plus3eta;

    // h41
    if (dM2)
    {
        coeffs->h41T0ff21 = -6.;
        coeffs->h41T0ff12 = 10.;
        coeffs->h41T0ff03 = -11.;
        coeffs->h41T0ff01 = 12.;
        coeffs->h41T1ff20 = (chiA/dM - chiS) * 15. * eta / (1.-2.*eta);
        coeffs->h41T1ff11 = (-chiA/dM + chiS) * 25. * eta / (1.-2.*eta);
        coeffs->h41T1ff02 = (chiA/dM - chiS) * 55. * eta / 2. / (1.-2.*eta);
        coeffs->h41T1ff00 = (-chiA/dM + chiS) * 30. * eta / (1-2.*eta);
    }
    else
    {
        coeffs->h41T1ff20 = 15.*chiA/2.;
        coeffs->h41T1ff11 = -25.*chiA/2.;
        coeffs->h41T1ff02 = 55.*chiA/4.;
        coeffs->h41T1ff00 = -15.*chiA;
    }
	/* All relevant coefficients should be set, so we return */
	return CEV_SUCCESS;
}


INT CalculateFacWaveformAmpResV3(FacWaveformCoeffs *const hCoeffs,
                                 const REAL8 eta,
                                 const INT l,
                                 const INT m,
                                 const REAL8 v,
                                 const REAL8 vPhi,
                                 const REAL8 vh3,
                                 const REAL8 eulerlogxabs,
                                 const REAL8 x0,
                                 const REAL8 x1,
                                 const REAL8 x2,
                                 COMPLEX16 *facAmpRes)
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
    v2 = v*v;
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
    g2 = pow(f2, -4./3.);
    switch(l)
    {
        case 2:
            switch(ABS(m))
            {
                case 2:
                    rhoHigh =
                        v2*v2*v*(hCoeffs->rho22v5 +
                        v * (hCoeffs->rho22v6 + hCoeffs->rho22v6l * eulerlogxabs +
                        v * (hCoeffs->rho22v7 + 
                        v * (hCoeffs->rho22v8 + hCoeffs->rho22v8l * eulerlogxabs +
                        v2 * (hCoeffs->rho22v10 + hCoeffs->rho22v10l * eulerlogxabs)))));
                    rhoLow = 1. + v2 * (hCoeffs->rho22v2 +
                        v * (hCoeffs->rho22v3 +
                        v * (hCoeffs->rho22v4)));

                    A0lm = hCoeffs->h22T0ff00 * ff00 + I*hCoeffs->h22T0ff11 * ff11 + hCoeffs->h22T0ff02 * ff02 + hCoeffs->h22T0ff20 * ff20;
                    A2lm = hCoeffs->h22T2ff00 * ff00 + I*hCoeffs->h22T2ff11 * ff11 + hCoeffs->h22T2ff02 * ff02 + hCoeffs->h22T2ff20 * ff20 +
                        hCoeffs->h22T2ff40 * ff40 + I*hCoeffs->h22T2ff31 * ff31 + hCoeffs->h22T2ff22 * ff22 + I*hCoeffs->h22T2ff13 * ff13 + 
                        hCoeffs->h22T2ff04 * ff04;
                    A3lm = hCoeffs->h22T3ff01 * ff01 + I*hCoeffs->h22T3ff10 * ff10 + 
                        hCoeffs->h22T3ff03 * ff03 + I*hCoeffs->h22T3ff30 * ff30 + 
                        hCoeffs->h22T3ff21 * ff21 + I*hCoeffs->h22T3ff12 * ff12;
                    A4lm = hCoeffs->h22T4ff00 * ff00 + hCoeffs->h22T4ff11 * ff11 + hCoeffs->h22T4ff02 * ff02 + hCoeffs->h22T4ff20 * ff20 +
                        hCoeffs->h22T4ff40 * ff40 + I*hCoeffs->h22T4ff31 * ff31 + hCoeffs->h22T4ff22 * ff22 + I*hCoeffs->h22T4ff13 * ff13 + 
                        hCoeffs->h22T4ff04 * ff04 + hCoeffs->h22T4ff60 * ff60 + I*hCoeffs->h22T4ff51 * ff51 + hCoeffs->h22T4ff42 * ff42 +
                        I*hCoeffs->h22T4ff33 * ff33 + hCoeffs->h22T4ff24 * ff24 + I*hCoeffs->h22T4ff15 * ff15 + hCoeffs->h22T4ff60 * ff60;
                    break;
                case 1:
                    rhoHigh =
                        v2*v2 * (hCoeffs->rho21v4 +
                        v * (hCoeffs->rho21v5 +
                        v * (hCoeffs->rho21v6 + hCoeffs->rho21v6l * eulerlogxabs +
                        v * (hCoeffs->rho21v7 + hCoeffs->rho21v7l * eulerlogxabs +
                        v * (hCoeffs->rho21v8 + hCoeffs->rho21v8l * eulerlogxabs +
                        (hCoeffs->rho21v10 + hCoeffs->rho21v10l * eulerlogxabs) * v2)))));
                    rhoLow =  1. + v * (hCoeffs->rho21v1 +v * (hCoeffs->rho21v2 + v * (hCoeffs->rho21v3)));
                    A0lm = g2*hCoeffs->h21T0ff01 * ff01;
                    A2lm = g2*(hCoeffs->h21T2ff21 * ff21 + hCoeffs->h21T2ff12 *I* ff12 + 
                        hCoeffs->h21T2ff03 * ff03 + hCoeffs->h21T2ff01 * ff01);
                    // (v*hCoeffs->f21v1 + v2*v*hCoeffs->f21v3)*(vPhi*vPhi2/x0_3) + 
                    Flm0 = v2*v2*(hCoeffs->f21v4 + 
                        v * (hCoeffs->f21v5 +
                        v * (hCoeffs->f21v6))) * (vPhi*vPhi2/x0_3);
                    // Flm0 = v2*v*(hCoeffs->f21v3 + 
                    //        v*(hCoeffs->f21v4 + 
                    //     v * (hCoeffs->f21v5 +
                    //     v * (hCoeffs->f21v6)))) * (vPhi*vPhi2/x0_3);

                    Flm1 = g2*hCoeffs->h21T1ff00 * ff00;
                    Flm3 = g2*(hCoeffs->h21T3ff20 * ff20 + hCoeffs->h21T3ff11 * I* ff11 + 
                        hCoeffs->h21T3ff02 * ff02 + hCoeffs->h21T3ff00 * ff00);
                    auxflm = Flm0 + Flm1 * x0 + Flm3 * x0_3;
                    // if (eta == 0.25 && m % 2)
                    // {
                        // auxflm = Flm1 * x0 + Flm0;
                        // auxflm = Flm0 + Flm1 * x0;
                        // auxflm += pow(x0, 7) * (hCoeffs->f21v7cEff00 + 
                        //     I*hCoeffs->f21v7cEff10*ff10 + I*hCoeffs->f21v7cEff11*ff11+
                        //     hCoeffs->f21v7cEff01 * ff01 + hCoeffs->f21v7cEff02*ff02);
                        auxflm += pow(x0,12) * hCoeffs->f21v7cEff00;
                        // auxflm = Flm0;
                        // auxflm = g2 * hCoeffs->f21v1 * x0 + 
                        //          v2 * v * (hCoeffs->f21v3 + 
                        //          v * (hCoeffs->f21v4 + 
                        //          v * (hCoeffs->f21v5 + 
                        //          v * (hCoeffs->f21v6 + 
                        //          v * (hCoeffs->f21v7c))))) * (vPhi*vPhi2/x0_3);
        // print_debug("f21v7c = %.16f\n", hCoeffs->f21v7c);
    // auxflm = Flm0 * cexp( x0* (Flm1/Flm0 - x0*cpow(Flm1/Flm0,2)/2. + x0_2*(cpow(Flm1/Flm0,3)/3. + Flm3/Flm0) ) );
                    // }
                    // else
                    // {
                    //     A1lm = Flm1;
                    //     A3lm = Flm3;
                    //     auxflm = Flm0;
                    // }
                    // print_debug("Flm0 = %g + i%g, Flm1 = %g + i%g, Flm3 = %g + i%g\n", Flm0, Flm1, Flm3);
                    // print_debug("v = %g, vPhi = %g, x0 = %g\n", v, vPhi, x0);
                    // print_debug("f21v4 = %g, f21v5 = %g, f21v6 = %g, f21v7 = %g\n", 
                    //     hCoeffs->f21v4, hCoeffs->f21v5, hCoeffs->f21v6, hCoeffs->f21v7c);
                    break;
                default:
                    PRINT_LOG_INFO(LOG_CRITICAL, "Currently mode (%d, %d) is unsupported in new factorized waveform.", l, m);
                    return CEV_FAILURE;

            }
            break;
        case 3:
            switch(ABS(m))
            {
                case 3:
                    rhoHigh =
                        v2*v2 * (hCoeffs->rho33v4 +
                        v * (hCoeffs->rho33v5 +
                        v * (hCoeffs->rho33v6 + hCoeffs->rho33v6l * eulerlogxabs +
                        v * (hCoeffs->rho33v7 +
                        v * (hCoeffs->rho33v8 + hCoeffs->rho33v8l * eulerlogxabs +
                        v2 * (hCoeffs->rho33v10 + hCoeffs->rho33v10l * eulerlogxabs))))));
                    rhoLow = 1. + v2 * (hCoeffs->rho33v2 + v * hCoeffs->rho33v3);
                    A0lm = I*hCoeffs->h33T0ff30 * ff30 + hCoeffs->h33T0ff21 * ff21 + 
                        I*hCoeffs->h33T0ff12 * ff12 + I*hCoeffs->h33T0ff10 * ff10 + 
                        hCoeffs->h33T0ff03 * ff03 + hCoeffs->h33T0ff01 * ff01;

                    A2lm = hCoeffs->h33T2ff50 *I*ff50 + hCoeffs->h33T2ff41*ff41 +
                        hCoeffs->h33T2ff32*I*ff32 + hCoeffs->h33T2ff30*I*ff30 +
                        hCoeffs->h33T2ff23*ff23 + hCoeffs->h33T2ff21*ff21 +
                        hCoeffs->h33T2ff14*I*ff14 + hCoeffs->h33T2ff12*I*ff12 +
                        hCoeffs->h33T2ff10*I*ff10 + hCoeffs->h33T2ff05*ff05 + 
                        hCoeffs->h33T2ff03*ff03 + hCoeffs->h33T2ff01*ff01;
                    Flm0 = (v2*v2 * (hCoeffs->f33v4 + 
                        v * (hCoeffs->f33v5 +
                        v * (hCoeffs->f33v6))) + 
                        I * vh3 * vh3 * hCoeffs->f33vh6) * (vPhi*vPhi2/x0_3);
                    
                    Flm3 = ( hCoeffs->h33T3ff20 * ff20 + 
                            hCoeffs->h33T3ff11 * I * ff11 +
                            hCoeffs->h33T3ff02 * ff02 + 
                            hCoeffs->h33T3ff00 * ff00 );
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
                    rhoHigh =
                        v2 * v * (hCoeffs->rho32v3 +
                        v * (hCoeffs->rho32v4 +
                        v * (hCoeffs->rho32v5 +
                        v * (hCoeffs->rho32v6 + hCoeffs->rho32v6l * eulerlogxabs +
                            (hCoeffs->rho32v8 + hCoeffs->rho32v8l * eulerlogxabs) * v2))));
                    rhoLow = 1. + v * (hCoeffs->rho32v + v * (hCoeffs->rho32v2));
                    A0lm = g2 * (hCoeffs->h32T0ff11 *I* ff11 + hCoeffs->h32T0ff02 * ff02);
                    A1lm = g2 * (hCoeffs->h32T1ff10 *I* ff10 + hCoeffs->h32T1ff01 * ff01);
                    A2lm = g2 * (hCoeffs->h32T2ff31 *I* ff31 + hCoeffs->h32T2ff22 * ff22 +
                        hCoeffs->h32T2ff13 *I* ff13 + hCoeffs->h32T2ff11 *I* ff11 + 
                        hCoeffs->h32T2ff04 * ff04 + hCoeffs->h32T2ff02 * ff02);
                    break;
                case 1:
                    Use_V1 = 1;
                    rhoHigh =
                        v2*v2 * (hCoeffs->rho31v4 +
                        v * (hCoeffs->rho31v5 +
                        v * (hCoeffs->rho31v6 + hCoeffs->rho31v6l * eulerlogxabs +
                        v * (hCoeffs->rho31v7 +
                            (hCoeffs->rho31v8 + hCoeffs->rho31v8l * eulerlogxabs) * v))));
                    rhoLow = 1. + v2 * (hCoeffs->rho31v2 + v * (hCoeffs->rho31v3));
                    A0lm = hCoeffs->h31T0ff30 *I* ff30 + hCoeffs->h31T0ff21*ff21 + 
                        hCoeffs->h31T0ff12*I*ff12 + hCoeffs->h31T0ff10*I*ff10 + 
                        hCoeffs->h31T0ff03*ff03 + hCoeffs->h31T0ff01 * ff01;
                    A2lm = hCoeffs->h31T2ff50 *I*ff50 + hCoeffs->h31T2ff41 * ff41+
                        hCoeffs->h31T2ff32 *I*ff32 + hCoeffs->h31T2ff30*I*ff30 +
                        hCoeffs->h31T2ff23*ff23 + hCoeffs->h31T2ff21*ff21 +
                        hCoeffs->h31T2ff14*I*ff14 + hCoeffs->h31T2ff12*I*ff12 +
                        hCoeffs->h31T2ff10*I*ff10 + hCoeffs->h31T2ff05*ff05 +
                        hCoeffs->h31T2ff03*ff03 + hCoeffs->h31T2ff01*ff01;
                    auxflm = (hCoeffs->h31T3ff20 * ff20 + 
                        hCoeffs->h31T3ff11 * I* ff11 +
                        hCoeffs->h31T3ff02 * ff02 + 
                        hCoeffs->h31T3ff00 * ff00) * x0_3;
                    break;
                default:
                    PRINT_LOG_INFO(LOG_CRITICAL, "Currently mode (%d, %d) is unsupported in new factorized waveform.", l, m);
                    return CEV_FAILURE;
            }
            break;
        case 4:
            switch(ABS(m))
            {
                case 4:
                    rhoHigh = v * v2 * (hCoeffs->rho44v3 + 
                        v * (hCoeffs->rho44v4 +
                        v * (hCoeffs->rho44v5 + 
                        v * ((hCoeffs->rho44v6 + hCoeffs->rho44v6l * eulerlogxabs) + 
                        v2 * ((hCoeffs->rho44v8 + hCoeffs->rho44v8l * eulerlogxabs) + 
                        v2 * ((hCoeffs->rho44v10 + hCoeffs->rho44v10l * eulerlogxabs)))))));
                    rhoLow = 1. + v2 * hCoeffs->rho44v2;
                    A0lm = hCoeffs->h44T0ff40 * ff40 + hCoeffs->h44T0ff31*I*ff31 +
                        hCoeffs->h44T0ff22*ff22 + hCoeffs->h44T0ff20*ff20 +
                        hCoeffs->h44T0ff13*I*ff13 + hCoeffs->h44T0ff11*I*ff11 +
                        hCoeffs->h44T0ff04*ff04 + hCoeffs->h44T0ff02*ff02 + 
                        hCoeffs->h44T0ff00 * ff00;
                    A2lm = hCoeffs->h44T2ff60*ff60 + hCoeffs->h44T2ff51*I*ff51 +
                        hCoeffs->h44T2ff42*ff42 + hCoeffs->h44T2ff40*ff40 + 
                        hCoeffs->h44T2ff31*I*ff31 + hCoeffs->h44T2ff24*ff24+
                        hCoeffs->h44T2ff22*ff22+ hCoeffs->h44T2ff20*ff20 + 
                        hCoeffs->h44T2ff15*I*ff15 + hCoeffs->h44T2ff13*I*ff13 + 
                        hCoeffs->h44T2ff11*I*ff11 + hCoeffs->h44T2ff06*ff06 + 
                        hCoeffs->h44T2ff04*ff04 + hCoeffs->h44T2ff02*ff02 + 
                        hCoeffs->h44T2ff00*ff00;

                    break;
                case 3:
                    rhoHigh =
                        v2 * (hCoeffs->rho43v2 +
                        v2 * (hCoeffs->rho43v4 +
                        v * (hCoeffs->rho43v5 + 
                        (hCoeffs->rho43v6 + hCoeffs->rho43v6l * eulerlogxabs) * v)));
                    rhoLow = 1.;
                    A0lm = g2*(hCoeffs->h43T0ff21*ff21 + hCoeffs->h43T0ff12*I*ff12 + 
                        hCoeffs->h43T0ff03*ff03 + hCoeffs->h43T0ff01*ff01);
                    
                    auxflm = x0 * g2 * (hCoeffs->h43T1ff20*ff20 + hCoeffs->h43T1ff11 *I*ff11 +
                        hCoeffs->h43T1ff02*ff02 + hCoeffs->h43T1ff00*ff00);
                    break;
                case 2:
                    Use_V1 = 1;
                    rhoHigh = 
                        v2*v * (hCoeffs->rho42v3 +
                        v * (hCoeffs->rho42v4 +
                        v * (hCoeffs->rho42v5 +
                            (hCoeffs->rho42v6 + hCoeffs->rho42v6l * eulerlogxabs) * v)));
                    rhoLow = 1. + v2 * hCoeffs->rho42v2;
                    A0lm = hCoeffs->h42T0ff40 * ff40 + hCoeffs->h42T0ff31*I*ff31 + 
                        hCoeffs->h42T0ff20*ff20 + hCoeffs->h42T0ff13*I*ff13 + 
                        hCoeffs->h42T0ff11*I*ff11 + hCoeffs->h42T0ff04*ff04 + 
                        hCoeffs->h42T0ff02*ff02 + hCoeffs->h42T0ff00 * ff00;
                    A2lm = hCoeffs->h42T2ff60*ff60 + hCoeffs->h42T2ff51*I*ff51 +
                        hCoeffs->h42T2ff42*ff42 + hCoeffs->h42T2ff40*ff40 + 
                        hCoeffs->h42T2ff33*I*ff33 + hCoeffs->h42T2ff31*I*ff31 +
                        hCoeffs->h42T2ff24*ff24 + hCoeffs->h42T2ff22*ff22 +
                        hCoeffs->h42T2ff20*ff20 + hCoeffs->h42T2ff15*I*ff15 + 
                        hCoeffs->h42T2ff13*I*ff13 + hCoeffs->h42T2ff11*I*ff11 +
                        hCoeffs->h42T2ff06*ff06 + hCoeffs->h42T2ff04 * ff04 +
                        hCoeffs->h42T2ff02*ff02 + hCoeffs->h42T2ff00 * ff00;
                    break;
                case 1:
                    Use_V1 = 1;
                    rhoHigh =
                        v2 * (hCoeffs->rho41v2 +
                        v2 * (hCoeffs->rho41v4 +
                        v * (hCoeffs->rho41v5 +
                            (hCoeffs->rho41v6 + hCoeffs->rho41v6l * eulerlogxabs) * v)));
                    rhoLow = 1.;
                    A0lm = g2 * (hCoeffs->h41T0ff21*ff21 + hCoeffs->h41T0ff12*I*ff12 +
                        hCoeffs->h41T0ff03*ff03 + hCoeffs->h41T0ff01*ff01);
                    
                    auxflm = x0 * g2 * (hCoeffs->h41T1ff20*ff20 + hCoeffs->h41T1ff11*I*ff11 + 
                        hCoeffs->h41T1ff02*ff02 + hCoeffs->h41T1ff00*ff00);
                    break;
                default:
                    PRINT_LOG_INFO(LOG_CRITICAL, "Currently mode (%d, %d) is unsupported in new factorized waveform.", l, m);
                    return CEV_FAILURE;
            }
            break;

        case 5:
            Use_V1 = 1;
            switch (ABS(m))
            {
                case 5:
                    //RC: This terms are in Eq.A9 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
                    rhoLow =
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
                    rhoLow = 1. + v2 * (hCoeffs->rho54v2 + v * (hCoeffs->rho54v3
                                            + hCoeffs->rho54v4 * v));
                    break;
                case 3:
                    rhoLow = 1. + v2 * (hCoeffs->rho53v2
                                + v * (hCoeffs->rho53v3 +
                                    v * (hCoeffs->rho53v4 +
                                    hCoeffs->rho53v5 * v)));
                    break;
                case 2:
                    rhoLow = 1. + v2 * (hCoeffs->rho52v2 + v * (hCoeffs->rho52v3
                                            + hCoeffs->rho52v4 * v));
                    break;
                case 1:
                    rhoLow = 1. + v2 * (hCoeffs->rho51v2
                                + v * (hCoeffs->rho51v3 +
                                    v * (hCoeffs->rho51v4 +
                                    hCoeffs->rho51v5 * v)));
                    break;
                default:
                    PRINT_LOG_INFO(LOG_CRITICAL, "Currently mode (%d, %d) is unsupported in new factorized waveform.", l, m);
                    return CEV_FAILURE;
            }
            break;
        case 6:
            Use_V1 = 1;
            switch (ABS(m))
            {
                case 6:
                    rhoLow = 1. + v2 * (hCoeffs->rho66v2 + v * (hCoeffs->rho66v3
                                            + hCoeffs->rho66v4 * v));
                    break;
                case 5:
                    rhoLow = 1. + v2 * (hCoeffs->rho65v2 + hCoeffs->rho65v3 * v);
                    break;
                case 4:
                    rhoLow = 1. + v2 * (hCoeffs->rho64v2 + v * (hCoeffs->rho64v3
                                        + hCoeffs->rho64v4 * v));
                    break;
                case 3:
                    rhoLow = 1. + v2 * (hCoeffs->rho63v2 + hCoeffs->rho63v3 * v);
                    break;
                case 2:
                    rhoLow = 1. + v2 * (hCoeffs->rho62v2 + v * (hCoeffs->rho62v3
                                        + hCoeffs->rho62v4 * v));
                    break;
                case 1:
                    rhoLow = 1. + v2 * (hCoeffs->rho61v2 + hCoeffs->rho61v3 * v);
                    break;
                default:
                    PRINT_LOG_INFO(LOG_CRITICAL, "Currently mode (%d, %d) is unsupported in new factorized waveform.", l, m);
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
                    PRINT_LOG_INFO(LOG_CRITICAL, "Currently mode (%d, %d) is unsupported in new factorized waveform.", l, m);
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
                    PRINT_LOG_INFO(LOG_CRITICAL, "Currently mode (%d, %d) is unsupported in new factorized waveform.", l, m);
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
        REAL8 LL = (REAL8) l;
        T0 = 1.;
        T1 = -A1lm/A0lm;
        T2 = cpow(A1lm/A0lm, 2.) - A2lm/A0lm;
        T3 = -cpow(A1lm/A0lm, 3) + 2.*A1lm*A2lm/A0lm/A0lm - A3lm/A0lm;
        T4 = cpow(A1lm/A0lm, 4.) - 3.*A1lm*A1lm*A2lm/cpow(A0lm, 3.) + cpow(A2lm/A0lm, 2) + 2.*A1lm*A3lm/A0lm - A4lm/A0lm;
        rho2PN = A0lm / (T0 + T1*x0 + T2 * x0_2 + T3 * x0_3 + T4 * x0_4);
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
        rholmPwrl = rho2PN + (rholmPwrl - rholmPwrlSub) * pow(vPhi/x0, l + (l + m) % 2);
        // rholmPwrl = rholmPwrl * pow(vPhi/x0, l + (l + m) % 2);
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
INT
XLALSimIMRSpinEOBGetSpinFactorizedWaveformV3 (COMPLEX16 * ret,
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
    REAL8 r, pp, Omega, v2, vh, vh3, k, hathatk, eulerlogxabs;	//pr
    REAL8 Slm, deltalm, rholm;
    COMPLEX16 auxflm = 0.0;
    COMPLEX16 Tlm, rholmPwrl;
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
    //pr    = values->data[2];
    pp = values->data[3];

    v2 = v * v;
    Omega = v2 * v;
    vh3 = Hreal * Omega;
    vh = cbrt (vh3);
    eulerlogxabs = CST_GAMMA + log (2.0 * (REAL8) m * v);

    /* Calculate the non-Keplerian velocity */
    if (params->alignedSpins)
    {
        // YP: !!!!! SEOBNRv3devel temporary change !!!!!
        vPhi =
            XLALSimIMRSpinAlignedEOBNonKeplerCoeff (values->data, params);

        if (IS_REAL8_FAIL_NAN (vPhi))
        {
            return CEV_FAILURE;
        }

        vPhi = r * cbrt (vPhi);
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
        status = XLALSimIMRSpinEOBCalculateNewtonianMultipole (&hNewton, vPhi2, r,
                                    values->data[1],
                                    (UINT) l, m,
                                    params);
    else
        status = XLALSimIMRSpinEOBFluxCalculateNewtonianMultipole (&hNewton,
							     vPhi2, r,
							     values->data[1],
							     (UINT) l, m,
							     params);
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
    //printf( "Hreal = %e, Slm = %e, eta = %e\n", Hreal, Slm, eta );

    /* Calculate the Tail term, 3rd term in Eq. 17, given by Eq. A6 */
    k = m * Omega;
    hathatk = Hreal * k;
    if (ret_type <= 0)
    {
        status = gsl_sf_lngamma_complex_e (l + 1.0, -2.0 * hathatk, &lnr1, &arg1);
        if (status != GSL_SUCCESS)
        {
            PRINT_LOG_INFO (LOG_CRITICAL, "Error in GSL function: %s", gsl_strerror (status));
            return CEV_FAILURE;
        }
        status = gsl_sf_fact_e (l, &z2);
        if (status != GSL_SUCCESS)
        {
            PRINT_LOG_INFO (LOG_CRITICAL, "Error in GSL function: %s", gsl_strerror (status));
            return CEV_FAILURE;
        }
        Tlm = cexp ((lnr1.val + CST_PI * hathatk) + I * (arg1.val + 2.0 * hathatk * log (4.0 * k / sqrt (CST_E))));
        Tlm /= z2.val;
        switch (l)
        {
            case 2:
                switch (abs (m))
                {
                    case 2:
                        deltalm = 
                            vh3 * (hCoeffs->delta22vh3 + 
                            vh3 * (hCoeffs->delta22vh6 +
                            vh * vh *(hCoeffs->delta22vh9 * vh))) +
                            hCoeffs->delta22v5 * v * v2 * v2 +
                            hCoeffs->delta22v6 * v2 * v2 * v2 +
                            hCoeffs->delta22v8 * v2 * v2 * v2 * v2;
                        break;
                    case 1:
                        deltalm = 
                            vh3 * (hCoeffs->delta21vh3 + 
                            vh3 * (hCoeffs->delta21vh6 +
                            vh * (hCoeffs->delta21vh7 +
                            vh * vh*(hCoeffs->delta21vh9)))) +
                            hCoeffs->delta21v5 * v * v2 * v2 +
                            hCoeffs->delta21v7 * v2 * v2 * v2 * v;
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
                        deltalm =
                            vh3 * (hCoeffs->delta33vh3 +
                            vh3 * (hCoeffs->delta33vh6 + hCoeffs->delta33vh9 * vh3)) +
                            hCoeffs->delta33v5 * v * v2 * v2;
                        break;
                    case 2:
                        deltalm =
                            vh3 * (hCoeffs->delta32vh3 +
                            vh * (hCoeffs->delta32vh4 +
                            vh * vh * (hCoeffs->delta32vh6 +
                            hCoeffs->delta32vh9 * vh3)));
                        break;
                    case 1:
                        deltalm = 
                            vh3 * (hCoeffs->delta31vh3 + 
                            vh3 * (hCoeffs->delta31vh6 +
                            vh * (hCoeffs->delta31vh7 + 
                            vh * vh * hCoeffs->delta31vh9))) +
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
                        if(use_hm)
                        {
                            //RC: This terms are in Eq.A15 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028
                            deltalm = 
                                vh3 * (hCoeffs->delta44vh3 + 
                                vh3 * (hCoeffs->delta44vh6 + 
                                vh3 * hCoeffs->delta44vh9))
                                + hCoeffs->delta44v5 * v2 * v2 * v;
                        }
                        else
                        {
                            deltalm = 
                                vh3 * (hCoeffs->delta44vh3 + 
                                vh3 * hCoeffs->delta44vh6) + 
                                hCoeffs->delta44v5 * v2 * v2 * v;
                        }
                        break;
                case 3:
                    deltalm = 
                        vh3 * (hCoeffs->delta43vh3 + 
                        vh * (hCoeffs->delta43vh4 +
                        vh * vh * hCoeffs->delta43vh6));
                    break;
                case 2:
                    deltalm = vh3 * (hCoeffs->delta42vh3 + hCoeffs->delta42vh6 * vh3);
                    break;
                case 1:
                    deltalm = 
                        vh3 * (hCoeffs->delta41vh3 + 
                        vh * (hCoeffs->delta41vh4 +
                        vh * vh * hCoeffs->delta41vh6));
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
                        if(use_hm)
                        {
                            //RC: This terms are in Eq.A16 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028
                            deltalm =
                                vh3 *(hCoeffs->delta55vh3 +
                                vh3*(hCoeffs->delta55vh6 +
                                vh3 *(hCoeffs->delta55vh9))) + 
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
        status = gsl_sf_fact_e (l, &z2);
        if (status != GSL_SUCCESS)
        {
            PRINT_LOG_INFO (LOG_CRITICAL, "Error in GSL function: %s", gsl_strerror (status));
            return CEV_FAILURE;
        }
        /* Calculating the prefactor of Tlm, outside the multiple product */
        Tlmprefac = sqrt (hathatk4pi / (1. - exp (-hathatk4pi))) / z2.val;

        /* Calculating the multiple product factor */
        for (Tlmprodfac = 1., i = 1; i <= l; i++)
        {
            Tlmprodfac *= (hathatksq4 + (REAL8) i * i);
        }

        Tlm = Tlmprefac * sqrt (Tlmprodfac);
        deltalm = 0.0;
    }
    if (CalculateFacWaveformAmpResV3(hCoeffs, eta, l, m, v, vPhi, vh3, eulerlogxabs, 1/sqrt(r), dr, ncrv, &rholmPwrl) != CEV_SUCCESS)
        return CEV_FAILURE;
    if (ret_type >1)
    {
        *ret = rholmPwrl;
        return CEV_SUCCESS;
    }
    /* Put all factors in Eq. 17 together */
    *ret = Tlm * cexp (I * deltalm) * Slm * rholmPwrl;
    *ret *= hNewton;
    /*if (r > 8.5)
    {
        printf("YP::FullWave: Reh = %.16e, Imh = %.16e, hAmp = %.16e, hPhi = %.16e\n",creal(*hlm),cimag(*hlm),cabs(*hlm),carg(*hlm));
    } */
    return CEV_SUCCESS;
}

REAL8
XLALInspiralSpinFactorizedFluxV3 (REAL8Vector * values,	/**< dynamical variables */
				EOBNonQCCoeffs * nqcCoeffs,
							/**< pre-computed NQC coefficients */
				const REAL8 omega,	/**< orbital frequency */
                const REAL8 dr, /**< radial velocity */
                const REAL8 ncrv, /**< nVec cross vVec */
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
            if (XLALSimIMRSpinEOBGetSpinFactorizedWaveformV3
                (&hLM, values, v, dr, ncrv, H, l, m, ak, 1) != CEV_SUCCESS)
            {
                return REAL8_FAIL_NAN;
            }

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



/*--------------------------------------------------------------*/
/*                                                              */
/*                                                              */
/*                                                              */
/*                            PREC                              */
/*                                                              */
/*                                                              */
/*                                                              */
/*--------------------------------------------------------------*/
INT
XLALSimIMRSpinEOBFluxGetPrecSpinFactorizedWaveformV3(
                        COMPLEX16 *  hlmTab,	/**< OUTPUT, hlm waveforms */
                        REAL8Vector *  values,	/**< dyanmical variables: \f$(r,\phi,p_r,p_\phi)\f$ */
                        REAL8Vector *  cartvalues,	/**< dyanmical variables */
                        const REAL8 v,	/**< velocity */
                        const REAL8 dr,
                        const REAL8 ncrv,
                        const REAL8 Hreal,	/**< real Hamiltonian */
                        const INT4 lMax,	/**< maximum l mode to compute, compute 0 < m <= lMax */
                        SpinEOBParams *  params	/**< Spin EOB parameters */
)
{
    // int		debugPK = 0;
	const	REAL8 vPhiKepler = params->alignedSpins ?
					XLALSimIMRSpinAlignedEOBNonKeplerCoeff(values->data, params) :
					XLALSimIMRSpinPrecEOBNonKeplerCoeff(cartvalues->data, params);
	if (IS_REAL8_FAIL_NAN(vPhiKepler)) {
		return CEV_FAILURE;
	}
    INT l, m;
	for ( l = 2; l <= lMax; l++)
	{
		for ( m = 1; m <= l; m++)
		{
			COMPLEX16 *hlm = &hlmTab[l*(lMax+1)+m];
			/* Status of function calls */
			INT		status;
			INT		i;

			REAL8		eta;
			REAL8 r , pp, Omega, v2, vh3,/* vh, vh3, */ k, hathatk, eulerlogxabs;
			//pr
			REAL8  rcrossp_x, rcrossp_y, rcrossp_z;
			REAL8		Slm     , rholm;
            COMPLEX16 rholmPwrl;
			REAL8		auxflm = 0.0;
			REAL8		hathatksq4, hathatk4pi, Tlmprefac, Tlmprodfac;
			REAL8		Tlm;
			COMPLEX16	hNewton;
			gsl_sf_result	z2;

			/* Non-Keplerian velocity */
			REAL8		vPhi    , vPhi2;

			/* Pre-computed coefficients */
			FacWaveformCoeffs *hCoeffs = params->hCoeffs;

			eta = params->eta;

			/*
			 * else if ( eta == 0.25 && m % 2 ) { // If m is odd and dM = 0, hLM
			 * will be zero memset( hlm, 0, sizeof( COMPLEX16 ) ); return
			 * XLAL_SUCCESS; }
			 */

			//r = sqrt(values->data[0] * values->data[0] + values->data[1] * values->data[1] + values->data[2] * values->data[2]);
			//pr = values->data[2];
			r = values->data[0];
			pp = values->data[3];

			rcrossp_x = cartvalues->data[1] * cartvalues->data[5] - cartvalues->data[2] * cartvalues->data[4];
			rcrossp_y = cartvalues->data[2] * cartvalues->data[3] - cartvalues->data[0] * cartvalues->data[5];
			rcrossp_z = cartvalues->data[0] * cartvalues->data[4] - cartvalues->data[1] * cartvalues->data[3];

			//pp = sqrt(rcrossp_x * rcrossp_x + rcrossp_y * rcrossp_y + rcrossp_z * rcrossp_z);

			v2 = v * v;
			Omega = v2 * v;
			vh3 = Hreal * Omega;
			//vh = cbrt(vh3);
			eulerlogxabs = CST_GAMMA + log(2.0 * (REAL8) m * v);

			/* Calculate the non-Keplerian velocity */
			vPhi = vPhiKepler;

			vPhi = r * cbrt(vPhi);

			// if (debugPK)
			// 	XLAL_PRINT_INFO("In XLALSimIMRSpinEOBFluxCalculateNewtonianMultipole, getting rW = %.12e\n", vPhi);
			vPhi *= Omega;
			vPhi2 = vPhi * vPhi;

			/*
			 * Calculate the newtonian multipole, 1st term in Eq. 17, given by
			 * Eq. A1
			 */
			//debugPK
			// if (debugPK) {
			// 	XLAL_PRINT_INFO("\nValues inside XLALSimIMRSpinEOBFluxGetPrecSpinFactorizedWaveform:\n");
			// 	for (i = 0; i < 14; i++)
			// 		XLAL_PRINT_INFO("values[%d] = %.12e\n", i, values->data[i]);

			// 	XLAL_PRINT_INFO("Calculating hNewton, with v = %.12e, vPhi = %.12e, r = %.12e, Phi = %.12e, l = %d, m = %d\n",
			// 	       v, vPhi, r, values->data[1], (UINT4) l, (UINT4) m);
			// }
			status = XLALSimIMRSpinEOBFluxCalculateNewtonianMultipole(&hNewton,
				vPhi2, r, values->data[1], (UINT) l, m, params);
			if (status == CEV_FAILURE) {
				return CEV_FAILURE;
			}
			/* Calculate the source term, 2nd term in Eq. 17, given by Eq. A5 */
			if (((l + m) % 2) == 0) {
				Slm = (Hreal * Hreal - 1.) / (2. * eta) + 1.;
			} else {
				Slm = v * pp;
				//Slm = v * sqrt(rcrossp_x * rcrossp_x + rcrossp_y * rcrossp_y + rcrossp_z * rcrossp_z);
			}
			// if (debugPK)
			// 	XLAL_PRINT_INFO("In XLALSimIMRSpinEOBFluxGetSpinFactorizedWaveform: Hreal = %e, Slm = %e, eta = %e\n", Hreal, Slm, eta);

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
				Tlmprodfac *= (hathatksq4 + (REAL8) i * i);
			}

			Tlm = Tlmprefac * sqrt(Tlmprodfac);

            if (CalculateFacWaveformAmpResV3(hCoeffs, eta, l, m, v, vPhi, vh3, eulerlogxabs, 1/sqrt(r), dr, ncrv, &rholmPwrl) != CEV_SUCCESS)
                return CEV_FAILURE;
			// if (r > 0.0 && debugPK) {
			// 	XLAL_PRINT_INFO("YP::dynamics variables in waveform: %i, %i, r = %.12e, v = %.12e\n", l, m, r, v);
			// 	XLAL_PRINT_INFO("rholm^l = %.16e, Tlm = %.16e + i %.16e, \nSlm = %.16e, hNewton = %.16e + i %.16e, delta = %.16e\n", rholmPwrl, Tlm, 0.0, Slm, creal(hNewton), cimag(hNewton), 0.0);
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
 * described in ￼the section "Factorized waveforms" of https://dcc.ligo.org/T1400476
 */
INT
XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveformV3(
            COMPLEX16 * hlm,	/**< OUTPUT, hlm waveforms */
            REAL8Vector * values,	/**< dyanmical variables: \f$(r,\phi,p_r,p_\phi)\f$ */
            REAL8Vector * cartvalues,	/**< dyanmical variables */
            const REAL8 v,	/**< velocity */
            const REAL8 dr, /**< radial velocity */
            const REAL8 ncrv, /**< angular velocity */
            const REAL8 Hreal,	/**< real Hamiltonian */
            const INT l,	/**< l mode index */
            const INT m,	/**< m mode index */
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
    switch (l) {
    case 2:
        switch (abs(m)) {
        case 2:
            deltalm = vh3 * (hCoeffs->delta22vh3 + vh3 * (hCoeffs->delta22vh6
                    + vh * vh * (hCoeffs->delta22vh9 * vh)))
                + Omega*(hCoeffs->delta22v5 * v2 + Omega*(hCoeffs->delta22v6  + hCoeffs->delta22v8 *v2));
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
            break;
        case 2:
            deltalm = vh3 * (hCoeffs->delta32vh3 + vh * (hCoeffs->delta32vh4 + vh * vh * (hCoeffs->delta32vh6
                            + hCoeffs->delta32vh9 * vh3)));
            break;
        case 1:
            deltalm = vh3 * (hCoeffs->delta31vh3 + vh3 * (hCoeffs->delta31vh6
                                        + vh * (hCoeffs->delta31vh7 + hCoeffs->delta31vh9 * vh * vh)))
                + hCoeffs->delta31v5 * v * v2 * v2;
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
            break;
        case 3:
            deltalm = vh3 * (hCoeffs->delta43vh3 + vh * (hCoeffs->delta43vh4
                        + hCoeffs->delta43vh6 * vh * vh));
            break;
        case 2:
            deltalm = vh3 * (hCoeffs->delta42vh3 + hCoeffs->delta42vh6 * vh3);
            break;
        case 1:
            deltalm = vh3 * (hCoeffs->delta41vh3 + vh * (hCoeffs->delta41vh4
                        + hCoeffs->delta41vh6 * vh * vh));
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
        switch (m) {
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
        switch (m) {
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
        switch (m) {
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

    if (CalculateFacWaveformAmpResV3(hCoeffs, eta, l, m, v, vPhi, vh3, eulerlogxabs, 1/sqrt(r), dr, ncrv, &rholmPwrl) != CEV_SUCCESS)
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Failure in CalculateFacWaveformAmpResV3");
        return CEV_FAILURE;
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


INT XLALSimIMRSpinEOBGetAmplitudeResidualPrecV3(COMPLEX16 *rholmpwrl, 
                REAL8Vector *  values,	/**< dyanmical variables: \f$(r,\phi,p_r,p_\phi)\f$ */
                REAL8Vector *  cartvalues,	/**< dyanmical variables */
                const REAL8 v, 
                const REAL8 dr,
                const REAL8 ncrv,
                const REAL8 Hreal, 
                const INT modeL, const INT modeM, 
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
    eulerlogxabs = CST_GAMMA + log (2.0 * (REAL8) modeM * v);
	REAL8 vPhi;

    if (abs (modeM) > (INT4) modeL)
    {
        return CEV_FAILURE;
    }
    if (modeM == 0)
    {
        return CEV_FAILURE;
    }
    eta = params->eta;

    /* Check our eta was sensible */
    if (eta > 0.25 && eta < 0.25 +1e-4) {
        eta = 0.25;
    }
    if (eta > 0.25)
    {
        PRINT_LOG_INFO
        (LOG_CRITICAL, "Eta seems to be > 0.25 - this isn't allowed!");
        return CEV_FAILURE;
    }
    r = values->data[0];
    Omega = v2 * v;
    vh3 = Hreal * Omega;
    //vh = cbrt(vh3);

    vPhi = params->alignedSpins ?
            XLALSimIMRSpinAlignedEOBNonKeplerCoeff(values->data, params) :
            XLALSimIMRSpinPrecEOBNonKeplerCoeff(cartvalues->data, params);
    /* Calculate the non-Keplerian velocity */

    vPhi = r * cbrt(vPhi);

    // if (debugPK)
    // 	XLAL_PRINT_INFO("In XLALSimIMRSpinEOBFluxCalculateNewtonianMultipole, getting rW = %.12e\n", vPhi);
    vPhi *= Omega;

    if (CalculateFacWaveformAmpResV3(hCoeffs, eta, modeL, modeM, v, vPhi, vh3, eulerlogxabs, 1/sqrt(r), dr, ncrv, &rholmPwrl) != CEV_SUCCESS)
        return CEV_FAILURE;
    *rholmpwrl = rholmPwrl;
    return CEV_SUCCESS;
}


REAL8
XLALInspiralPrecSpinFactorizedFluxV3(
				   REAL8Vector * polvalues,	/**< \f$(r,\phi,p_r,p_\phi)\f$ */
				   REAL8Vector * values,	/**< dynamical variables */
				   EOBNonQCCoeffs * nqcCoeffs,	/**< pre-computed NQC coefficients */
				   const REAL8 omega,	/**< orbital frequency */
				   SpinEOBParams * ak,	/**< physical parameters */
				   const REAL8 H,	/**< real Hamiltonian */
				   const INT4 lMax,	/**< upper limit of the summation over l */
				   const UINT SpinAlignedEOBversion	/**< 1 for SEOBNRv1, 2 for SEOBNRv2 */
)
{
    // int	debugPK = 0;
    int i = 0;
    double radius = sqrt(values->data[0]*values->data[0] + values->data[1] *values->data[1]  + values->data[2] *values->data[2]  );
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