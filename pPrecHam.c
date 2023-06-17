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
#include "pPrecHam.h"
#include "pPrecRRForce.h"
#include "pEnergyFlux.h"
#include "pFactorizedWaveform.h"
#include "newFactorizedWaveformPrec.h"
#include <gsl/gsl_deriv.h>

/*--------------------------------------------------------------*/
/*                                                              */
/*                                                              */
/*                                                              */
/*                            PREC                              */
/*                                                              */
/*                                                              */
/*                                                              */
/*--------------------------------------------------------------*/
#define KRONECKER_DELTA(i,j) ((REAL8)((i)==(j)))

/** --------------------------------------------------- **/
/**                                                     **/
/**                                                     **/
/**                    deriv function                   **/
/**                                                     **/
/**                                                     **/
/** --------------------------------------------------- **/
#define KRONECKER_DELTA(i,j) ((REAL8)((i)==(j)))

int PrecHcapNumericalDerivative(double t,
                                const REAL8 values[],
                                REAL8 dvalues[],
                                void *funcParams)
{
	// int		debugPK = 0;
    /** lMax: l index up to which h_{lm} modes are included in the computation of the GW enegy flux: see Eq. in 13 in PRD 86,  024011 (2012) */
    static const INT lMax = 8;

	HcapDerivParams	params;

	/* Since we take numerical derivatives wrt dynamical variables */
	/* but we want them wrt time, we use this temporary vector in  */
	/* the conversion */
	REAL8	tmpDValues[14];

	REAL8	H;
	//Hamiltonian
    REAL8 flux;

	gsl_function	F;
	INT4		gslStatus;
	UINT		SpinAlignedEOBversion = 4;
	/* This is needed because SpinAlignedEOBversion is set to 2 for v3 */
	/* while XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients requires 3 ... */
	UINT		SpinAlignedEOBversionForWaveformCoefficients;

	UINT		i       , j, k, l, jj;

	REAL8Vector	rVec, pVec;
	REAL8		rData[3],  pData[3];

	/* We need r, phi, pr, pPhi to calculate the flux */
	REAL8		r;
	REAL8Vector	polarDynamics, cartDynamics;
	REAL8		polData  [4];

	REAL8		mass1   , mass2, eta;
	REAL8 	rrTerm2, pDotS1, pDotS2;
	REAL8Vector	s1 , s2, s1norm, s2norm, sKerr, sStar;
	REAL8		s1Data   [3], s2Data[3], s1DataNorm[3], s2DataNorm[3];
	REAL8		sKerrData[3], sStarData[3];
	REAL8   chiS, chiA, a, tplspin;
	REAL8 	s1dotLN, s2dotLN;


	/* Orbital angular momentum */
	REAL8		Lx      , Ly, Lz, magL;
	REAL8		Lhatx   , Lhaty, Lhatz;
	REAL8		dLx     , dLy, dLz;
	REAL8		dLhatx  , dLhaty, dMagL;

	REAL8		alphadotcosi;

	REAL8		rCrossV_x, rCrossV_y, rCrossV_z, omega;

	/* The error in a derivative as measured by GSL */
	REAL8		absErr;

	REAL8		tmpP[3], rMag, rMag2, prT;
	REAL8		u, u2, u3, u4, u5, w2, a2;
	REAL8		D, m1PlusetaKK, bulk, logTerms, deltaU, deltaT, deltaR;
    REAL8  		eobD_r, deltaU_u, deltaU_r;
	REAL8		dcsi, csi;

	REAL8		tmpValues[12];
	REAL8		Tmatrix  [3][3], invTmatrix[3][3], dTijdXk[3][3][3];
	REAL8		tmpPdotT1[3], tmpPdotT2[3], tmpPdotT3[3];
	//3 terms of Eq.A5

	/* NQC coefficients container */
    EOBNonQCCoeffs * nqcCoeffs = NULL;
	SEOBPrecVariables vars;
	prec_SetSEOBPrecVariables(&vars, values, values+3, values+6, values+9);
	/* Set up pointers for GSL */
	params.values = values;
	params.params = (SpinEOBParams *) funcParams;
	nqcCoeffs = params.params->nqcCoeffs;

	F.function = &GSLSpinPrecHamiltonianWrapper;
	F.params = &params;

	mass1 = params.params->m1;
	mass2 = params.params->m2;
	// SO: Rescale the masses so that the total mass is 1
	REAL8 m_total = mass1 + mass2;
	mass1 /=m_total;
	mass2 /=m_total;
	eta = params.params->eta;
	SpinAlignedEOBversion = 4;
	SpinEOBHCoeffs *coeffs = (SpinEOBHCoeffs *) params.params->seobCoeffs;
	REAL8 STEP_SIZE; // The step size passed to GSL to compute derivatives
    STEP_SIZE = 2.0e-3; //Allow a different step size for v4P 
	/*
	 * For precessing binaries, the effective spin of the Kerr background
	 * evolves with time. The coefficients used to compute the
	 * Hamiltonian depend on the Kerr spin, and hence need to be updated
	 * for the current spin values
	 */

	/*
	 * Set the position/momenta vectors to point to the appropriate
	 * things
	 */
    /* Here pvec is the reduced tortoise p^* vector of Pan et al. PRD 81, 084041 (2010) */
	rVec.length = pVec.length = 3;
	rVec.data = rData;
	pVec.data = pData;
	memcpy(rData, values, sizeof(rData));
	memcpy(pData, values + 3, sizeof(pData));

	/*
	 * We need to re-calculate the parameters at each step as precessing
	 * spins will not be constant
	 */
	/* TODO: Modify so that only spin terms get re-calculated */

	/*
	 * We cannot point to the values vector directly as it leads to a
	 * warning
	 */
	s1.length = s2.length = s1norm.length = s2norm.length = 3;
	s1.data = s1Data;
	s2.data = s2Data;
	s1norm.data = s1DataNorm;
	s2norm.data = s2DataNorm;

	memcpy(s1Data, values + 6, 3 * sizeof(REAL8));
	memcpy(s2Data, values + 9, 3 * sizeof(REAL8));
	memcpy(s1DataNorm, values + 6, 3 * sizeof(REAL8));
	memcpy(s2DataNorm, values + 9, 3 * sizeof(REAL8));
	for (i = 0; i < 3; i++) 
    {
		s1.data[i] *= (mass1 + mass2) * (mass1 + mass2);
		s2.data[i] *= (mass1 + mass2) * (mass1 + mass2);
	}
	sKerr.length = 3;
	sKerr.data = sKerrData;
	EOBCalculateSigmaKerr(&sKerr, &s1norm, &s2norm);

	sStar.length = 3;
	sStar.data = sStarData;
	EOBCalculateSigmaStar(&sStar, mass1, mass2, &s1norm, &s2norm);

	a = sqrt(sKerr.data[0] * sKerr.data[0] + sKerr.data[1] * sKerr.data[1]
		 + sKerr.data[2] * sKerr.data[2]);

	/* Convert momenta to p Eq. A3 of PRD 81, 084041 (2010) */
	// rMag = sqrt(rVec.data[0] * rVec.data[0] + rVec.data[1] * rVec.data[1] + rVec.data[2] * rVec.data[2]);
	rMag = vars.r;
    /* This is p^*.r/|r| */
	// prT = pVec.data[0] * (rVec.data[0] / rMag) + pVec.data[1] * (rVec.data[1] / rMag)
	// 	+ pVec.data[2] * (rVec.data[2] / rMag);
	prT = vars.prT;

	rMag2 = rMag * rMag;
	u = 1. / rMag;
	u2 = u * u;
	u3 = u2 * u;
	u4 = u2 * u2;
	u5 = u4 * u;
	a2 = a * a;
	w2 = rMag2 + a2;
	/* Eq. 5.83 of BB1, inverse */
	D = 1. + log(1. + 6. * eta * u2 + 2. * (26. - 3. * eta) * eta * u3);
     /* d(Eq. 5.83 of BB1)/dr */
	eobD_r = (u2 / (D * D)) * (12. * eta * u + 6. * (26. - 3. * eta) * eta * u2) / (1.
			+ 6. * eta * u2 + 2. * (26. - 3. * eta) * eta * u3);
	m1PlusetaKK = -1. + eta * coeffs->KK;
	/* Eq. 5.75 of BB1 */
    /* a as returned by XLALSimIMRSpinEOBCalculateSigmaKerr is S/M^2 so that a is unitless, i.e. 1/M^2 is absorbed in a2. Also, u = M/|r| is unitless */
	bulk = 1. / (m1PlusetaKK * m1PlusetaKK) + (2. * u) / m1PlusetaKK + a2 * u2;
	/* Eq. 5.73 of BB1 */
    /* The 4PN term with coefficients k5 and k5l are defined in the SEOBNRv2 review document https://dcc.ligo.org/T1400476 */
	logTerms = 1. + eta * coeffs->k0 + eta * log(1. + coeffs->k1 * u
		       + coeffs->k2 * u2 + coeffs->k3 * u3 + coeffs->k4 * u4
			     + coeffs->k5 * u5 + coeffs->k5l * u5 * log(u));
	/* Eq. 5.73 of BB1 */
	deltaU = bulk * logTerms;
	/* Eq. 5.71 of BB1 */
	deltaT = rMag2 * deltaU;
	/* ddeltaU/du */
    /* The log(u) is treated as a constant when taking the derivative wrt u */
	deltaU_u = 2. * (1. / m1PlusetaKK + a2 * u) * logTerms +
		bulk * (eta * (coeffs->k1 + u * (2. * coeffs->k2 + u * (3. * coeffs->k3
									+ u * (4. * coeffs->k4 + 5. * (coeffs->k5 + coeffs->k5l * log(u)) * u)))))
		/ (1. + coeffs->k1 * u + coeffs->k2 * u2 + coeffs->k3 * u3
	      + coeffs->k4 * u4 + (coeffs->k5 + coeffs->k5l * log(u)) * u5);
	deltaU_r = -u2 * deltaU_u;
	/* Eq. 5.38 of BB1 */
	deltaR = deltaT * D;
	if (params.params->tortoise)
		csi = sqrt(deltaT * deltaR) / w2;	/* Eq. 28 of Pan et al.
							 * PRD 81, 084041 (2010) */
	else
		csi = 1.0;

    /* This is A3 of PRD 81, 084041 (2010) explicitly evaluated */
	for (i = 0; i < 3; i++) 
    {
		tmpP[i] = pVec.data[i] - (rVec.data[i] / rMag) * prT * (csi - 1.) / csi;
	}
	/*
	 * Calculate the T-matrix, required to convert P from tortoise to
	 * non-tortoise coordinates, and/or vice-versa. This is given
	 * explicitly in Eq. A3 of 0912.3466
	 */
	for (i = 0; i < 3; i++)
		for (j = 0; j <= i; j++) 
        {
			Tmatrix[i][j] = Tmatrix[j][i] = (rVec.data[i] * rVec.data[j] / rMag2)
				* (csi - 1.);

			invTmatrix[i][j] = invTmatrix[j][i] =
				-(csi - 1) / csi * (rVec.data[i] * rVec.data[j] / rMag2);

			if (i == j) 
            {
				Tmatrix[i][j]++;
				invTmatrix[i][j]++;
			}
		}

     /* This is dcsi/dr: this is needed in the last term of Eq. A5 of PRD 81, 084041 (2010) */
	dcsi = csi * (2. / rMag + deltaU_r / deltaU) + csi * csi * csi
		/ (2. * rMag2 * rMag2 * deltaU * deltaU) * (rMag * (-4. * w2) / D - eobD_r * (w2 * w2));

	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			for (k = 0; k < 3; k++) 
            {
				dTijdXk[i][j][k] =
					(rVec.data[i] * KRONECKER_DELTA(j, k) + KRONECKER_DELTA(i, k) * rVec.data[j])
					* (csi - 1.) / rMag2
					+ rVec.data[i] * rVec.data[j] * rVec.data[k] / rMag2 / rMag * (-2. / rMag * (csi - 1.) + dcsi);
			}

    INT4   updateHCoeffsOld =  params.params->seobCoeffs->updateHCoeffs;
	/* Now calculate derivatives w.r.t. each parameter */
	for (i = 0; i < 3; i++)
    {
		params.varyParam = i;
		params.params->seobCoeffs->updateHCoeffs = 1;
        params.params->tortoise = 2;
        memcpy(tmpValues, params.values, sizeof(tmpValues));
        tmpValues[3] = tmpP[0];
        tmpValues[4] = tmpP[1];
        tmpValues[5] = tmpP[2];
        params.values = tmpValues;
        /* We need derivatives of H wrt to P (and not P^*) */
        /* Note that in the 1st term on the last line of Eq. A5 of PRD 81, 084041 (2010) one needs
         * dH/dX @ fixed P, not P^*, hence the need for what follows  */
        gslStatus = gsl_deriv_central(&F, values[i],
        STEP_SIZE, &tmpDValues[i], &absErr);
        params.values = values;
        params.params->tortoise = 1;
        if (gslStatus != GSL_SUCCESS) 
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "Failure in GSL function");
            return CEV_FAILURE;
        }
    }

    for (i = 3; i < 6; i++) 
    {
        params.varyParam = i;
        params.params->seobCoeffs->updateHCoeffs = 1;
        gslStatus = gsl_deriv_central(&F, values[i],
        STEP_SIZE, &tmpDValues[i], &absErr);
        if (gslStatus != GSL_SUCCESS) 
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "Failure in GSL function");
            return CEV_FAILURE;
        }
    }
    
    for (i = 6; i < 9; i++) 
    {
        params.varyParam = i;
        params.params->seobCoeffs->updateHCoeffs = 1;
        gslStatus = gsl_deriv_central(&F, values[i],
        STEP_SIZE * mass1 * mass1, &tmpDValues[i], &absErr);
        if (gslStatus != GSL_SUCCESS) 
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "Failure in GSL function");
            return CEV_FAILURE;
        }
    }
    for (i = 9; i < 12; i++) 
    {
        params.varyParam = i;
        params.params->seobCoeffs->updateHCoeffs = 1;
        gslStatus = gsl_deriv_central(&F, values[i],
        STEP_SIZE * mass2 * mass2, &tmpDValues[i], &absErr);
        if (gslStatus != GSL_SUCCESS) 
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "Failure in GSL function");
            return CEV_FAILURE;
        }
    }
    params.params->seobCoeffs->updateHCoeffs = updateHCoeffsOld;

	/* Now make the conversion */
	/* rVectorDot */
    // Eq. A4 of PRD 81, 084041 (2010).  Note that dvalues[i] = \dot{X^i} but tmpDValues[j] = dH/dvalues[j]
	for (i = 0; i < 3; i++)
		for (j = 0, dvalues[i] = 0.; j < 3; j++)
			dvalues[i] += tmpDValues[j + 3] * Tmatrix[i][j];

	/* Calculate the orbital angular momentum */
	// Lx = values[1] * values[5] - values[2] * values[4];
	// Ly = values[2] * values[3] - values[0] * values[5];
	// Lz = values[0] * values[4] - values[1] * values[3];
	Lx = vars.LVec[0];
	Ly = vars.LVec[1];
	Lz = vars.LVec[2];
	magL = sqrt(Lx * Lx + Ly * Ly + Lz * Lz);

	Lhatx = Lx / magL;
	Lhaty = Ly / magL;
	Lhatz = Lz / magL;

	/* Calculate the polar data */
	polarDynamics.length = 4;
	polarDynamics.data = polData;

	// r = polarDynamics.data[0] = sqrt(values[0] * values[0] + values[1] * values[1]
	// 		      + values[2] * values[2]);
	r = polarDynamics.data[0] = rMag;
	polarDynamics.data[1] = 0;
    /* Note that poldata[2] and poldata[3] differ from v2 in which one is normalized by polData[0]. Behaviour is reverted. */
	polarDynamics.data[2] = (values[0] * values[3] + values[1] * values[4]
		      + values[2] * values[5]) / polData[0];
	polarDynamics.data[3] = magL;


	/* Now calculate rCrossRdot and omega */
	rCrossV_x = values[1] * dvalues[2] - values[2] * dvalues[1];
	rCrossV_y = values[2] * dvalues[0] - values[0] * dvalues[2];
	rCrossV_z = values[0] * dvalues[1] - values[1] * dvalues[0];

	omega = sqrt(rCrossV_x * rCrossV_x + rCrossV_y * rCrossV_y + rCrossV_z * rCrossV_z) / (r * r);
    REAL8 dr, ncrv;
    ncrv = omega * r;
    dr = (values[0]*dvalues[0] + values[1]*dvalues[1] + values[2]*dvalues[2]) / r;
#if 0
    REAL8 magY, cax, cay, caz, csx, csy, csz;
    REAL8 c1x, c1y, c1z, c2x, c2y, c2z;
    REAL8 m1sq = mass1*mass1;
    REAL8 m2sq = mass2*mass2;
    REAL8 yhat[3] = {0,0,0};
    REAL8 dr, ncrv;
    ncrv = omega * r;
    dr = (values[0]*dvalues[0] + values[1]*dvalues[1] + values[2]*dvalues[2]) / r;
    // Calculate chiAVec, chiSVec
    // cross_product3d(LNhat, rvec, yhat);
    yhat[0] = rCrossV_y*values[2] - rCrossV_z*values[1];
    yhat[1] = rCrossV_z*values[0] - rCrossV_x*values[2];
    yhat[2] = rCrossV_x*values[1] - rCrossV_y*values[0];
    magY = sqrt(yhat[0]*yhat[0] + yhat[1]*yhat[1] + yhat[2]*yhat[2]);
    // magY = sqrt(dvalues[0]*dvalues[0] + dvalues[1]*dvalues[1] + dvalues[2]*dvalues[2]);
    c1x = (s1.data[0] * values[0] + s1.data[1] * values[1] + s1.data[2] * values[2]) / (r * m1sq);
    c1z = (s1.data[0] * rCrossV_x + s1.data[1] * rCrossV_y + s1.data[2] * rCrossV_z) / (r * r * omega * m1sq);
    // c1y = sqrt( (s1.data[0]*s1.data[0] + s1.data[1]*s1.data[1] + s1.data[2]*s1.data[2])/(m1sq*m1sq) - c1x*c1x - c1z*c1z);
    c1y = (s1.data[0] * yhat[0] + s1.data[1] * yhat[1] + s1.data[2] * yhat[2]) / m1sq / magY;

    c2x = (s2.data[0] * values[0] + s2.data[1] * values[1] + s2.data[2] * values[2]) / (r * m2sq);
    c2z = (s2.data[0] * rCrossV_x + s2.data[1] * rCrossV_y + s2.data[2] * rCrossV_z) / (r * r * omega * m2sq);
    // c2y = sqrt( (s2.data[0]*s2.data[0] + s2.data[1]*s2.data[1] + s2.data[2]*s2.data[2])/(m2sq*m2sq) - c2x*c2x - c2z*c2z);
    c2y = (s2.data[0] * yhat[0] + s2.data[1] * yhat[1] + s2.data[2] * yhat[2]) / m2sq / magY;
    cax = 0.5 * (c1x - c2x);
    cay = 0.5 * (c1y - c2y);
    caz = 0.5 * (c1z - c2z);

    csx = 0.5 * (c1x + c2x);
    csy = 0.5 * (c1y + c2y);
    csz = 0.5 * (c1z + c2z);
#endif
    /* Eq. 16 of PRD 89, 084006 (2014): it's S_{1,2}/m_{1,2}^2.LNhat */
	if (SpinAlignedEOBversion == 4)
	{
		s1dotLN = (s1.data[0] * Lhatx + s1.data[1] * Lhaty + s1.data[2] * Lhatz) /
				  (mass1 * mass1);
		s2dotLN = (s2.data[0] * Lhatx + s2.data[1] * Lhaty + s2.data[2] * Lhatz) /
				  (mass2 * mass2);
	}
	else
	{
		s1dotLN = (s1.data[0] * rCrossV_x + s1.data[1] * rCrossV_y + s1.data[2] * rCrossV_z) /
				  (r * r * omega * mass1 * mass1);
		s2dotLN = (s2.data[0] * rCrossV_x + s2.data[1] * rCrossV_y + s2.data[2] * rCrossV_z) /
				  (r * r * omega * mass2 * mass2);
	}

	chiS = 0.5 * (s1dotLN + s2dotLN);
	chiA = 0.5 * (s1dotLN - s2dotLN);
	REAL8 Lhat[3] = {Lhatx, Lhaty, Lhatz};
	REAL8 tempS1_p = inner_product3d(s1.data, Lhat);
	REAL8 tempS2_p = inner_product3d(s2.data, Lhat);
	REAL8 S1_perp[3] = {0, 0, 0};
	REAL8 S2_perp[3] = {0, 0, 0};
	for ( jj = 0; jj < 3; jj++)
	{
		S1_perp[jj] = s1.data[jj] - tempS1_p * Lhat[jj];
		S2_perp[jj] = s2.data[jj] - tempS2_p * Lhat[jj];
	}
	REAL8 sKerr_norm = sqrt(inner_product3d(sKerr.data, sKerr.data));
	REAL8 S_con = 0.0;
	if (sKerr_norm > 1e-6){
		S_con = sKerr.data[0] * Lhat[0] + sKerr.data[1] * Lhat[1] + sKerr.data[2] * Lhat[2];
		S_con /= (1 - 2 * eta);
		S_con += (inner_product3d(S1_perp, sKerr.data) + inner_product3d(S2_perp, sKerr.data)) / sKerr_norm / (1 - 2 * eta) / 2.;
	}
	REAL8 chi = S_con;
	/*
	 * Compute the test-particle limit spin of the deformed-Kerr
	 * background
	 */
    tplspin = (1. - 2. * eta) * chiS + (mass1 - mass2) / (mass1 + mass2) * chiA;

    memcpy(params.params->s1Vec->data, s1norm.data, 3*sizeof(*params.params->s1Vec->data));
    memcpy(params.params->s2Vec->data, s2norm.data, 3*sizeof(*params.params->s2Vec->data));
    memcpy(params.params->sigmaStar->data, sStar.data, 3*sizeof(*params.params->sigmaStar->data));
    memcpy(params.params->sigmaKerr->data, sKerr.data, 3*sizeof(*params.params->sigmaKerr->data));
    
	params.params->a = a;
    /* This is needed because SpinAlignedEOBversion is set to 2 for v3 */
    /* while XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients requires 3 ... */
    SpinAlignedEOBversionForWaveformCoefficients = 4;

    // if ( params.params->use_hm ) SpinAlignedEOBversionForWaveformCoefficients = 451;
    // else SpinAlignedEOBversionForWaveformCoefficients = SpinAlignedEOBversion;
    XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients(params.params->hCoeffs, mass1, mass2, eta, tplspin,
        chiS, chiA, SpinAlignedEOBversionForWaveformCoefficients);
	// if (SpinAlignedEOBversion == 4)
	// {
    XLALSimIMRCalculateSpinPrecEOBHCoeffs_v2(params.params->seobCoeffs, eta, a, chi,
                                                SpinAlignedEOBversion, params.params->hParams);
	// }
	// else
	// {
	// 	XLALSimIMRCalculateSpinPrecEOBHCoeffs(params.params->seobCoeffs, eta, a,
	// 										  SpinAlignedEOBversion);
	// }

	H = XLALSimIMRSpinPrecEOBHamiltonian(eta, &rVec, &pVec, &s1norm, &s2norm,
					     &sKerr, &sStar, params.params->tortoise, params.params->seobCoeffs, params.params->hParams);
    params.params->cache[0] = H;
	H = H * (mass1 + mass2);

	/* Now we have the ingredients to compute the flux */
	memcpy(tmpValues, values, 12 * sizeof(REAL8));
	cartDynamics.data = tmpValues;

    /* Eq. 13 of PRD 86, 024011 (2012) */
    if ( params.params->ignoreflux == 0 ) 
    {
        flux = XLALInspiralPrecSpinFactorizedFlux(&polarDynamics, &cartDynamics,
						  nqcCoeffs, omega, dr, ncrv, params.params, H / (mass1 + mass2), lMax, SpinAlignedEOBversion);
    }
    else if ( params.params->ignoreflux  == 1) 
    {
            flux = 0.;
    }
    else 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Wrong ignorflux option in XLALSpinPrecHcapNumericalDerivative!");
        return CEV_FAILURE;
    }
    params.params->cache[1] = flux;
	/*
	 * Looking at consistency with the non-spinning model, we have to divide the
	 * flux by eta
	 */
	flux = flux / eta;


#if 0
	pDotS1 = pData[0] * s1.data[0] + pVec.data[1] * s1.data[1] + pVec.data[2] * s1.data[2];
	pDotS2 = pVec.data[0] * s2.data[0] + pVec.data[1] * s2.data[1] + pVec.data[2] * s2.data[2];
	rrTerm2 = 8. / 15. * eta * eta * pow(omega, 8. / 3.) / (magL * magL * r) * ((61. + 48. * mass2 / mass1) * pDotS1 + (61. + 48. * mass1 / mass2) * pDotS2);
    REAL8 corrForce[3] = {1., 1., 1.}, cFr, cFf;
    // Calculate pr, prDot
    REAL8 c_pr, c_prDot, c_nDot[3] = {0., 0., 0.}, c_ndtmp;
    c_pr = (values[0] * tmpP[0] + values[1] * tmpP[1] + values[2] * tmpP[2]) / polData[0];
    c_ndtmp = dr / r / r;
    c_nDot[0] = dvalues[0]/r - c_ndtmp;
    c_nDot[1] = dvalues[1]/r - c_ndtmp;
    c_nDot[2] = dvalues[2]/r - c_ndtmp;
    c_prDot = c_nDot[0]*values[3]+c_nDot[1]*values[4]+c_nDot[2]*values[5] - 
        (values[0]*tmpDValues[0] + values[1]*tmpDValues[1] + values[2]*tmpDValues[2])/r;
    CalculateEccCorrectionToFlux(eta, s1dotLN, s2dotLN, r, c_pr, c_prDot, &cFr, &cFf);
    REAL8 c_vec[3] = {1., 1., 1.};
    c_vec[0] = r*cFf*values[3] + (cFr-cFf)*c_pr*values[0]/r;
    c_vec[1] = r*cFf*values[4] + (cFr-cFf)*c_pr*values[1]/r;
    c_vec[2] = r*cFf*values[5] + (cFr-cFf)*c_pr*values[2]/r;
    // print_debug("pr = %f\n", c_pr);
#endif
	/* Now pDot */
	/* Compute the first and second terms in eq. A5 of PRD 81, 084041 (2010) */
# if 1
	for (i = 0; i < 3; i++) 
    {
		for (j = 0, tmpPdotT1[i] = 0.; j < 3; j++)
			tmpPdotT1[i] += -tmpDValues[j] * Tmatrix[i][j];
	}
#else
    // with ecc correction
	for (i = 0; i < 3; i++) 
    {
		for (j = 0, tmpPdotT1[i] = 0.; j < 3; j++)
        {
			tmpPdotT1[i] += -tmpDValues[j] * Tmatrix[i][j];
		    tmpPdotT2[i] = -flux * c_vec[j] * Tmatrix[i][j] / (omega * magL);
        }
	}
#endif
	/* Compute the third term in eq. A5 */
	REAL8	tmpPdotT3T11[3][3][3], tmpPdotT3T12[3][3], tmpPdotT3T2[3];

	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			for (l = 0; l < 3; l++)
				for (k = 0, tmpPdotT3T11[i][j][l] = 0.; k < 3; k++)
					tmpPdotT3T11[i][j][l] += dTijdXk[i][k][j] * invTmatrix[k][l];


	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			for (k = 0, tmpPdotT3T12[i][j] = 0.; k < 3; k++)
				tmpPdotT3T12[i][j] += tmpPdotT3T11[i][j][k] * pData[k];

	for (i = 0; i < 3; i++)
		for (j = 0, tmpPdotT3T2[i] = 0.; j < 3; j++)
			tmpPdotT3T2[i] += tmpDValues[j + 3] * Tmatrix[i][j];

	for (i = 0; i < 3; i++)
		for (j = 0, tmpPdotT3[i] = 0.; j < 3; j++)
			tmpPdotT3[i] += tmpPdotT3T12[i][j] * tmpPdotT3T2[j];

#if 1
	// calculate prTDot, for ecc RR force
	REAL8 prTDot, nDoti;
	prTDot = 0.0;
	// rDotv = 0.0;
	// for (i=0; i<3; i++)
	// 	rDotv += values[i]*dvalues[i];
	for (i=0; i<3; i++)
	{
		nDoti = -values[i]*dr/(r*r) + dvalues[i]/r;
		prTDot += (tmpPdotT1[i] + tmpPdotT3[i])*values[i]/r + nDoti * values[i+3];
	}
	
	params.params->cache[2] = prTDot;
	REAL8 fRRVec[3] = {0}, fRRs1Vec[3] = {0}, fRRs2Vec[3] = {0};
	// if (PREC_FLAG == 2 && params.params->e0 != 0.0)
	if (PREC_FLAG == 2 || (PREC_FLAG == 3 && params.params->e0 != 0.0))
		prec_CalculateRRForce(&vars, eta, prTDot, fRRVec, fRRs1Vec, fRRs2Vec);
	else
	{
		memcpy(fRRVec, values+3, 3*sizeof(REAL8));
	}
	for (i = 0; i < 3; i++)
		tmpPdotT2[i] = -flux * fRRVec[i] / (omega * magL);

#endif

	/* Add them to obtain pDot */
	for (i = 0; i < 3; i++)
		dvalues[i + 3] = tmpPdotT1[i] + tmpPdotT2[i] + tmpPdotT3[i];

    /* Eqs. 11c-11d of PRD 89, 084006 (2014) */
	/* spin1 */
     /* The factor eta is there becasue Hreal is normalized by eta */
	dvalues[6] = eta * (tmpDValues[7] * values[8] - tmpDValues[8] * values[7]) + fRRs1Vec[0];
	dvalues[7] = eta * (tmpDValues[8] * values[6] - tmpDValues[6] * values[8]) + fRRs1Vec[1];
	dvalues[8] = eta * (tmpDValues[6] * values[7] - tmpDValues[7] * values[6]) + fRRs1Vec[2];

	/* spin2 */
    /* The factor eta is there becasue Hreal is normalized by eta */
	dvalues[9] = eta * (tmpDValues[10] * values[11] - tmpDValues[11] * values[10]) + fRRs2Vec[0];
	dvalues[10] = eta * (tmpDValues[11] * values[9] - tmpDValues[9] * values[11]) + fRRs2Vec[1];
	dvalues[11] = eta * (tmpDValues[9] * values[10] - tmpDValues[10] * values[9]) + fRRs2Vec[2];

	/* phase and precessing bit */
	dLx = dvalues[1] * values[5] - dvalues[2] * values[4]
		+ values[1] * dvalues[5] - values[2] * dvalues[4];

	dLy = dvalues[2] * values[3] - dvalues[0] * values[5]
		+ values[2] * dvalues[3] - values[0] * dvalues[5];

	dLz = dvalues[0] * values[4] - dvalues[1] * values[3]
		+ values[0] * dvalues[4] - values[1] * dvalues[3];

	dMagL = (Lx * dLx + Ly * dLy + Lz * dLz) / magL;

	dLhatx = (dLx * magL - Lx * dMagL) / (magL * magL);
	dLhaty = (dLy * magL - Ly * dMagL) / (magL * magL);

	/*
	 * Finn Chernoff convention is used here.
     */
    /* Eqs. 19-20 of PRD 89, 084006 (2014) */
	if (Lhatx == 0.0 && Lhaty == 0.0)
    {
		alphadotcosi = 0.0;
	} else {
		alphadotcosi = Lhatz * (Lhatx * dLhaty - Lhaty * dLhatx) / (Lhatx * Lhatx + Lhaty * Lhaty);
	}
    /* These are ODEs for the phase that enters the h_{lm}: see Eq. 3.11 of PRD 79, 104023 (2009) */
	dvalues[12] = omega - alphadotcosi;
	dvalues[13] = alphadotcosi;


    return CEV_SUCCESS;
}

#undef KRONECKER_DELTA