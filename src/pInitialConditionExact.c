/**
 * Writer: Xiaolin.liu
 * xiaolin.liu@mail.bnu.edu.cn
 *
 * This module contains basic functions for  calculation.
 * Functions list:
 * Kernel:
 * 20xx.xx.xx, LOC
 **/

#include "pInitialConditionExact.h"
#include "myLog.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>

typedef struct tagICEParams
{
    REAL8 r;
    REAL8 pphi;
    REAL8 Mfini;
    REAL8 ecc;
    REAL8 newMfini;
    REAL8 newecc;
    SpinEOBParams *core;
} ICEParams;

static int GSLSearchMfiniEcc(const gsl_vector *x, void *params, gsl_vector *f)
{
    REAL8 r, pphi;
    REAL8 ecc, MfMin;
    ICEParams *pms = (ICEParams *)params;
    r = gsl_vector_get(x, 0);
    pphi = gsl_vector_get(x, 1);
    if (r < 6.)
    {
        print_debug("bad\n");
        gsl_vector_set(f, 0, r);
        gsl_vector_set(f, 1, pphi);
        return CEV_SUCCESS;
    }
    MfMin = pms->Mfini;
    ecc = pms->ecc;
    REAL8 recMfMin, rececc;
    print_debug("r = %.16e, pphi = %.16e, e0 = %.16e, Mf = %.16e\n", r, pphi, pms->newecc, pms->newMfini);
    EstimateEquatorialEccentricity(r, pphi, &recMfMin, &rececc, pms->newecc, pms->newMfini, pms->core);
    pms->newMfini = recMfMin;
    pms->newecc = rececc;
    print_debug("tgtecc = %.16e, tgtMfMin = %.16e\n\n", ecc, MfMin);
    gsl_vector_set(f, 0, rececc - ecc);
    gsl_vector_set(f, 1, 1e3 * (recMfMin - MfMin));
    return CEV_SUCCESS;
}

INT SEOBComputeExactEquatorialInitialCondition(REAL8Vector *ICvalues, REAL8 MfMin, REAL8 ecc, SpinEOBParams *seobParams)
{
    if (!seobParams->alignedSpins)
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "This is not a spin-aligned system.");
        return CEV_FAILURE;
    }
    REAL8 eta = seobParams->eta;
    REAL8 m1, m2, mTotal;
    m1 = seobParams->m1;
    m2 = seobParams->m2;
    mTotal = m1 + m2;
    REAL8 mSpin1data[3] = {0., 0., 0.};
    REAL8 mSpin2data[3] = {0., 0., 0.};

    mSpin1data[2] = seobParams->chi1 * m1 * m1;
    mSpin2data[2] = seobParams->chi2 * m2 * m2;

    REAL8 xSph[3] = {0};
    REAL8 pSph[3] = {0};
    REAL8 xCart[3] = {0};
    REAL8 pCart[3] = {0};

    REAL8Vector *tmpInitvals = NULL;
    tmpInitvals = CreateREAL8Vector(14);
    if (EOBInitialConditionsPrec_Conserve(tmpInitvals, m1, m2, MfMin, ecc, 0, mSpin1data, mSpin2data, seobParams) !=
        CEV_SUCCESS)
        return CEV_FAILURE;
    memcpy(xCart, tmpInitvals->data, 3 * sizeof(REAL8));
    memcpy(pCart, tmpInitvals->data + 3, 3 * sizeof(REAL8));
    STRUCTFREE(tmpInitvals, REAL8Vector);
    CartesianToSpherical(xSph, pSph, xCart, pCart);
    // print_debug("xIni = (%.16e, %.16e, %.16e)\n", xSph[0], xSph[1], xSph[2]);
    // print_debug("pIni = (%.16e, %.16e, %.16e)\n", pSph[0], pSph[1], pSph[2]);
    REAL8 r0, pphi0;
    r0 = xSph[0];
    pphi0 = pSph[2];
    STRUCTFREE(tmpInitvals, REAL8Vector);

    ICEParams rootParams;
    const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver *rootSolver = NULL;

    gsl_multiroot_function rootFunction;
    gsl_vector *initValues = NULL;
    gsl_vector *finalValues = NULL;
    int gslStatus;
    const int maxIter = 100;

    memset(&rootParams, 0, sizeof(rootParams));
    rootSolver = gsl_multiroot_fsolver_alloc(T, 2);
    initValues = gsl_vector_calloc(2);

    rootFunction.f = GSLSearchMfiniEcc;
    rootFunction.n = 2;
    rootFunction.params = &rootParams;

    rootParams.r = r0;
    rootParams.pphi = pphi0;
    rootParams.ecc = rootParams.newecc = ecc;
    rootParams.Mfini = rootParams.newMfini = MfMin;

    rootParams.core = seobParams;

    gsl_vector_set(initValues, 0, r0);
    gsl_vector_set(initValues, 1, pphi0);

    gsl_multiroot_fsolver_set(rootSolver, &rootFunction, initValues);
    REAL8 RES = 1.e-3;
    INT i = 0;
    do
    {
        gslStatus = gsl_multiroot_fsolver_iterate(rootSolver);
        if (gslStatus == GSL_ENOPROG || gslStatus == GSL_ENOPROGJ)
        {
            PRINT_LOG_INFO(LOG_ERROR, "NO PROGRESS being made by Spherical orbit root solver\n");
            finalValues = gsl_multiroot_fsolver_f(rootSolver);
            PRINT_LOG_INFO(LOG_ERROR, "Function value here given by the following:\n");
            PRINT_LOG_INFO(LOG_ERROR, " F1 = %.16e, F2 = %.16e\n", gsl_vector_get(finalValues, 0),
                           gsl_vector_get(finalValues, 1));
        }
        else if (gslStatus == GSL_EBADFUNC)
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "Inf or Nan encountered in evaluluation of spherical orbit Eqn");
            gsl_multiroot_fsolver_free(rootSolver);
            gsl_vector_free(initValues);
            return CEV_FAILURE;
        }
        else if (gslStatus != GSL_SUCCESS)
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "Error in GSL iteration function!");
            gsl_multiroot_fsolver_free(rootSolver);
            gsl_vector_free(initValues);
            return CEV_FAILURE;
        }
        gslStatus = gsl_multiroot_test_residual(rootSolver->f, RES);
    } while (gslStatus == GSL_CONTINUE && i <= maxIter);

    if (i > maxIter && gslStatus != GSL_SUCCESS)
    {
        gsl_multiroot_fsolver_free(rootSolver);
        gsl_vector_free(initValues);
        return CEV_FAILURE;
    }
    finalValues = gsl_multiroot_fsolver_root(rootSolver);

    gsl_multiroot_fsolver_free(rootSolver);
    gsl_vector_free(initValues);
    return CEV_SUCCESS;
}

INT EstimateEquatorialEccentricity(REAL8 r, REAL8 pphi, REAL8 *MfOrb, REAL8 *ecc, REAL8 inputfMin, REAL8 inputecc,
                                   SpinEOBParams *core)
{
    REAL8 EPS_REL = 1.0e-9;
    REAL8 EPS_ABS = 1.0e-10;
    INT retLenAdaS, failed = 0, status;
    SEOBdynamics *seobdynamicsAdaS = NULL;
    REAL8Array *dynamicsAdaS = NULL;
    REAL8Vector *ICvalues = NULL;
    ICvalues = CreateREAL8Vector(14);

    // Set initial Conditions
    REAL8 xSph[3] = {r, 0., 0.};
    REAL8 pSph[3] = {0, 0, pphi};
    REAL8 xCart[3] = {0, 0, 0};
    REAL8 pCart[3] = {0, 0, 0};
    SphericalToCartesian(xCart, pCart, xSph, pSph);
    memcpy(ICvalues->data, xCart, sizeof(xCart));
    memcpy(ICvalues->data + 3, pCart, sizeof(pCart));
    memcpy(ICvalues->data + 6, core->s1Vec->data, sizeof(pCart));
    memcpy(ICvalues->data + 9, core->s2Vec->data, sizeof(pCart));

    REAL8 deltaT, deltaT_min, tstartAdaS, tendAdaS;
    REAL8 tthresh;
    deltaT = 0.5;
    deltaT_min = 8.0e-5;
    tstartAdaS = 0.;
    tendAdaS = 12. / (inputfMin * pow(1. - inputecc, 2.));
    tthresh = tendAdaS / 3.;
    print_debug("r = %.16e, pphi = %.16e, deltaT = %.16e, tthresh = %.16e\n", r, pphi, deltaT, tthresh);
    status = SEOBIntegrateDynamics_Conserve(&dynamicsAdaS, &retLenAdaS, ICvalues, EPS_ABS, EPS_REL, deltaT, deltaT_min,
                                            tstartAdaS, tendAdaS, core, core->alignedSpins);
    if (status != CEV_SUCCESS)
    {
        failed = 1;
        goto QUIT;
    }

    status = SEOBComputeExtendedSEOBdynamics_Conserve(&seobdynamicsAdaS, dynamicsAdaS, retLenAdaS, core);
    if (status != CEV_SUCCESS)
    {
        failed = 1;
        goto QUIT;
    }

    REAL8 rmax = 0.0, rmin = seobdynamicsAdaS->polarrVec[0];
    REAL8 phi0 = seobdynamicsAdaS->polarphiVec[0], t0 = seobdynamicsAdaS->tVec[0];
    REAL8 sumDphi = 0.0, avgDphi = 1.;
    REAL8Vector *dpY = CreateREAL8Vector(retLenAdaS - 1);
    REAL8 xavg = 0.0, yavg = 0.0;
    // REAL8 xNoverFac = 1e4;
    INT i, index_thresh = floor(tthresh / deltaT);
    REAL8 *tVec = seobdynamicsAdaS->tVec;
    for (i = 0; i < retLenAdaS; i++)
    {
        if (i > 0)
        {
            sumDphi = sumDphi + (seobdynamicsAdaS->polarphiVec[i] - phi0) / (seobdynamicsAdaS->tVec[i] - t0);
            dpY->data[i] = sumDphi / i;
        }

        if (i > index_thresh)
        {
            xavg = xavg + 1. / tVec[i];
            yavg = yavg + dpY->data[i];
        }
        rmin = GET_MIN(rmin, seobdynamicsAdaS->polarrVec[i]);
        rmax = GET_MAX(rmax, seobdynamicsAdaS->polarrVec[i]);
    }
    xavg = xavg / (retLenAdaS - index_thresh - 1);
    yavg = yavg / (retLenAdaS - index_thresh - 1);

    REAL8 kN = 0.0, kD = 0.0;
    for (i = index_thresh; i < retLenAdaS; i++)
    {
        kN = kN + (1. / tVec[i] - xavg) * (dpY->data[i] - yavg);
        kD = kD + (1. / tVec[i] - xavg) * (1. / tVec[i] - xavg);
    }
    REAL8 kk = kN / kD;
    avgDphi = yavg - kk * xavg;
    *ecc = (1. - (rmin / rmax)) / (1. + (rmin / rmax));
    *MfOrb = avgDphi / CST_2PI;
    // print_debug("kk = %.16e\n", kk);
    // print_debug("retLen = %d, tend = %.16e\n", retLenAdaS,
    // seobdynamicsAdaS->tVec[retLenAdaS-1]); print_debug("dt0 = %.16e, dt1 =
    // %.16e\n",
    //     seobdynamicsAdaS->tVec[1]-seobdynamicsAdaS->tVec[0],
    //     seobdynamicsAdaS->tVec[retLenAdaS-1]-seobdynamicsAdaS->tVec[retLenAdaS-2]);
    // print_debug("rmin = %.16e, rmax = %.16e\n", rmin, rmax);
    print_debug("e0 = %.16e, forb = %.16e\n", *ecc, avgDphi / CST_2PI);

QUIT:
    STRUCTFREE(dynamicsAdaS, REAL8Array);
    STRUCTFREE(seobdynamicsAdaS, SEOBdynamics);
    STRUCTFREE(ICvalues, REAL8Vector);
    STRUCTFREE(dpY, REAL8Vector);
    return CEV_SUCCESS;
}

typedef struct tagSEOBEPIRootParams
{
    REAL8 p0;
    REAL8 e0;
    REAL8 ra;
    REAL8 rp;
    REAL8 values_ra[12];
    REAL8 values_rp[12];
    SpinEOBParams *params;
} SEOBEPIRootParams;

/**
 * @brief Calculate dH/dpr of Hamiltonian in spherical coordinates
 *
 * @param x
 * @param params
 * @return double
 */
static double GSLSpinAlignedHamiltonianDpr(double x, void *params)
{
    HcapSphDeriv2Params *dParams = (HcapSphDeriv2Params *)params;

    REAL8 sphValues[12];
    REAL8 cartValues[12];

    REAL8 dHdpr, dHdx, dHdpy, dHdpz;
    REAL8 r, ptheta, pphi;

    return 0.0;
}

static double TMPSpinAlignedHamiltonian(const REAL8 values[], SpinEOBParams *params)
{
    REAL8 Ham;
    SpinEOBHCoeffs *seobHamCoeffs = params->seobCoeffs;
    REAL8 tmpVec[12];
    REAL8 s1normData[3], s2normData[3], sKerrData[3], sStarData[3];

    /* These are the vectors which will be used in the call to the Hamiltonian */
    REAL8Vector r, p, spin1, spin2, spin1norm, spin2norm;
    REAL8Vector sigmaKerr, sigmaStar;

    int i;
    REAL8 a;
    REAL8 m1 = params->m1;
    REAL8 m2 = params->m2;
    REAL8 mT2 = (m1 + m2) * (m1 + m2);

    memcpy(tmpVec, values, sizeof(tmpVec));

    /* Set the LAL-style vectors to point to the appropriate things */
    r.length = p.length = spin1.length = spin2.length = spin1norm.length = spin2norm.length = 3;
    sigmaKerr.length = sigmaStar.length = 3;
    r.data = tmpVec;
    p.data = tmpVec + 3;
    spin1.data = tmpVec + 6;
    spin2.data = tmpVec + 9;
    spin1norm.data = s1normData;
    spin2norm.data = s2normData;
    sigmaKerr.data = sKerrData;
    sigmaStar.data = sStarData;

    memcpy(s1normData, tmpVec + 6, 3 * sizeof(REAL8));
    memcpy(s2normData, tmpVec + 9, 3 * sizeof(REAL8));

    for (i = 0; i < 3; i++)
    {
        s1normData[i] /= mT2;
        s2normData[i] /= mT2;
    }

    /* Calculate various spin parameters */
    EOBCalculateSigmaKerr(&sigmaKerr, &spin1norm, &spin2norm);
    EOBCalculateSigmaStar(&sigmaStar, m1, m2, &spin1norm, &spin2norm);
    a = sqrt(sigmaKerr.data[0] * sigmaKerr.data[0] + sigmaKerr.data[1] * sigmaKerr.data[1] +
             sigmaKerr.data[2] * sigmaKerr.data[2]);
    // EOBCalculateSpinEOBHamCoeffs(seobHamCoeffs, seobParams->eta, a, params);
    Ham = EOBHamiltonian(params->eta, &r, &p, &spin1norm, &spin2norm, &sigmaKerr, &sigmaStar, params->tortoise,
                         seobHamCoeffs) /
          params->eta;

    return Ham;
}

static int FindSphericalOrbit(const gsl_vector *x, void *params, gsl_vector *f)
{
    SEOBEPIRootParams *rootParams = (SEOBEPIRootParams *)params;

    REAL8 pphi, px_a, px_p;
    REAL8 dHdpx_a, dHdpx_p;
    REAL8 dHdpr_a, dHdpr_p;
    rootParams->values_ra[3] = px_a = gsl_vector_get(x, 0);
    rootParams->values_rp[3] = px_p = gsl_vector_get(x, 1);
    pphi = gsl_vector_get(x, 2);
    rootParams->values_ra[4] = pphi / rootParams->ra;
    rootParams->values_rp[4] = pphi / rootParams->rp;
    /* dHdpx */
    dHdpx_a = XLALSpinHcapNumDerivWRTParam(3, rootParams->values_ra, rootParams->params);
    dHdpx_p = XLALSpinHcapNumDerivWRTParam(3, rootParams->values_rp, rootParams->params);
    if (IS_REAL8_FAIL_NAN(dHdpx_a) || IS_REAL8_FAIL_NAN(dHdpx_p))
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "dHdpx is nan");
        return CEV_FAILURE;
    }

    REAL8 Ham_a, Ham_p;
    Ham_a = TMPSpinAlignedHamiltonian(rootParams->values_ra, rootParams->params);
    Ham_p = TMPSpinAlignedHamiltonian(rootParams->values_rp, rootParams->params);

    // print_debug("dHdpx_a = %.16e, dHdpx_p = %.16e, Ham_a - Ham_p = %.16e\n",
    // dHdpx_a, dHdpx_p, Ham_a - Ham_p); print_debug("Ham_a = %.16e, Ham_p =
    // %.16e\n", Ham_a, Ham_p); print_debug("px_a = %.16e, px_p = %.16e, pphi =
    // %.16e\n", px_a, px_p, pphi);
    gsl_vector_set(f, 0, dHdpx_a);
    gsl_vector_set(f, 1, dHdpx_p);
    gsl_vector_set(f, 2, Ham_a - Ham_p);

    return CEV_SUCCESS;
}

INT EOBInitialConditionsPrec_epi(REAL8Vector *initConds,         /**<< OUTPUT, Initial dynamical variables */
                                 const REAL8 mass1,              /**<< mass 1 */
                                 const REAL8 mass2,              /**<< mass 2 */
                                 const REAL8 p0,                 /**<< Initial semi-latus */
                                 const REAL8 e0, const REAL8 x0, /**<< Inclination */
                                 const REAL8 spin1[],            /**<< Initial spin vector 1 */
                                 const REAL8 spin2[],            /**<< Initial spin vector 2 */
                                 SpinEOBParams *params           /**<< Spin EOB parameters */
)
{
    PRINT_LOG_INFO(LOG_INFO, "Compute initial conditions via semi-latus and eccentricity.");
    /* solving initial conditions */
    /* only spin-aligned cases are supported (x0 == 1, spin1x,y = 0 = spin2x,y)
     */
    if (x0 != 1. || spin1[0] != 0 || spin1[1] != 0 || spin2[0] != 0 || spin2[1] != 0)
    {
        PRINT_LOG_INFO(LOG_WARNING, "The input parameters denote spin-precessing "
                                    "orbits, which is un-supported");
    }
    int i, failed = 0;
    REAL8 mTotal = mass1 + mass2;
    REAL8 eta = mass1 * mass2 / (mTotal * mTotal);

    REAL8 tmpS1[3];
    REAL8 tmpS2[3];
    REAL8 tmpS1Norm[3];
    REAL8 tmpS2Norm[3];

    memcpy(tmpS1, spin1, sizeof(tmpS1));
    memcpy(tmpS2, spin2, sizeof(tmpS2));
    memcpy(tmpS1Norm, spin1, sizeof(tmpS1Norm));
    memcpy(tmpS2Norm, spin2, sizeof(tmpS2Norm));
    for (i = 0; i < 3; i++)
    {
        tmpS1Norm[i] /= mTotal * mTotal;
        tmpS2Norm[i] /= mTotal * mTotal;
    }
    tmpS1[0] = tmpS1[1] = tmpS2[0] = tmpS2[1] = tmpS1Norm[0] = tmpS1Norm[1] = tmpS2Norm[0] = tmpS2Norm[1] = 0.0;
    /* We will need a full values vector for calculating derivs of Hamiltonian */
    REAL8 sphValues[12];
    REAL8 cartValues[12];

    /* compute initial conditions */
    SEOBEPIRootParams rootParams;

    const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver *rootSolver = NULL;

    gsl_multiroot_function rootFunction;
    gsl_vector *initValues = NULL;
    gsl_vector *finalValues = NULL;
    int gslStatus;
    const int maxIter = 100;

    memset(&rootParams, 0, sizeof(rootParams));
    rootParams.params = params;
    rootParams.e0 = e0;
    rootParams.p0 = p0;

    REAL8 sma, ra, rp;
    rootParams.ra = ra = p0 / (1. - e0);
    rootParams.rp = rp = p0 / (1. + e0);
    sma = p0 / (1. - e0 * e0);
    print_debug("ra = %.16e, rp = %.16e, p = %.16e, semi-major_a = %.16e\n", ra, rp, p0, sma);
    rootParams.values_ra[0] = rootParams.ra;
    rootParams.values_rp[0] = rootParams.rp;

    REAL8 pphi0 = sqrt(p0);
    // print_debug("p0 = %.16e, e0 = %.16e, pphi0 = %.16e\n", p0, e0, pphi0);
    rootParams.values_ra[4] = pphi0 / rootParams.ra;
    rootParams.values_rp[4] = pphi0 / rootParams.rp;

    /* Initialise the gsl stuff */
    rootSolver = gsl_multiroot_fsolver_alloc(T, 3);
    if (!rootSolver)
    {
        failed = 1;
        goto QUIT;
    }

    initValues = gsl_vector_calloc(3);
    if (!initValues)
    {
        failed = 1;
        goto QUIT;
    }

    // pr_a, pr_p, pphi
    gsl_vector_set(initValues, 0, 0.0);
    gsl_vector_set(initValues, 1, 0.0);
    gsl_vector_set(initValues, 2, pphi0);

    rootFunction.f = FindSphericalOrbit;
    rootFunction.n = 3;
    rootFunction.params = &rootParams;

    gsl_multiroot_fsolver_set(rootSolver, &rootFunction, initValues);
    i = 0;

    do
    {
        gslStatus = gsl_multiroot_fsolver_iterate(rootSolver);
        if (gslStatus != GSL_SUCCESS)
        {
            print_warning("Error in GSL iteration function!\n");
            gsl_multiroot_fsolver_free(rootSolver);
            gsl_vector_free(initValues);
            return CEV_FAILURE;
        }

        gslStatus = gsl_multiroot_test_residual(rootSolver->f, 1.0e-16);
        i++;
    } while (gslStatus == GSL_CONTINUE && i <= maxIter);

    if (i > maxIter && gslStatus != GSL_SUCCESS)
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Exceed maximun iteration steps, quit");
        failed = 1;
        goto QUIT;
    }

    finalValues = gsl_multiroot_fsolver_root(rootSolver);

    REAL8 qCart[3], pCart[3];
    REAL8 qSph[3], pSph[3];

    memset(qCart, 0, sizeof(qCart));
    memset(pCart, 0, sizeof(pCart));

    qCart[0] = rootParams.ra;
    pCart[0] = gsl_vector_get(finalValues, 0);
    pCart[1] = gsl_vector_get(finalValues, 2) / qCart[0];

    memcpy(initConds->data, qCart, sizeof(qCart));
    memcpy(initConds->data + 3, pCart, sizeof(pCart));
    memcpy(initConds->data + 6, tmpS1Norm, sizeof(tmpS1Norm));
    memcpy(initConds->data + 9, tmpS2Norm, sizeof(tmpS2Norm));

QUIT:
    if (rootSolver)
        gsl_multiroot_fsolver_free(rootSolver);
    if (initValues)
        gsl_vector_free(initValues);
    if (failed)
        return CEV_FAILURE;
    return CEV_SUCCESS;
}
