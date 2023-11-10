/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_PPREC__
#define __INCLUDE_PPREC__

#include "pUtils.h"
#include "pCore.h"

INT SEOBIntegrateDynamics_prec(REAL8Array **dynamics,
                          INT *retLenOut,
                          REAL8Vector *ICvalues,
                          REAL8 EPS_ABS,
                          REAL8 EPS_REL,
                          REAL8 deltaT,
                          REAL8 deltaT_min,
                          REAL8 tstart,
                          REAL8 tend ,
                          SpinEOBParams *seobParams,
                          INT flagConstantSampling);

INT SEOBIntegrateDynamics_prec_inverse(REAL8Array **dynamics,
                          INT *retLenOut,
                          REAL8Vector *ICvalues,
                          REAL8 EPS_ABS,
                          REAL8 EPS_REL,
                          REAL8 deltaT,
                          REAL8 deltaT_min,
                          REAL8 tstart,
                          REAL8 tend ,
                          SpinEOBParams *seobParams,
                          INT flagConstantSampling);
void SEOBConcactInverseDynToAdaSDynPrec(REAL8Array **dyn_out, REAL8Array *dyn_inv, 
        INT *retLen_out, INT retLen_inv);


INT SEOBComputeExtendedSEOBPrecdynamics(SEOBPrecdynamics **seobdynamics,
                                    REAL8Array *dynamics,
                                    INT retLen,
                                    SpinEOBParams *seobParams);

int SEOBInterpolatePrecDynamicsAtTime(
    REAL8Vector **seobdynamics_values, /**<< Output: pointer to vector for
                                          seobdynamics interpolated values */
    REAL8 t,                           /**<< Input: time at which to evaluate */
    SEOBPrecdynamics *seobdynamics         /**<< Input: SEOB dynamics */
);

int SEOBPrecLocateTimePeakOmega(
    REAL8 *tPeakOmega, /**<< Output: time of peak of Omega if found (see inside
        XLALSimLocateOmegaTime for what is returned otherwise) */
    INT *foundPeakOmega, /**<< Output: flag indicating wether tPeakOmega has been found */
    REAL8Array *dynamics,      /**<< Input: array for dynamics */
    SEOBPrecdynamics *seobdynamics,       /**<< Input: SEOB dynamics object */
    UINT retLen,                     /**<< Input: length of dynamics */
    SpinEOBParams *seobParams /**<< SEOB params */
);

INT SEOBPrecCalculateNQCWindowFactorsFromDyn(SEOBPrecdynamics *dyn,
                                         REAL8 thPeak,
                                         REAL8 tr6M,
                                         REAL8 tHiStart,
                                         REAL8 tThresh,
                                         INT is_first,
                                         SpinEOBParams *seobParams);


int SEOBPrecCalculateSphHarmListNQCCoefficientsV4(
    SphHarmListEOBNonQCCoeffs *
        *nqcCoeffsList, /**<< Output: non-quasi-circular coefficients as a list
                           for each mode */
    INT modes[][2],    /**<< Input: array of modes (l,m) */
    UINT nmodes,       /**<< Input: number of modes (l,m) */
    REAL8 tPeakOmega,   /**<< Input: time of peak of Omega */
    SEOBPrecdynamics *seobdynamics,  /**<< Input: SEOB dynamics */
    SpinEOBParams *seobParams,   /**<< Input: SEOB params */
    REAL8Vector *chi1_omegaPeak, /**<< Input: dimensionless spin 1 at peak of
                                    omega in L_N frame */
    REAL8Vector *chi2_omegaPeak /**<< Input: dimensionless spin 2 at peak of
                                    omega in L_N frame */
);

int SEOBPrecCalculateSphHarmListhlmAmpPhase(
    SphHarmListCAmpPhaseSequence **listhlm,               /**<< Output: list of modes for hlm */
    INT modes[][2],            /**<< Input: array of modes (l,m) */
    UINT nmodes,               /**<< Input: number of modes (l,m) */
    SEOBPrecdynamics *seobdynamics, /**<< Input: SEOB dynamics */
    SphHarmListEOBNonQCCoeffs *listnqcCoeffs, /**<< Input: list of NQCs */
    SpinEOBParams *seobParams,                /**<< SEOB params */
    UINT flagNQC /**<< flag to choose wether or not to include NQC */
);

INT SEOBPrecAttachRDToSphHarmListhPlm(
    SphHarmListCAmpPhaseSequence **listhPlm_RDattached, /**<< Output: list of extended modes hlm with RD
                                 attached */
    COMPLEX16Vector **sigmaQNMlm0, /**<< Output: list of QNM complex frequency for modes lm,
                         0th overtone (dimensionless) */
    INT4 modes[][2],  /**<< Input: array of modes (l,m) */
    UINT nmodes,     /**<< Input: number of modes (l,m) */
    REAL8 finalMass,  /**<< Input: final mass computed from fit (scaled by M) */
    REAL8
        finalSpin, /**<< Input: final spin computed from fit (dimensionless) */
    SphHarmListCAmpPhaseSequence *listhPlm, /**<< Input: list of modes hlm */
    REAL8 deltaT,                           /**<< Input: time step */
    UINT retLen,        /**<< Input: length of the input modes and dynamics */
    UINT retLenRDPatch, /**<< Input: length of the ringdown patch */
    REAL8 tAttach,       /**<< Input: time of attachment */
    // REAL8 tStart, /**<< Input: starting time (of the HiS part) */
    REAL8Vector *seobvalues, /**<< Input: vector for dynamics values at time of
                                peak of omega */
    SEOBPrecdynamics *seobdynamics,            /**<< Input: SEOB dynamics */
    SpinEOBParams *seobParams             /**<< SEOB params */
);


int SEOBPrecJoinTimeVector(
    REAL8Vector **tVecPmodes, /**<< Output: vector of times for P-modes
                                 (AdaS+HiS+RDpatch) */
    UINT *retLenPmodes,      /**<< Output: length of output vector of times for
                                 P-modes */
    REAL8 *tJoinHiS,          /**<< Output: first time >= tstartHiS */
    UINT *indexJoinHiS,      /**<< Output: first index >= tstartHiS */
    REAL8 *tJoinAttach,       /**<< Output: first time >= tAttach */
    UINT *indexJoinAttach,   /**<< Output: first index >= tAttach */
    UINT retLenHiSRDpatch,   /**<< Input: length of RD patch to be added at the
                                 end of HiS with the same constant sampling */
    REAL8 deltaTHiS,          /**<< Input: time step for the high sampling */
    REAL8 tstartHiS,          /**<< Input: time of start of HiS */
    REAL8 tAttach,            /**<< Input: time of attachment */
    SEOBPrecdynamics
        *seobdynamicsAdaS, /**<< Input: SEOB dynamics with adaptive-sampling */
    SEOBPrecdynamics
        *seobdynamicsHiS /**<< Input: SEOB dynamics with high-sampling */
);

int SEOBPrecJoinDynamics(
    SEOBPrecdynamics **seobdynamicsJoined,     /**<< Output: pointer to joined dynamics */
    SEOBPrecdynamics *seobdynamics1, /**<< Input: first dynamics */
    SEOBPrecdynamics *seobdynamics2, /**<< Input: second dynamics */
    UINT indexJoin12, /**<< Input: index where to join the two dynamics */
    UINT indexEnd2    /**<< Input: index of the joined dynamics where to stop
                          dynamics 2 (excluded) */
);

int SEOBPrecEulerJ2PFromDynamics(
    REAL8Vector **alphaJ2P, /**<< Output: pointer to vector for alpha J2P */
    REAL8Vector **betaJ2P,  /**<< Output: pointer to vector for beta J2P */
    REAL8Vector **gammaJ2P, /**<< Output: pointer to vector for gamma J2P */
    REAL8Vector *e1J,       /**<< Input: unit Jframe vector e1J */
    REAL8Vector *e2J,       /**<< Input: unit Jframe vector e2J */
    REAL8Vector *e3J,       /**<< Input: unit Jframe vector e3J */
    UINT retLen, /**<< Input: total length of Euler angles data to be allocated
                     (length of P-modes) */
    UINT indexStop, /**<< Input: index where we stop the computation (excluded,
                        index of time of attachment) */
    SEOBPrecdynamics *seobdynamics, /**<<Input: SEOB dynamics (joined AdaS+HiS, up
                                   to tAttach) */
    SpinEOBParams *seobParams  /**<< SEOB params */
);

#endif

