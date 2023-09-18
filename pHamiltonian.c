/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for calculating EOB real SpinHamiltonian.
* Functions list:
    Calculate_EOBSpinHamiltonian()
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "pHamiltonian.h"
#include "pEnergyFlux.h"
#include "pFactorizedWaveform.h"
#include "newFactorizedWaveformPrec.h"
#include "myLog.h"
#include <gsl/gsl_deriv.h>

// Calibrationv21_Sep8a
// K
static const REAL8 coeff00K = 1.7336;
static const REAL8 coeff01K = -1.62045;
static const REAL8 coeff02K = -1.38086;
static const REAL8 coeff03K = 1.43659;
static const REAL8 coeff10K = 10.2573;
static const REAL8 coeff11K = 2.26831;
static const REAL8 coeff12K = 0;
static const REAL8 coeff13K = -0.426958;
static const REAL8 coeff20K = -126.687;
static const REAL8 coeff21K = 17.3736;
static const REAL8 coeff22K = 6.16466;
static const REAL8 coeff23K = 0;
static const REAL8 coeff30K = 267.788;
static const REAL8 coeff31K = -27.5201;
static const REAL8 coeff32K = 31.1746;
static const REAL8 coeff33K = -59.1658;

// dSO
static const REAL8 coeff00dSO = -44.5324;
static const REAL8 coeff01dSO = 0;
static const REAL8 coeff02dSO = 0;
static const REAL8 coeff03dSO = 66.1987;
static const REAL8 coeff10dSO = 0;
static const REAL8 coeff11dSO = 0;
static const REAL8 coeff12dSO = -343.313;
static const REAL8 coeff13dSO = -568.651;
static const REAL8 coeff20dSO = 0;
static const REAL8 coeff21dSO = 2495.29;
static const REAL8 coeff22dSO = 0;
static const REAL8 coeff23dSO = 147.481;
static const REAL8 coeff30dSO = 0;
static const REAL8 coeff31dSO = 0;
static const REAL8 coeff32dSO = 0;
static const REAL8 coeff33dSO = 0;

// dSS
static const REAL8 coeff00dSS = 6.06807;
static const REAL8 coeff01dSS = 0;
static const REAL8 coeff02dSS = 0;
static const REAL8 coeff03dSS = 0;
static const REAL8 coeff10dSS = -36.0272;
static const REAL8 coeff11dSS = 37.1964;
static const REAL8 coeff12dSS = 0;
static const REAL8 coeff13dSS = -41.0003;
static const REAL8 coeff20dSS = 0;
static const REAL8 coeff21dSS = 0;
static const REAL8 coeff22dSS = -326.325;
static const REAL8 coeff23dSS = 528.511;
static const REAL8 coeff30dSS = 706.958;
static const REAL8 coeff31dSS = 0;
static const REAL8 coeff32dSS = 1161.78;
static const REAL8 coeff33dSS = 0.;

INT EOBCalculateSpinEOBHamCoeffs(SpinEOBHCoeffs *coeffs,
                                 const REAL8 eta,
                                 const REAL8 a,
                                 HyperParams *params)
{
    if (!coeffs)
        return CEV_FAILURE;
    memset (coeffs, 0, sizeof (SpinEOBHCoeffs));


    REAL8 KK, k0, k1, k2, k3, k4, k5, k5l, k1p2, k1p3;
    REAL8 m1PlusEtaKK;

    REAL8 chi = a / (1. - 2. * eta);
    REAL8 eta2 = eta * eta, eta3 = eta2 * eta;
    REAL8 chi2 = chi * chi, chi3 = chi2 * chi;
    static const REAL8 ln2 = 0.6931471805599453094172321214581765680755;
    coeffs->KK = KK =
        coeff00K + coeff01K * chi + coeff02K * chi2 + coeff03K * chi3 +
        coeff10K * eta + coeff11K * eta * chi + coeff12K * eta * chi2 +
        coeff13K * eta * chi3 + coeff20K * eta2 + coeff21K * eta2 * chi +
        coeff22K * eta2 * chi2 + coeff23K * eta2 * chi3 + coeff30K * eta3 +
        coeff31K * eta3 * chi + coeff32K * eta3 * chi2 + coeff33K * eta3 * chi3;
    if (params && params->flagTuning)
        coeffs->KK = params->KK;
    m1PlusEtaKK = -1. + eta * KK;
    /* Eqs. 5.77 - 5.81 of BB1 */
    coeffs->k0 = k0 = KK * (m1PlusEtaKK - 1.);
    coeffs->k1 = k1 = -2. * (k0 + KK) * m1PlusEtaKK;
    k1p2 = k1 * k1;
    k1p3 = k1 * k1p2;
    coeffs->k2 = k2 =
        (k1 * (k1 - 4. * m1PlusEtaKK)) / 2. -
        a * a * k0 * m1PlusEtaKK * m1PlusEtaKK;
    coeffs->k3 = k3 =
        -k1 * k1 * k1 / 3. + k1 * k2 + k1 * k1 * m1PlusEtaKK - 2. * (k2 -
                                                                 m1PlusEtaKK)
        * m1PlusEtaKK - a * a * k1 * m1PlusEtaKK * m1PlusEtaKK;
    coeffs->k4 = k4 =
        (24. * k1 * k1 * k1 * k1 - 96. * k1 * k1 * k2 + 48. * k2 * k2 -
        64. * k1 * k1 * k1 * m1PlusEtaKK + 48. * a * a * (k1 * k1 -
                                                        2. * k2) *
        m1PlusEtaKK * m1PlusEtaKK + 96. * k1 * (k3 + 2. * k2 * m1PlusEtaKK) -
        m1PlusEtaKK * (192. * k3 +
                        m1PlusEtaKK * (-3008. + 123. * CST_PI * CST_PI))) / 96.;
        /* Include eta^2 terms at 4PN from arXiv:1305.4884 */
    coeffs->k5 = k5 = m1PlusEtaKK * m1PlusEtaKK
        * (-4237. / 60. + 128. / 5. * CST_GAMMA +
           2275. * CST_PI * CST_PI / 512. - 1. / 3. * a * a * (k1p3 -
                                                               3. * k1 * k2 +
                                                               3. * k3) -
           (k1p3 * k1p2 - 5. * k1p3 * k2 + 5. * k1 * k2 * k2 +
            5. * k1p2 * k3 - 5. * k2 * k3 -
            5. * k1 * k4) / 5. / m1PlusEtaKK / m1PlusEtaKK + (k1p2 * k1p2 -
                                                              4. * k1p2 * k2 +
                                                              2. * k2 * k2 +
                                                              4. * k1 * k3 -
                                                              4. * k4) / 2. /
           m1PlusEtaKK + 256. / 5. * ln2 + (41. * CST_PI * CST_PI / 32. -
                                                 221. / 6.) * eta);
    coeffs->k5l = k5l = m1PlusEtaKK * m1PlusEtaKK * 64. / 5.;
    coeffs->d1v2 =
            coeff00dSO + coeff01dSO * chi + coeff02dSO * chi2 + coeff03dSO * chi3 +
            coeff10dSO * eta + coeff11dSO * eta * chi + coeff12dSO * eta * chi2 +
            coeff13dSO * eta * chi3 + coeff20dSO * eta2 + coeff21dSO * eta2 * chi +
            coeff22dSO * eta2 * chi2 + coeff23dSO * eta2 * chi3 + coeff30dSO * eta3 +
            coeff31dSO * eta3 * chi + coeff32dSO * eta3 * chi2 + coeff33dSO * eta3 * chi3;
            
            // dSS
    coeffs->dheffSSv2 =
            coeff00dSS + coeff01dSS * chi + coeff02dSS * chi2 + coeff03dSS * chi3 +
            coeff10dSS * eta + coeff11dSS * eta * chi + coeff12dSS * eta * chi2 +
            coeff13dSS * eta * chi3 + coeff20dSS * eta2 + coeff21dSS * eta2 * chi +
            coeff22dSS * eta2 * chi2 + coeff23dSS * eta2 * chi3 + coeff30dSS * eta3 +
            coeff31dSS * eta3 * chi + coeff32dSS * eta3 * chi2 + coeff33dSS * eta3 * chi3;
    if (params && params->flagTuning)
    {
        coeffs->d1v2 = params->dSO;
        coeffs->dheffSSv2 = params->dSS;
    }

    return CEV_SUCCESS;
}

INT XLALSimIMRCalculateSpinPrecEOBHCoeffs_v2(
    SpinEOBHCoeffs *coeffs,           /**<< OUTPUT, EOB parameters including pre-computed coefficients */
    const REAL8 eta,                  /**<< symmetric mass ratio */
    REAL8 a,                          /**<< Normalized deformed Kerr spin */
    REAL8 chi,                        /**<< The augmented spin, with correct aligned-spin limit */
    const UINT SpinAlignedEOBversion, /**<< 4 for SEOBNRv4P; Possible to extend this later later */
    HyperParams *params
)
{

    REAL8 KK, k0, k1, k2, k3, k4, k5, k5l, k1p2, k1p3;
    REAL8 m1PlusEtaKK;

    if (!coeffs)
    {
        return CEV_FAILURE;
    }

    // coeffs->SpinAlignedEOBversion = SpinAlignedEOBversion;
    const int debugPK = 0;
    // if (debugPK)
    // {
    //     PRINT_LOG_INFO(LOG_CRITICAL, "In XLALSimIMRCalculateSpinPrecEOBHCoeffs: SpinAlignedEOBversion = %d,%d\n",
    //                     (int)SpinAlignedEOBversion, (int)coeffs->SpinAlignedEOBversion);
    //     fflush(NULL);
    // }

    static const REAL8 third = 1. / 3.;
    static const REAL8 fifth = 1. / 5.;
    static const REAL8 ln2 = 0.6931471805599453094172321214581765680755; // log(2.)
    coeffs->b3 = 0.;
    coeffs->bb3 = 0.;
    REAL8 eta2 = eta * eta, eta3 = eta2 * eta;
    REAL8 chi2 = chi * chi, chi3 = chi2 * chi;
    coeffs->KK = KK = 0;
    if (SpinAlignedEOBversion == 4)
    {
        // See Eq.(4.8) and Eq.(4.12) in Bohe et al
        coeffs->KK = KK =
            coeff00K + coeff01K * chi + coeff02K * chi2 + coeff03K * chi3 +
            coeff10K * eta + coeff11K * eta * chi + coeff12K * eta * chi2 +
            coeff13K * eta * chi3 + coeff20K * eta2 + coeff21K * eta2 * chi +
            coeff22K * eta2 * chi2 + coeff23K * eta2 * chi3 + coeff30K * eta3 +
            coeff31K * eta3 * chi + coeff32K * eta3 * chi2 + coeff33K * eta3 * chi3;
    }
    if (params && params->flagTuning)
        coeffs->KK = params->KK;

    m1PlusEtaKK = -1. + eta * KK;
    const REAL8 invm1PlusEtaKK = 1. / m1PlusEtaKK;
    /* Eqs.(5.77 - 5.81) of Baruasse and Buonnano PRD 81, 084024 (2010) */
    coeffs->k0 = k0 = KK * (m1PlusEtaKK - 1.);
    coeffs->k1 = k1 = -2. * (k0 + KK) * m1PlusEtaKK;
    k1p2 = k1 * k1;
    k1p3 = k1 * k1p2;
    coeffs->k2 = k2 = (k1 * (k1 - 4. * m1PlusEtaKK)) * 0.5 - a * a * k0 * m1PlusEtaKK * m1PlusEtaKK;
    coeffs->k3 = k3 = -(k1 * k1) * k1 * third + k1 * k2 + (k1 * k1) * m1PlusEtaKK - 2. * (k2 - m1PlusEtaKK) * m1PlusEtaKK - a * a * k1 * (m1PlusEtaKK * m1PlusEtaKK);
    coeffs->k4 = k4 = ((24. / 96.) * (k1 * k1) * (k1 * k1) - (96. / 96.) * (k1 * k1) * k2 + (48. / 96.) * k2 * k2 - (64. / 96.) * (k1 * k1) * k1 * m1PlusEtaKK + (48. / 96.) * (a * a) * (k1 * k1 - 2. * k2) * (m1PlusEtaKK * m1PlusEtaKK) +
                        (96. / 96.) * k1 * (k3 + 2. * k2 * m1PlusEtaKK) - m1PlusEtaKK * ((192. / 96.) * k3 + m1PlusEtaKK * (-(3008. / 96.) + (123. / 96.) * CST_PI * CST_PI)));
    if (SpinAlignedEOBversion == 4)
    {
        // Look at Eq.(A2C) in Steinhoff et al. The last term from Eq.(2.3) in Bohe et al.
        coeffs->k5 = k5 = m1PlusEtaKK * m1PlusEtaKK * (-4237. / 60. + 128. / 5. * CST_GAMMA + 2275. * CST_PI * CST_PI / 512. - third * (a * a) * (k1p3 - 3. * (k1 * k2) + 3. * k3) - ((k1p3 * k1p2) - 5. * (k1p3 * k2) + 5. * k1 * k2 * k2 + 5. * k1p2 * k3 - 5. * k2 * k3 - 5. * k1 * k4) * fifth * invm1PlusEtaKK * invm1PlusEtaKK + ((k1p2 * k1p2) - 4. * (k1p2 * k2) + 2. * k2 * k2 + 4. * k1 * k3 - 4. * k4) * 0.5 * invm1PlusEtaKK + (256. / 5.) * ln2 + (41. * CST_PI * CST_PI / 32. - 221. / 6.) * eta);
        // This is the first term in the brackets in Eq.(A2C) in Steinhoff et al.
        coeffs->k5l = k5l = (m1PlusEtaKK * m1PlusEtaKK) * (64. / 5.);
    }

    /* Now calibrated parameters for spin models */
    if (SpinAlignedEOBversion == 4)
    {

        coeffs->d1 = 0.;
        coeffs->dheffSS = 0.;
        // dSO: Eq.(4.13) in Bohe et al
        coeffs->d1v2 =
            coeff00dSO + coeff01dSO * chi + coeff02dSO * chi2 + coeff03dSO * chi3 +
            coeff10dSO * eta + coeff11dSO * eta * chi + coeff12dSO * eta * chi2 +
            coeff13dSO * eta * chi3 + coeff20dSO * eta2 + coeff21dSO * eta2 * chi +
            coeff22dSO * eta2 * chi2 + coeff23dSO * eta2 * chi3 + coeff30dSO * eta3 +
            coeff31dSO * eta3 * chi + coeff32dSO * eta3 * chi2 + coeff33dSO * eta3 * chi3;

        // dSS: Eq.(4.14) in Bohe et al
        coeffs->dheffSSv2 =
            coeff00dSS + coeff01dSS * chi + coeff02dSS * chi2 + coeff03dSS * chi3 +
            coeff10dSS * eta + coeff11dSS * eta * chi + coeff12dSS * eta * chi2 +
            coeff13dSS * eta * chi3 + coeff20dSS * eta2 + coeff21dSS * eta2 * chi +
            coeff22dSS * eta2 * chi2 + coeff23dSS * eta2 * chi3 + coeff30dSS * eta3 +
            coeff31dSS * eta3 * chi + coeff32dSS * eta3 * chi2 + coeff33dSS * eta3 * chi3;
        //          printf("dSO %.16e, dSS %.16e\n", coeffs->d1v2,coeffs->dheffSSv2);
    }
    if (params && params->flagTuning)
    {
        coeffs->d1v2 = params->dSO;
        coeffs->dheffSSv2 = params->dSS;
    }

// print_debug("Hparams:\n");
// print_err("a = %.16e, eta = %.16e, chi = %.16e\n", a, eta, chi);
// print_err("K = %.16e, k0 = %.16e, k1 = %.16e\n", KK, k0, k1);
// print_err("k2 = %.16e, k3 = %.16e, k4 = %.16e\n", k2, k3, k4);
// print_err("k5 = %.16e, k5l = %.16e, d1v2 = %.16e, dheffSSv2 = %.16e\n", k5, k5l, coeffs->d1v2, coeffs->dheffSSv2);

    return CEV_SUCCESS;
}

void CalculateSpinEOBHSACoeffs(REAL8 m1, REAL8 m2, REAL8 s1z, REAL8 s2z, SpinEOBHSACoeffs *coeffs)
{
    REAL8 mtot = m1 + m2;
    REAL8 eta = m1*m2/mtot/mtot;
    REAL8 s1N = s1z * m1*m1 / mtot / mtot;
    REAL8 s2N = s2z * m2*m2 / mtot / mtot;
    coeffs->eta = eta;
    REAL8 sigmaStar, sigmaKerr;
    sigmaStar = coeffs->sigmaStar = m2*s1N/m1 + m1*s2N/m2;
    sigmaKerr = coeffs->sigmaKerr = s1N + s2N;
    if(sigmaKerr < 0.0)
        coeffs->sign = -1.;
    else
        coeffs->sign = 1.;
    REAL8 a;
    a = coeffs->a = fabs(coeffs->sigmaKerr);
    coeffs->a2 = a*a;
    REAL8 chi, chi2, chi3;
    REAL8 eta2, eta3;
    chi = sigmaKerr / (1. - eta*2.);
    chi2 = chi*chi;
    chi3 = chi2*chi;
    eta2 = eta*eta;
    eta3 = eta2*eta;
    REAL8 KK, m1PlusEtaKK;
    KK = 1.7336 + (-1.62045)*chi + (-1.38086)*chi2 + (1.43659)*
        chi3 + (10.2573)*eta + (2.26831)*eta*chi + (-0.426958)*eta*
        chi3 + (-126.687)*eta2 + (17.3736)*eta2*chi + (6.16466)*eta2*
        chi2 + (267.788)*eta3 + (-27.5201)*eta3*chi + (31.1746)*eta3*
        chi2 + (-59.1658)*eta3*chi3;
    REAL8 k0, k1, k2, k3, k4;
    m1PlusEtaKK = -1. + eta*KK;
    coeffs->k0 = k0 = KK*(m1PlusEtaKK - 1.);
    coeffs->k1 = k1 = -2.*(k0 + KK)*m1PlusEtaKK;
    coeffs->k2 = k2 = (k1*(k1 - 4.*m1PlusEtaKK))*0.5 - a*a*k0*m1PlusEtaKK*m1PlusEtaKK;
    coeffs->k3 = k3 = -(k1*k1)*k1*(1./3.) + k1*k2 + (k1*k1)*m1PlusEtaKK - 
        2.*(k2 - m1PlusEtaKK)*m1PlusEtaKK - 
        a*a*k1*(m1PlusEtaKK*m1PlusEtaKK);
    coeffs->k4 = k4 = ((24./96.)*(k1*k1)*(k1*k1) - (96./96.)*(k1*k1)*k2 + (48./96.)*
            k2*k2 - (64./96.)*(k1*k1)*k1*m1PlusEtaKK + (48./96.)*(a*a)*(k1*k1 - 2.*k2)*(m1PlusEtaKK*
            m1PlusEtaKK) + (96./96.)*k1*(k3 + 2.*k2*m1PlusEtaKK) - 
            m1PlusEtaKK*((192./96.)*k3 + m1PlusEtaKK*(-(3008./96.) + (123./96.)*CST_PI*CST_PI)));

    REAL8 k1p2 = k1*k1;
    REAL8 k1p3 = k1*k1p2;
    REAL8 invm1PlusEtaKK = 1./m1PlusEtaKK;
    coeffs->k5 = m1PlusEtaKK*
        m1PlusEtaKK*(-4237./60. + 128./5.*CST_GAMMA + 
        2275.*CST_PI*CST_PI/
            512. - (1./3.)*(a*a)*(k1p3 - 3.*(k1*k2) + 
            3.*k3) - ((k1p3*k1p2) - 5.*(k1p3*k2) + 5.*k1*k2*k2 + 
            5.*k1p2*k3 - 5.*k2*k3 - 5.*k1*k4)*(1./5.)*invm1PlusEtaKK*
            invm1PlusEtaKK + ((k1p2*k1p2) - 4.*(k1p2*k2) + 2.*k2*k2 + 
            4.*k1*k3 - 4.*k4)*0.5*invm1PlusEtaKK + (256./5.)*
            CST_LN2 + (41.*CST_PI*CST_PI/32. - 221./6.)*eta);
    coeffs->k5l = (m1PlusEtaKK*m1PlusEtaKK)*(64./5.);
    REAL8 d1v2 = -44.5324 + 
        66.1987*chi3 + (-343.313)*eta*chi2 + (-568.651)*eta*
        chi3 + (2495.29)*eta2*chi + (147.481)*eta2*chi3;
    REAL8 dheffSSv2 = 
        6.06807 + (-36.0272)*eta + 
        37.1964*eta*chi + (-41.0003)*eta*chi3 + (-326.325)*eta2*
        chi2 + (528.511)*eta2*chi3 + (706.958)*eta3 + (1161.78)*eta3*
        chi2;
// print_debug("k1 = %.16e, k2 = %.16e, k3 = %.16e, k4 = %.16e, k5 = %.16e, k5l = %.16e\n",
//         k1, k2, k3, k4, coeffs->k5, coeffs->k5l);
    coeffs->b3 = 0;
    coeffs->bb3 = 0;
    coeffs->invm1PlusEtaKK = invm1PlusEtaKK;

    coeffs->PQ4 = 2.*eta*(4.-3.*eta);
    coeffs->Pds0u = eta*(7.*sigmaStar - 4.*sigmaKerr)/6.;
    coeffs->Pds0p2 = eta * (3.*sigmaKerr + 4.*sigmaStar)/12.;
    coeffs->Pds0pn2 = -eta*(6.*sigmaKerr + 5.*sigmaStar)/2.;

    coeffs->PdsSu2 = eta*(353. - 27.*eta)/36.;
    coeffs->PdsSup2 = (-103. + 60.*eta)*eta / 36.;
    coeffs->PdsSpn4 = 5.*eta*eta;
    coeffs->PdsSp4 = (-23.-3.*eta)*eta/72.;
    coeffs->PdsSupn2 = (47. - 54.*eta)*eta/12.;
    coeffs->PdsSp2pn2 = (16. - 21.*eta)*eta/12.;
    
    coeffs->PdsKu2 = (-56. - 21.*eta)*eta/9.;
    coeffs->PdsKpn4 = 5.*27*eta*eta/24.;
    coeffs->PdsKp4 = -45.*eta/144.;
    coeffs->PdsKup2 = (-109.+51.*eta)*eta/36.;
    coeffs->PdsKp2pn2 = (6. - 39.*eta)*eta/24.;
    coeffs->PdsKupn2 = (-16.-147.*eta)*eta/24.;
    coeffs->PdsKu3 = d1v2 * eta;
    coeffs->Peffss = dheffSSv2 * eta * (s1N*s1N + s2N*s2N);
// print_debug("PdsKu2 = %.16e, PdsKpn4 = %.16e, PdsKup2 = %.16e, PdsKp2pn2 = %.16e, PdsKupn2 = %.16e, PdsKu3 = %.16e\n",
//         coeffs->PdsKu2, coeffs->PdsKpn4, coeffs->PdsKup2, coeffs->PdsKp2pn2, coeffs->PdsKupn2, coeffs->PdsKu3);
    coeffs->D2 = 6.*eta;
    coeffs->D3 = 2.*(26.-3.*eta)*eta;
    if (IS_DEBUG)
    {
        print_debug("eta = %.16e, chi = %.16e, KK = %.16e\n", eta, chi, KK);
    }
    return;
}

REAL8 EOBSAHamiltonian(REAL8 r, REAL8 prT, REAL8 pphi, SpinEOBHSACoeffs *coeffs, REAL8 *csi)
{
    REAL8 r2 = r*r;
    REAL8 u = 1./r, u2, u3, u4, u5, logu;
    REAL8 deltaU, deltaU0, deltaULog, deltaU_u;
    REAL8 deltaT, deltaT_r, DD, deltaR, B, CC, w2, Lambda, Lambda_r, invLambda;
    REAL8 invxia, pr, wfd, wfd_r;
    REAL8 alpha, beta, gammar, gammaphi;
    REAL8 pnBar2, pBar2, Xi;
    REAL8 deltaSigma0, deltaSigmaS, deltaSigmaK;
    REAL8 HNS, SStar, wr, BR, nur, sqQ, QTilde;
    REAL8 HTwr, HSOL, HSONL, HSS, HSSeff, HS, Heff, Hreal;
    REAL8 pf = coeffs->sign * pphi;
    u2 = u*u;
    u3 = u2*u;
    u4 = u3*u;
    u5 = u4*u;
    deltaU0 = coeffs->invm1PlusEtaKK * (coeffs->invm1PlusEtaKK + 2.*u) + coeffs->a2 * u2;
    // logu = log2(u)*CST_INV_LOG2E;
    // logu = log(u);
    const REAL8 invlog_2e = 0.69314718055994530941723212145817656807550013436026;
    logu = log2(u)*invlog_2e;
    // const REAL8 logarg = coeffs->k1*u + coeffs->k2*u2 + coeffs->k3*u3 + coeffs->k4*u4
    //                                             + coeffs->k5*u5 + coeffs->k5l*u5*logu;
    const REAL8 logarg = u*(coeffs->k1 + u*(coeffs->k2 + u*(coeffs->k3 + u*(coeffs->k4 + u*(coeffs->k5 + coeffs->k5l*logu)))));
    // deltaULog = 1. + coeffs->eta*coeffs->k0 + coeffs->eta*log1p(fabs(1. + logarg) - 1.);
    deltaULog = 1. + coeffs->eta*coeffs->k0 + coeffs->eta*log(1. + logarg);
    deltaU = deltaU0 * deltaULog;
    // deltaU_u = 2.*(coeffs->invm1PlusEtaKK + coeffs->a2*u) * deltaULog + 
    //     deltaU0 * coeffs->eta * (coeffs->k1 + 2.*coeffs->k2*u + 3.*coeffs->k3*u2 + 4.*coeffs->k4*u3 + 5.*(coeffs->k5 + coeffs->k5l*logu)*u4) / (1.+logarg);
    deltaU_u = 2.*(coeffs->invm1PlusEtaKK + coeffs->a2*u) * deltaULog + 
        deltaU0 * coeffs->eta * (coeffs->k1 + 
            u*(2.*coeffs->k2 + 
            u*(3.*coeffs->k3 + 
            u*(4.*coeffs->k4 + 
            u*5.*(coeffs->k5 + coeffs->k5l*logu))))) / (1.+logarg);

    deltaT = r2 * deltaU;
    deltaT_r = 2.*r*deltaU - deltaU_u;
    DD = 1. + log(1. + coeffs->D2 * u2 + coeffs->D3 * u3);
    // DD = 1. + log1p(coeffs->D2 * u2 + coeffs->D3 * u3);
    deltaR = deltaT * DD;
    B = sqrt(deltaT);
    CC = sqrt(deltaR);
    w2 = r2 + coeffs->a2;
    Lambda = w2*w2 - coeffs->a2*deltaT;
    Lambda_r = 4.*r*w2 - coeffs->a2*deltaT_r;
    invLambda = 1./Lambda;
    invxia = w2 / (B*CC);
    *csi = B*CC / w2;
    pr = prT * invxia;
    wfd = 2.*coeffs->a*r + coeffs->b3 / r + coeffs->bb3 / r;
    wfd_r = 2.*coeffs->a - coeffs->b3 * u2 - coeffs->bb3*u2;
    alpha = r * B * sqrt(invLambda);
    beta = wfd * invLambda;
    gammar = deltaR * u2;
    gammaphi = r2 * invLambda;
    pnBar2 = gammar * pr*pr;
    pBar2 = pnBar2 + gammaphi * pf*pf;
    Xi = 1. + coeffs->PQ4 * prT * prT *prT *prT * u2 + pBar2;
    HNS = pf * beta + alpha * sqrt(Xi);
    // if(IS_DEBUG){
    // print_debug( "term 1 in Hns: %.16e\n",  coeffs->PQ4 * prT * prT *prT *prT * u2 );
    // print_debug( "term 2 in Hns: %.16e\n", 0.0 );
    // print_debug( "term 3 in Hns = %.16e\n", gammaphi * pf*pf );
    // print_debug( "term 4 in Hns = %.16e\n", pnBar2 );
    // print_debug( "term 5 in Hns = %.16e\n", 1./(alpha*alpha) );
    // print_debug( "term 6 in Hns = %.16e\n", pf );}
    deltaSigma0 = coeffs->Pds0u * u + coeffs->Pds0p2 * pBar2 + coeffs->Pds0pn2 * pnBar2;
    deltaSigmaS = coeffs->PdsSu2 * u2 + coeffs->PdsSup2 * u*pBar2 + coeffs->PdsSpn4*pnBar2*pnBar2 + coeffs->PdsSp4*pBar2*pBar2 + 
        coeffs->PdsSupn2*u*pnBar2+coeffs->PdsSp2pn2*pBar2*pnBar2;
    deltaSigmaK = coeffs->PdsKu2*u2 + coeffs->PdsKup2*u*pBar2 + coeffs->PdsKpn4*pnBar2*pnBar2 + coeffs->PdsKp4*pBar2*pBar2 + 
        coeffs->PdsKupn2 * u*pnBar2 + coeffs->PdsKp2pn2*pBar2*pnBar2 + coeffs->PdsKu3 * u3;
    SStar = coeffs->sign*(coeffs->sigmaStar + deltaSigma0 + deltaSigmaS * coeffs->sigmaStar + deltaSigmaK * coeffs->sigmaKerr);
    wr = (-Lambda_r * wfd + Lambda * wfd_r) * invLambda*invLambda;
    BR = -sqrt(1./DD) + deltaT_r/(2.*B);
    nur = u + 0.5*w2*(-4.*r*deltaT + w2*deltaT_r)*invLambda/deltaT;
    sqQ = sqrt(1. + pBar2);
    QTilde = sqQ * (1. + sqQ);
    REAL8 alpha2 = alpha*alpha;
    HTwr = 0.5 * CC * SStar * (r2 * alpha2*pf*pf + deltaT*(r2*QTilde - pr*pr*deltaR))*u3 / (B*QTilde*alpha);
    HSOL = alpha2*u*(-B+r*alpha)*pf*SStar/(deltaT*sqQ);
    HSONL = alpha2*SStar*CC*pf*(-BR*(1.+sqQ) + B*nur*(1.+2.*sqQ))*u/(deltaT*QTilde);
    HSS = -0.5*u3*SStar*SStar;
    HS = beta*SStar + wr*HTwr + HSOL + HSONL;
    HSSeff = coeffs->Peffss * u4;
    Heff = HNS + HS + HSS + HSSeff;
    Hreal = sqrt(1. + 2.*coeffs->eta*(Heff-1.));

if (IS_DEBUG)
{
    print_debug("u = %.16e, \n", u);
    print_debug("a = %.16e, eta = %.16e\n", coeffs->a, coeffs->eta);
    print_debug("logarg = %.16e, deltaU0 = %.16e, deltaULog = %.16e\n", logarg, deltaU0, deltaULog);
    print_debug("deltaU = %.16e, deltaU_u = %.16e\n", deltaU, deltaU_u);
    print_debug("deltaT = %.16e, deltaT_r = %.16e\n", deltaT, deltaT_r);
    print_debug("DD = %.16e, deltaR = %.16e\n", DD, deltaR);
    print_debug("wfd = %.16e, wfd_r = %.16e\n", wfd, wfd_r);
    print_debug("pr = %.16e, pnBar2 = %.16e, pBar2 = %.16e\n", pr, pnBar2, pBar2);
    print_debug("alpha = %.16e, beta = %.16e, gammar = %.16e, gammaphi = %.16e\n", alpha, beta, gammar, gammaphi);
    print_debug("deltasigma0 = %.16e, deltasigmaS = %.16e, deltasigmaK = %.16e\n", deltaSigma0, deltaSigmaS, deltaSigmaK);
    print_debug("BR = %.16e, nur = %.16e, wr = %.16e, SStar = %.16e\n", BR, nur, wr, SStar);
    print_debug("Hns = %.16e, HTwr = %.16e, HSOL = %.16e, HSONL = %.16e\n", HNS, HTwr, HSOL, HSONL);
    print_debug("HeffSS = %.16e, HS = %.16e, Heff = %.16e\n", HSSeff, HS, Heff);
    print_debug("Hreal = %.16e\n", Hreal);

}
    return Hreal;
}

REAL8 EOBHamiltonian(const REAL8 eta,
                     REAL8Vector * x,
                     REAL8Vector * p,
                     REAL8Vector * s1Vec,
                     REAL8Vector * s2Vec,
                     REAL8Vector * sigmaKerr,
                     REAL8Vector * sigmaStar,
                     INT    tortoise,
                     SpinEOBHCoeffs *coeffs)
{
    REAL8 r, r2, nx, ny, nz;
    REAL8 sKerr_x, sKerr_y, sKerr_z, a, a2;
    REAL8 sStar_x, sStar_y, sStar_z;
    REAL8 e3_x, e3_y, e3_z;
    REAL8 costheta;        /* Cosine of angle between Skerr and r */
    REAL8 xi2, xi_x, xi_y, xi_z;    /* Cross product of unit vectors in direction of Skerr and r */
    REAL8 vx, vy, vz, pxir, pvr, pn, prT, pr, pf, ptheta2;    /*prT is the tortoise pr */
    REAL8 w2, rho2;
    REAL8 u, u2, u3, u4, u5;
    REAL8 bulk, deltaT, deltaR, Lambda;
    REAL8 D, qq, ww, B, w, MU, nu, BR, wr, nur, mur;
    REAL8 wcos, nucos, mucos, ww_r, Lambda_r;
    REAL8 logTerms, deltaU, deltaU_u, Q, deltaT_r, pn2, pp;
    REAL8 deltaSigmaStar_x, deltaSigmaStar_y, deltaSigmaStar_z;
    REAL8 sx, sy, sz, sxi, sv, sn, s3;
    REAL8 H, Hns, Hs, Hss, Hreal, Hwcos, Hwr, HSOL, HSONL;
    REAL8 m1PlusetaKK;
    
    /* Terms which come into the 3.5PN mapping of the spins */
    //REAL8 aaa, bbb, a13P5, a23P5, a33P5, b13P5, b23P5, b33P5;
    REAL8 sMultiplier1, sMultiplier2;
    
    /*Temporary p vector which we will make non-tortoise */
    REAL8 tmpP[3];
    
    REAL8 csi;
    
    /* Spin gauge parameters. (YP) simplified, since both are zero. */
    // static const double aa=0., bb=0.;
    
    //printf( "In Hamiltonian:\n" );
    //printf( "x = %.16e\t%.16e\t%.16e\n", x->data[0], x->data[1], x->data[2] );
    //printf( "p = %.16e\t%.16e\t%.16e\n", p->data[0], p->data[1], p->data[2] );
    
    r2 =
    x->data[0] * x->data[0] + x->data[1] * x->data[1] +
    x->data[2] * x->data[2];
    r = sqrt (r2);
    nx = x->data[0] / r;
    ny = x->data[1] / r;
    nz = x->data[2] / r;
    
    sKerr_x = sigmaKerr->data[0];
    sKerr_y = sigmaKerr->data[1];
    sKerr_z = sigmaKerr->data[2];
    
    sStar_x = sigmaStar->data[0];
    sStar_y = sigmaStar->data[1];
    sStar_z = sigmaStar->data[2];
    
    a2 = sKerr_x * sKerr_x + sKerr_y * sKerr_y + sKerr_z * sKerr_z;
    a = GET_SQRT (a2);
    
    if (a != 0.)
    {
        e3_x = sKerr_x / a;
        e3_y = sKerr_y / a;
        e3_z = sKerr_z / a;
    }
    else
    {
        e3_x = 0.;
        e3_y = 0.;
        e3_z = 1.;
    }
    
    costheta = e3_x * nx + e3_y * ny + e3_z * nz;
    
    xi2 = 1. - costheta * costheta;
    
    xi_x = -e3_z * ny + e3_y * nz;
    xi_y = e3_z * nx - e3_x * nz;
    xi_z = -e3_y * nx + e3_x * ny;
    
    vx = -nz * xi_y + ny * xi_z;
    vy = nz * xi_x - nx * xi_z;
    vz = -ny * xi_x + nx * xi_y;
    
    w2 = r2 + a2;
    rho2 = r2 + a2 * costheta * costheta;
    
    u = 1. / r;
    u2 = u * u;
    u3 = u2 * u;
    u4 = u2 * u2;
    u5 = u4 * u;
    
    //printf( "KK = %.16e\n", coeffs->KK );
    m1PlusetaKK = -1. + eta * coeffs->KK;
    /* Eq. 5.75 of BB1 */
    bulk = 1. / (m1PlusetaKK * m1PlusetaKK) + (2. * u) / m1PlusetaKK + a2 * u2;
    /* Eq. 5.73 of BB1 */
    logTerms =
    1. + eta * coeffs->k0 + eta * GET_LN (1. + coeffs->k1 * u + coeffs->k2 * u2 +
                                       coeffs->k3 * u3 + coeffs->k4 * u4 +
                                       coeffs->k5 * u5 +
                                       coeffs->k5l * u5 * GET_LN (u));
    //printf( "bulk = %.16e, logTerms = %.16e\n", bulk, logTerms );
    /* Eq. 5.73 of BB1 */
    deltaU = bulk * logTerms;
    
    /* Eq. 5.71 of BB1 */
    deltaT = r2 * deltaU;
    /* ddeltaU/du */
    deltaU_u = 2. * (1. / m1PlusetaKK + a2 * u) * logTerms +
    bulk * (eta *
            (coeffs->k1 +
             u * (2. * coeffs->k2 +
                  u * (3. * coeffs->k3 +
                       u * (4. * coeffs->k4 +
                            5. * (coeffs->k5 +
                                  coeffs->k5l * GET_LN (u)) * u))))) / (1. +
                                                                     coeffs->
                                                                     k1 * u +
                                                                     coeffs->
                                                                     k2 * u2 +
                                                                     coeffs->
                                                                     k3 * u3 +
                                                                     coeffs->
                                                                     k4 * u4 +
                                                                     (coeffs->
                                                                      k5 +
                                                                      coeffs->
                                                                      k5l *
                                                                      GET_LN (u))
                                                                     * u5);
    /* ddeltaT/dr */
    deltaT_r = 2. * r * deltaU - deltaU_u;
    /* Eq. 5.39 of BB1 */
    Lambda = w2 * w2 - a2 * deltaT * xi2;
    /* Eq. 5.83 of BB1, inverse */
    D = 1. + log (1. + 6. * eta * u2 + 2. * (26. - 3. * eta) * eta * u3);
    /* Eq. 5.38 of BB1 */
    deltaR = deltaT * D;
    /* See Hns below, Eq. 4.34 of Damour et al. PRD 62, 084011 (2000) */
    qq = 2. * eta * (4. - 3. * eta);
    /* See Hns below. In Sec. II D of BB2 b3 and bb3 coeffs are chosen to be zero. */
    ww = 2. * a * r + coeffs->b3 * eta * a2 * a * u + coeffs->bb3 * eta * a * u;
    
    /* We need to transform the momentum to get the tortoise co-ord */
    if (tortoise)
    {
        csi = GET_SQRT (deltaT * deltaR) / w2;    /* Eq. 28 of Pan et al. PRD 81, 084041 (2010) */
    }
    else
    {
        csi = 1.0;
    }
    //printf( "csi(miami) = %.16e\n", csi );
    
    prT = p->data[0] * nx + p->data[1] * ny + p->data[2] * nz;
    /* p->data is BL momentum vector; tmpP is tortoise momentum vector */
    tmpP[0] = p->data[0] - nx * prT * (csi - 1.) / csi;
    tmpP[1] = p->data[1] - ny * prT * (csi - 1.) / csi;
    tmpP[2] = p->data[2] - nz * prT * (csi - 1.) / csi;
    
    pxir = (tmpP[0] * xi_x + tmpP[1] * xi_y + tmpP[2] * xi_z) * r;
    pvr = (tmpP[0] * vx + tmpP[1] * vy + tmpP[2] * vz) * r;
    pn = tmpP[0] * nx + tmpP[1] * ny + tmpP[2] * nz;
    
    pr = pn;
    pf = pxir;
    ptheta2 = pvr * pvr / xi2;
    
    //printf( "pr = %.16e, prT = %.16e\n", pr, prT );
    
    //printf( " a = %.16e, r = %.16e\n", a, r );
    //printf( "D = %.16e, ww = %.16e, rho = %.16e, Lambda = %.16e, xi = %.16e\npr = %.16e, pf = %.16e, deltaR = %.16e, deltaT = %.16e\n",
    //D, ww, sqrt(rho2), Lambda, sqrt(xi2), pr, pf, deltaR, deltaT );
    /* Eqs. 5.36 - 5.46 of BB1 */
    /* Note that the tortoise prT appears only in the quartic term, explained in Eqs. 14 and 15 of Tarrachini et al. */
    Hns =
    GET_SQRT (1. + prT * prT * prT * prT * qq * u2 + ptheta2 / rho2 +
          pf * pf * rho2 / (Lambda * xi2) +
          pr * pr * deltaR / rho2) / GET_SQRT (Lambda / (rho2 * deltaT)) +
    pf * ww / Lambda;
    
    //printf( "term 1 in Hns: %.16e\n",  prT*prT*prT*prT*qq*u2 );
    //printf( "term 2 in Hns: %.16e\n", ptheta2/rho2 );
    //printf( "term 3 in Hns = %.16e\n", pf*pf*rho2/(Lambda*xi2) );
    //printf( "term 4 in Hns = %.16e\n", pr*pr*deltaR/rho2 );
    //printf( "term 5 in Hns = %.16e\n", Lambda/(rho2*deltaT) );
    //printf( "term 6 in Hns = %.16e\n", pf*ww/Lambda );
    /* Eqs. 5.30 - 5.33 of BB1 */
    B = GET_SQRT (deltaT);
    w = ww / Lambda;
    nu = 0.5 * GET_LN (deltaT * rho2 / Lambda);
    MU = 0.5 * GET_LN (rho2);
    /* dLambda/dr */
    Lambda_r = 4. * r * w2 - a2 * deltaT_r * xi2;
    
    ww_r =
    2. * a - (a2 * a * coeffs->b3 * eta) * u2 - coeffs->bb3 * eta * a * u2;
    /* Eqs. 5.47a - 5.47d of BB1 */
    BR =
    (-2. * deltaT + GET_SQRT (deltaR) * deltaT_r) / (2. * GET_SQRT (deltaR * deltaT));
    wr = (-Lambda_r * ww + Lambda * ww_r) / (Lambda * Lambda);
    nur =
    (r / rho2 +
     (w2 * (-4. * r * deltaT + w2 * deltaT_r)) / (2. * deltaT * Lambda));
    mur = (r / rho2 - 1. / sqrt (deltaR));
    /* Eqs. 5.47f - 5.47h of BB1 */
    wcos = -2. * a2 * costheta * deltaT * ww / (Lambda * Lambda);
    nucos = a2 * costheta * w2 * (w2 - deltaT) / (rho2 * Lambda);
    mucos = a2 * costheta / rho2;
    /* Eq. 5.52 of BB1, (YP) simplified */
    //Q = 1. + pvr*pvr/(exp(2.*MU)*xi2) + exp(2.*nu)*pxir*pxir/(B*B*xi2) + pn*pn*deltaR/exp(2.*MU);
    Q =
    1. + pvr * pvr / (rho2 * xi2) +
    deltaT * rho2 / Lambda * pxir * pxir / (B * B * xi2) +
    pn * pn * deltaR / rho2;
    
    pn2 = pr * pr * deltaR / rho2;
    pp = Q - 1.;
    
    //printf( "pn2 = %.16e, pp = %.16e\n", pn2, pp );
    //printf( "sigmaKerr = %.16e, sigmaStar = %.16e\n", sKerr_z, sStar_z );
    /* Eq. 5.68 of BB1, (YP) simplified for aa=bb=0. */
    /*
     deltaSigmaStar_x=(- 8.*aa*(1. + 3.*pn2*r - pp*r)*sKerr_x - 8.*bb*(1. + 3.*pn2*r - pp*r)*sStar_x +
     eta*(-8.*sKerr_x - 36.*pn2*r*sKerr_x + 3.*pp*r*sKerr_x + 14.*sStar_x - 30.*pn2*r*sStar_x + 4.*pp*r*sStar_x))/(12.*r);
     
     deltaSigmaStar_y=(-8.*aa*(1. + 3.*pn2*r - pp*r)*sKerr_y - 8.*bb*(1. + 3.*pn2*r - pp*r)*sStar_y +
     eta*(-8.*sKerr_y - 36.*pn2*r*sKerr_y + 3.*pp*r*sKerr_y + 14.*sStar_y - 30.*pn2*r*sStar_y + 4.*pp*r*sStar_y))/(12.*r);
     
     deltaSigmaStar_z=(-8.*aa*(1. + 3.*pn2*r - pp*r)*sKerr_z - 8.*bb*(1. + 3.*pn2*r - pp*r)*sStar_z +
     eta*(-8.*sKerr_z - 36.*pn2*r*sKerr_z + 3.*pp*r*sKerr_z + 14.*sStar_z - 30.*pn2*r*sStar_z + 4.*pp*r*sStar_z))/(12.*r);
     */
    deltaSigmaStar_x =
    eta * (-8. * sKerr_x - 36. * pn2 * r * sKerr_x + 3. * pp * r * sKerr_x +
           14. * sStar_x - 30. * pn2 * r * sStar_x +
           4. * pp * r * sStar_x) / (12. * r);
    
    deltaSigmaStar_y =
    eta * (-8. * sKerr_y - 36. * pn2 * r * sKerr_y + 3. * pp * r * sKerr_y +
           14. * sStar_y - 30. * pn2 * r * sStar_y +
           4. * pp * r * sStar_y) / (12. * r);
    
    deltaSigmaStar_z =
    eta * (-8. * sKerr_z - 36. * pn2 * r * sKerr_z + 3. * pp * r * sKerr_z +
           14. * sStar_z - 30. * pn2 * r * sStar_z +
           4. * pp * r * sStar_z) / (12. * r);
    
    
    /* Now compute the additional 3.5PN terms. */
    /* The following gauge parameters correspond to those given by
     * Eqs. (69) and (70) of BB2 (aaa -> a0, bbb -> b0).
     * In SEOBNRv1 model, we chose to set all of them to zero,
     * described between Eqs. (3) and (4).
     */
    /*
     aaa = -3./2.*eta;
     bbb = -5./4.*eta;
     a1 = eta*eta/2.;
     a2 = -(1./8.)*eta*(-7. + 8.*eta);
     a3 = -((9.*eta*eta)/16.);
     b1 = 1./16.*eta*(9. + 5.*eta);
     b2 = -(1./8.)*eta*(-17. + 5.*eta);
     b3 = -3./8.*eta*eta;
     */
    /*aaa = 0.;
     bbb = 0.;
     a13P5 = 0.;
     a23P5 = 0.;
     a33P5 = 0.;
     b13P5 = 0.;
     b23P5 = 0.;
     b33P5 = 0.;
     */
    /* Eq. 52 of BB2, (YP) simplified for zero gauge parameters */
    /*
     sMultiplier1 =-(2.*(24.*b23P5 + eta*(-353. + 27.*eta) + bbb*(56. + 60.*eta)) +
     2.*(24.*b13P5 - 24.*b23P5 + bbb*(14. - 66.*eta) + 103.*eta - 60.*eta*eta)*pp*
     r + 120.*(2.*b33P5 - 3.*eta*(bbb + eta))*pn2*pn2*r*r +
     (-48.*b13P5 + 4.*bbb*(1. + 3.*eta) + eta*(23. + 3.*eta))*pp*pp*
     r*r + 6.*pn2*r*(16.*b13P5 + 32.*b23P5 + 24.*b33P5 - 47.*eta +
     54.*eta*eta + 24.*bbb*(1. + eta) +
     (24.*b13P5 - 24.*b33P5 - 16.*eta + 21.*eta*eta + bbb*(-2. + 30.*eta))*pp*
     r))/(72.*r*r);
     */
    sMultiplier1 =
    -(2. * eta * (-353. + 27. * eta) +
      2. * (103. * eta - 60. * eta * eta) * pp * r +
      120. * (-3. * eta * eta) * pn2 * pn2 * r * r +
      (eta * (23. + 3. * eta)) * pp * pp * r * r +
      6. * pn2 * r * (-47. * eta + 54. * eta * eta +
                      (-16. * eta +
                       21. * eta * eta) * pp * r)) / (72. * r * r);
    /* Eq. 52 of BB2, (YP) simplified for zero gauge parameters */
    /*
     sMultiplier2 = (-16.*(6.*a23P5 + 7.*eta*(8. + 3.*eta) + aaa*(14. + 15.*eta)) +
     4.*(-24.*a13P5 + 24.*a23P5 - 109.*eta + 51.*eta*eta + 2.*aaa*(-7. + 33.*eta))*
     pp*r + 30.*(-16.*a33P5 + 3.*eta*(8.*aaa + 9.*eta))*pn2*pn2*r*r +
     (96.*a13P5 - 45.*eta - 8.*aaa*(1. + 3.*eta))*pp*pp*r*r -
     6.*pn2*r*(32.*a13P5 + 64.*a23P5 + 48.*a33P5 + 16.*eta + 147.*eta*eta +
     48.*aaa*(1. + eta) + (48.*a13P5 - 48.*a33P5 - 6.*eta + 39.*eta*eta +
     aaa*(-4. + 60.*eta))*pp*r))/(144.*r*r);
     */
    sMultiplier2 =
    (-16. * (7. * eta * (8. + 3. * eta)) +
     4. * (-109. * eta + 51. * eta * eta) * pp * r +
     810. * eta * eta * pn2 * pn2 * r * r - 45. * eta * pp * pp * r * r -
     6. * pn2 * r * (16. * eta + 147. * eta * eta +
                     (-6. * eta +
                      39. * eta * eta) * pp * r)) / (144. * r * r);
    /* Eq. 52 of BB2 */
    deltaSigmaStar_x +=
    sMultiplier1 * sigmaStar->data[0] + sMultiplier2 * sigmaKerr->data[0];
    deltaSigmaStar_y +=
    sMultiplier1 * sigmaStar->data[1] + sMultiplier2 * sigmaKerr->data[1];
    deltaSigmaStar_z +=
    sMultiplier1 * sigmaStar->data[2] + sMultiplier2 * sigmaKerr->data[2];
    
    /* And now the (calibrated) 4.5PN term */
    deltaSigmaStar_x += coeffs->d1 * eta * sigmaStar->data[0] / (r * r * r);
    deltaSigmaStar_y += coeffs->d1 * eta * sigmaStar->data[1] / (r * r * r);
    deltaSigmaStar_z += coeffs->d1 * eta * sigmaStar->data[2] / (r * r * r);
    deltaSigmaStar_x += coeffs->d1v2 * eta * sigmaKerr->data[0] / (r * r * r);
    deltaSigmaStar_y += coeffs->d1v2 * eta * sigmaKerr->data[1] / (r * r * r);
    deltaSigmaStar_z += coeffs->d1v2 * eta * sigmaKerr->data[2] / (r * r * r);
    
    
    //printf( "deltaSigmaStar_x = %.16e, deltaSigmaStar_y = %.16e, deltaSigmaStar_z = %.16e\n",
    //   deltaSigmaStar_x, deltaSigmaStar_y, deltaSigmaStar_z );
    
    sx = sStar_x + deltaSigmaStar_x;
    sy = sStar_y + deltaSigmaStar_y;
    sz = sStar_z + deltaSigmaStar_z;
    
    
    sxi = sx * xi_x + sy * xi_y + sz * xi_z;
    sv = sx * vx + sy * vy + sz * vz;
    sn = sx * nx + sy * ny + sz * nz;
    
    s3 = sx * e3_x + sy * e3_y + sz * e3_z;
    /* Eq. 3.45 of BB1, second term */
    Hwr =
    (GET_EXP (-3. * MU - nu) * GET_SQRT (deltaR) *
     (GET_EXP (2. * (MU + nu)) * pxir * pxir * sv -
      B * GET_EXP (MU + nu) * pvr * pxir * sxi +
      B * B * xi2 * (GET_EXP (2. * MU) * (GET_SQRT (Q) + Q) * sv +
                     pn * pvr * sn * GET_SQRT (deltaR) -
                     pn * pn * sv * deltaR))) / (2. * B * (1. +
                                                           GET_SQRT (Q)) *
                                                 GET_SQRT (Q) * xi2);
    /* Eq. 3.45 of BB1, third term */
    Hwcos =
    (GET_EXP (-3. * MU - nu) *
     (sn *
      (-(GET_EXP (2. * (MU + nu)) * pxir * pxir) +
       B * B * (pvr * pvr - GET_EXP (2. * MU) * (GET_SQRT (Q) + Q) * xi2)) -
      B * pn * (B * pvr * sv -
                GET_EXP (MU +
                     nu) * pxir * sxi) * GET_SQRT (deltaR))) / (2. * B * (1. +
                                                                      GET_SQRT
                                                                      (Q)) *
                                                            GET_SQRT (Q));
    /* Eq. 3.44 of BB1, leading term */
    HSOL =
    (GET_EXP (-MU + 2. * nu) * (-B + GET_EXP (MU + nu)) * pxir * s3) / (B * B *
                                                                GET_SQRT (Q) *
                                                                xi2);
    /* Eq. 3.44 of BB1, next-to-leading term */
    HSONL =
    (GET_EXP (-2. * MU + nu) *
     (-(B * GET_EXP (MU + nu) * nucos * pxir * (1. + 2. * GET_SQRT (Q)) * sn * xi2) +
      (-(BR * GET_EXP (MU + nu) * pxir * (1. + GET_SQRT (Q)) * sv) +
       B * (GET_EXP (MU + nu) * nur * pxir * (1. + 2. * GET_SQRT (Q)) * sv +
            B * mur * pvr * sxi + B * sxi * (-(mucos * pn * xi2) +
                                             GET_SQRT (Q) * (mur * pvr -
                                                         nur * pvr + (-mucos +
                                                                      nucos) *
                                                         pn * xi2)))) *
      GET_SQRT (deltaR))) / (B * B * (GET_SQRT (Q) + Q) * xi2);
    /* Eq. 3.43 and 3.45 of BB1 */
    Hs = w * s3 + Hwr * wr + Hwcos * wcos + HSOL + HSONL;
    /* Eq. 5.70 of BB1, last term */
    Hss = -0.5 * u3 * (sx * sx + sy * sy + sz * sz - 3. * sn * sn);
    /* Add correction for leading-order spin-induced quadrupole, relevant for BNS - this is zero when kappa_1,2=1 */
    /* Eq. 5.70 of BB1 */
    H = Hns + Hs + Hss;
    
    /* Add the additional calibrated term */
    REAL8 Hess = 0;
    Hess = coeffs->dheffSS * eta * (sKerr_x * sStar_x + sKerr_y * sStar_y +
                             sKerr_z * sStar_z) / (r * r * r * r);
    
    /* One more calibrated term proportional to S1^2+S2^2. Note that we use symmetric expressions of m1,m2 and S1,S2 */
    /*H += coeffs->dheffSSv2 * eta / (r*r*r*r) / (1.-4.*eta)
     * ( (sKerr_x*sKerr_x + sKerr_y*sKerr_y + sKerr_z*sKerr_z)*(1.-4.*eta+2.*eta*eta)
     +(sKerr_x*sStar_x + sKerr_y*sStar_y + sKerr_z*sStar_z)*(-2.*eta+4.*eta*eta)
     +(sStar_x*sStar_x + sStar_y*sStar_y + sStar_z*sStar_z)*(2.*eta*eta) );*/
    Hess += coeffs->dheffSSv2 * eta / (r * r * r * r)
    * (s1Vec->data[0] * s1Vec->data[0] + s1Vec->data[1] * s1Vec->data[1] +
       s1Vec->data[2] * s1Vec->data[2] + s2Vec->data[0] * s2Vec->data[0] +
       s2Vec->data[1] * s2Vec->data[1] + s2Vec->data[2] * s2Vec->data[2]);
    H += Hess;
    //printf( "Hns = %.16e, Hs = %.16e, Hss = %.16e\n", Hns, Hs, Hss );
    //printf( "H = %.16e\n", H );
    /* Real Hamiltonian given by Eq. 2, ignoring the constant -1. */
    Hreal = GET_SQRT (1. + 2. * eta * (H - 1.));
    return Hreal;
}


double GSLSpinHamiltonianWrapper( double x, void *params )
{
    HcapDerivParams *dParams = (HcapDerivParams *)params;
    
    SpinEOBParams *seobParams = dParams->params;
    SpinEOBHCoeffs *seobHamCoeffs = seobParams->seobCoeffs;
    REAL8 tmpVec[12];
    REAL8 s1normData[3], s2normData[3], sKerrData[3], sStarData[3];
    
    /* These are the vectors which will be used in the call to the Hamiltonian */
    REAL8Vector r, p, spin1, spin2, spin1norm, spin2norm;
    REAL8Vector sigmaKerr, sigmaStar;
    
    int i;
    REAL8 a;
    REAL8 m1 = seobParams->m1;
    REAL8 m2 = seobParams->m2;
    REAL8 mT2 = (m1+m2)*(m1+m2);
    
    /* Use a temporary vector to avoid corrupting the main function */
    memcpy( tmpVec, dParams->values, sizeof(tmpVec) );
    
    /* Set the relevant entry in the vector to the correct value */
    tmpVec[dParams->varyParam] = x;
    
    /* Set the LAL-style vectors to point to the appropriate things */
    r.length = p.length = spin1.length = spin2.length = spin1norm.length = spin2norm.length = 3;
    sigmaKerr.length = sigmaStar.length = 3;
    r.data     = tmpVec;
    p.data     = tmpVec+3;
    spin1.data = tmpVec+6;
    spin2.data = tmpVec+9;
    spin1norm.data = s1normData;
    spin2norm.data = s2normData;
    sigmaKerr.data = sKerrData;
    sigmaStar.data = sStarData;
    
    memcpy( s1normData, tmpVec+6, 3*sizeof(REAL8) );
    memcpy( s2normData, tmpVec+9, 3*sizeof(REAL8) );
    
    for ( i = 0; i < 3; i++ )
    {
        s1normData[i] /= mT2;
        s2normData[i] /= mT2;
    }
    
    /* Calculate various spin parameters */
    EOBCalculateSigmaKerr( &sigmaKerr, &spin1norm, &spin2norm );
    EOBCalculateSigmaStar( &sigmaStar, m1, m2, &spin1norm, &spin2norm );
    a = sqrt( sigmaKerr.data[0]*sigmaKerr.data[0] + sigmaKerr.data[1]*sigmaKerr.data[1]
                 + sigmaKerr.data[2]*sigmaKerr.data[2] );
    // EOBCalculateSpinEOBHamCoeffs(seobHamCoeffs, seobParams->eta, a, params);
    return EOBHamiltonian( seobParams->eta, &r, &p, &spin1norm, &spin2norm, &sigmaKerr, &sigmaStar, seobParams->tortoise, seobHamCoeffs ) / seobParams->eta;
}

/* Wrapper for GSL to call the Hamiltonian function */
static double
GSLSpinAlignedHamiltonianWrapper (double x, void *params)
{
    HcapDerivParams *dParams = (HcapDerivParams *) params;

    SpinEOBParams *seobParams = dParams->params;

    REAL8 tmpVec[6];

    /* These are the vectors which will be used in the call to the Hamiltonian */
    REAL8Vector r, p;
    REAL8Vector *s1Vec = seobParams->s1Vec;
    REAL8Vector *s2Vec = seobParams->s2Vec;
    REAL8Vector *sigmaKerr = seobParams->sigmaKerr;
    REAL8Vector *sigmaStar = seobParams->sigmaStar;

    /* Use a temporary vector to avoid corrupting the main function */
    memcpy (tmpVec, dParams->values, sizeof (tmpVec));

    /* Set the relevant entry in the vector to the correct value */
    tmpVec[dParams->varyParam] = x;

    /* Set the LAL-style vectors to point to the appropriate things */
    r.length = p.length = 3;
    r.data = tmpVec;
    p.data = tmpVec + 3;

    REAL8 rX = tmpVec[0];
    REAL8 prT = tmpVec[3];
    REAL8 pphi = tmpVec[4]*rX;

    return EOBHamiltonian(seobParams->eta, &r, &p, s1Vec, s2Vec,
                        sigmaKerr, sigmaStar,
                        seobParams->tortoise,
                        seobParams->seobCoeffs) / seobParams->eta;
    // print_debug("Deriv....(r,prT,pphi)=(%.16e, %.16e, %.16e) H = %.16e\n", 
    //         rX, prT, pphi, H);
    // return H;
}

double
GSLSpinAlignedHamiltonianWrapper_SA (double x, void *params)
{
    HcapDerivParams *dParams = (HcapDerivParams *) params;

    SpinEOBParams *seobParams = dParams->params;

    REAL8 tmpVec[6];

    /* These are the vectors which will be used in the call to the Hamiltonian */
    // REAL8Vector r, p;
    // REAL8Vector *s1Vec = seobParams->s1Vec;
    // REAL8Vector *s2Vec = seobParams->s2Vec;
    // REAL8Vector *sigmaKerr = seobParams->sigmaKerr;
    // REAL8Vector *sigmaStar = seobParams->sigmaStar;

    /* Use a temporary vector to avoid corrupting the main function */
    
    memcpy (tmpVec, dParams->values, sizeof (tmpVec));

    /* Set the relevant entry in the vector to the correct value */
    tmpVec[dParams->varyParam] = x;

    /* Set the LAL-style vectors to point to the appropriate things */
    // r.length = p.length = 3;
    // r.data = tmpVec;
    // p.data = tmpVec + 3;
    REAL8 r = tmpVec[0];
    REAL8 prT = tmpVec[3];
    REAL8 pphi = tmpVec[4]*r;

    // return EOBHamiltonian(seobParams->eta, &r, &p, s1Vec, s2Vec,
    //                     sigmaKerr, sigmaStar,
    //                     seobParams->tortoise,
    //                     seobParams->seobCoeffs) / seobParams->eta;
    REAL8 tmp, H;
    return EOBSAHamiltonian(r, prT, pphi, seobParams->saCoeffs, &tmp) / seobParams->eta;
    // print_debug("Deriv....(r,prT,pphi)=(%.16e, %.16e, %.16e) H = %.16e\n", 
    //         r, prT, pphi, H);
    // return H;
}

/**
 * Calculate the derivative of the Hamiltonian w.r.t. a specific parameter
 * Used by generic spin EOB model, including initial conditions solver.
 */
REAL8 XLALSpinHcapNumDerivWRTParam(
                 const INT paramIdx,      /**<< Index of the parameters */
                 const REAL8 values[],     /**<< Dynamical variables */
                 SpinEOBParams *funcParams /**<< EOB Parameters */
)
{
    static const REAL8 STEP_SIZE = 1.0e-3;

    HcapDerivParams params;

    REAL8 result;

    gsl_function F;
    INT4         gslStatus;

    REAL8 mass1, mass2;

    /* The error in a derivative as measured by GSL */
    REAL8 absErr;
    /* Set up pointers for GSL */
    params.values  = values;
    params.params  = funcParams;

    F.function       = &GSLSpinHamiltonianWrapper;
    F.params         = &params;
    params.varyParam = paramIdx;
// print_debug("here 1\n");

    mass1 = params.params->m1;
    mass2 = params.params->m2;
// print_debug("here 2\n");

    /* Now calculate derivatives w.r.t. the required parameter */
    if ( paramIdx >=6 && paramIdx < 9 )
    {
        gslStatus = gsl_deriv_central( &F, values[paramIdx],
                        STEP_SIZE*mass1*mass1, &result, &absErr );
    }
    else if ( paramIdx >= 9 )
    {
        gslStatus = gsl_deriv_central( &F, values[paramIdx],
                        STEP_SIZE*mass2*mass2, &result, &absErr );
    }
    else
    {
        gslStatus = gsl_deriv_central( &F, values[paramIdx],
                        STEP_SIZE, &result, &absErr );
    }
    if ( gslStatus != GSL_SUCCESS )
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Failure in GSL function");
        return REAL8_FAIL_NAN;
    }
// print_debug("here 3\n");

    //printf( "Abserr = %e\n", absErr );

    return result;
}

/**
 * This function calculates the function \f$\Delta_t(r)\f$ which appears in the spinning EOB
 * potential function. Eqs. 7a and 8.
 */
REAL8
XLALSimIMRSpinEOBHamiltonianDeltaT (SpinEOBHCoeffs * coeffs,
				/**<< Pre-computed coefficients which appear in the function */
				    const REAL8 r,
				/**<< Current orbital radius (in units of total mass) */
				    const REAL8 eta,
				/**<< Symmetric mass ratio */
				    const REAL8 a
				/**<< Normalized deformed Kerr spin */
  )
{
    REAL8 a2;
    REAL8 u, u2, u3, u4, u5;
    REAL8 m1PlusetaKK;

    REAL8 bulk;
    REAL8 logTerms;
    REAL8 deltaU;
    REAL8 deltaT;

    u = 1. / r;
    u2 = u * u;
    u3 = u2 * u;
    u4 = u2 * u2;
    u5 = u4 * u;

    a2 = a * a;

    m1PlusetaKK = -1. + eta * coeffs->KK;

    bulk = 1. / (m1PlusetaKK * m1PlusetaKK) + (2. * u) / m1PlusetaKK + a2 * u2;

    logTerms =
        1. + eta * coeffs->k0 + eta * log (1. + coeffs->k1 * u + coeffs->k2 * u2 +
                        coeffs->k3 * u3 + coeffs->k4 * u4 +
                        coeffs->k5 * u5 +
                        coeffs->k5l * u5 * log (u));
    /*printf(" a = %.16e, u = %.16e\n",a,u);
        printf( "k0 = %.16e, k1 = %.16e, k2 = %.16e, k3 = %.16e , k4 = %.16e, k5 = %.16e, k5l = %.16e\n",coeffs->k0,
        coeffs->k1,coeffs->k2,coeffs->k3,coeffs->k4,coeffs->k5,coeffs->k5l);
        printf( "bulk = %.16e, logTerms = %.16e\n", bulk, logTerms ); */
    deltaU = bulk * logTerms;

        // if ( (coeffs->tidal1->lambda2Tidal != 0. && coeffs->tidal1->omega02Tidal != 0.) || (coeffs->tidal2->lambda2Tidal != 0. && coeffs->tidal2->omega02Tidal != 0.) ) {
        //     deltaU += XLALSimIMRTEOBdeltaUTidal(u, eta, coeffs->tidal1, coeffs->tidal2);
        // }

    deltaT = r * r * deltaU;


    return deltaT;
}


/**
 * This function calculates the function \f$\Delta_r(r)\f$ which appears in the spinning EOB
 * potential function. Eqs. 10a and 10b
 */
REAL8
XLALSimIMRSpinEOBHamiltonianDeltaR (SpinEOBHCoeffs * coeffs,
				/**<< Pre-computed coefficients which appear in the function */
				    const REAL8 r,
				/**<< Current orbital radius (in units of total mass) */
				    const REAL8 eta,
				/**<< Symmetric mass ratio */
				    const REAL8 a
				/**<< Normalized deformed Kerr spin */
  )
{
    REAL8 u2, u3;
    REAL8 D;
    REAL8 deltaT;			/* The potential function, not a time interval... */
    REAL8 deltaR;

    u2 = 1. / (r * r);
    u3 = u2 / r;

    D = 1. + log (1. + 6. * eta * u2 + 2. * (26. - 3. * eta) * eta * u3);

    deltaT = XLALSimIMRSpinEOBHamiltonianDeltaT (coeffs, r, eta, a);

    deltaR = deltaT * D;
    return deltaR;
}

/**
 * Function to calculate the value of omega for the spin-aligned EOB waveform.
 * Can NOT be used in precessing cases. This omega is defined as \f$\dot{y}/r\f$ by setting \f$y=0\f$.
 * The function calculates omega = v/r, by first converting (r,phi,pr,pphi) to Cartesian coordinates
 * in which rVec={r,0,0} and pVec={0,pphi/r,0}, i.e. the effective-test-particle is positioned at x=r,
 * and its velocity along y-axis. Then it computes omega, which is now given by dydt/r = (dH/dp_y)/r.
 */
REAL8
XLALSimIMRSpinAlignedEOBCalcOmega (const REAL8 values[],/**<< Dynamical variables */
				   SpinEOBParams * funcParams,
							/**<< EOB parameters */
                    REAL8 STEP_SIZE /**<< Step size for numerical derivation of H */
  )
{

    HcapDerivParams params;

    /* Cartesian values for calculating the Hamiltonian */
    REAL8 cartValues[6];

    gsl_function F;
    INT4 gslStatus;

    REAL8 omega;
    REAL8 r;

    /* The error in a derivative as measured by GSL */
    REAL8 absErr;

    /* Set up pointers for GSL */
    params.values = cartValues;
    params.params = funcParams;

    F.function = &GSLSpinAlignedHamiltonianWrapper;
    F.params = &params;

    /* Populate the Cartesian values vector */
    /* We can assume phi is zero wlog */
    memset (cartValues, 0, sizeof (cartValues));
    cartValues[0] = r = values[0];
    cartValues[3] = values[2];
    cartValues[4] = values[3] / values[0];

    /* Now calculate omega. In the chosen co-ordinate system, */
    /* we need dH/dpy to calculate this, i.e. varyParam = 4   */
    params.varyParam = 4;
    gslStatus = gsl_deriv_central (&F, cartValues[4], STEP_SIZE, &omega, &absErr);

    if (gslStatus != GSL_SUCCESS)
    {
        PRINT_LOG_INFO (LOG_CRITICAL, "Failure in GSL function");
        return REAL8_FAIL_NAN;
    }

    omega = omega / r;

    return omega;
}


/**
 * Function to calculate the non-Keplerian coefficient for the spin-aligned EOB model.
 * radius \f$r\f$ times the cuberoot of the returned number is \f$r_\Omega\f$ defined in Eq. A2.
 * i.e. the function returns \f$(r_{\Omega} / r)^3\f$.
 */
REAL8
XLALSimIMRSpinAlignedEOBNonKeplerCoeff (const REAL8 values[],
							/**<< Dynamical variables */
					SpinEOBParams * funcParams
							/**<< EOB parameters */
)
{
    REAL8 STEP_SIZE = 1.0e-4;

    REAL8 omegaCirc;

    REAL8 tmpValues[4];

    REAL8 r3;

    /* We need to find the values of omega assuming pr = 0 */
    memcpy (tmpValues, values, sizeof (tmpValues));
    tmpValues[2] = 0.0;

    omegaCirc = XLALSimIMRSpinAlignedEOBCalcOmega (tmpValues, funcParams, STEP_SIZE);
    if (IS_REAL8_FAIL_NAN (omegaCirc))
    {
        return REAL8_FAIL_NAN;
    }

    r3 = values[0] * values[0] * values[0];

    return 1.0 / (omegaCirc * omegaCirc * r3);
}


/**
 * Function to calculate R.H.S. of the ODEs, given dyanmical variables,
 * their derivatives and EOB parameters. Since SEOBNRv1 spin Hamiltonian
 * was implemented for Cartesean coordinates while dynamical evolution was
 * implemented in polar coordinates, we need to perform a transform.
 * This is done in a particular transform in which
 * x = r, y = z = 0, px = pr, py = pphi/r, pz = 0, and
 * omega = v/r = (dy/dt)/r = (dH/dpy)/r, dr/dt = dx/dt = dH/dpx, etc.
 */
int XLALSpinAlignedHcapDerivative(
                  double t,          /**< UNUSED */
                  const REAL8   values[],   /**< dynamical varables */
                  REAL8         dvalues[],  /**< time derivative of dynamical variables */
                  void         *funcParams  /**< EOB parameters */
                  )
{
// clock_t t0, te, t_p1, t_p2, t_p3;
// t0 = clock();

    static const REAL8 STEP_SIZE = 1.0e-4;

    static const INT lMax = 8;

    HcapDerivParams params;

    /* Since we take numerical derivatives wrt dynamical variables */
    /* but we want them wrt time, we use this temporary vector in  */
    /* the conversion */
    REAL8           tmpDValues[6];

    /* Cartesian values for calculating the Hamiltonian */
    REAL8           cartValues[6];

    REAL8           H; //Hamiltonian
    REAL8           flux;

    gsl_function F;
    INT4         gslStatus;
    // UINT SpinAlignedEOBversion;
    UINT i;

    REAL8Vector rVec, pVec;
    REAL8 rData[3], pData[3];

    /* We need r, phi, pr, pPhi to calculate the flux */
    REAL8       r;
    REAL8Vector polarDynamics;
    REAL8Vector cartDynamics;
    REAL8       polData[4];

    REAL8 mass1, mass2, eta;

    /* Spins */
    REAL8Vector *s1Vec = NULL;
    REAL8Vector *s2Vec = NULL;
    REAL8Vector *sKerr = NULL;
    REAL8Vector *sStar = NULL;

    REAL8 a;

    REAL8 omega;

    /* EOB potential functions */
    REAL8 DeltaT, DeltaR;
    REAL8 csi;

    /* The error in a derivative as measured by GSL */
    REAL8 absErr;

    /* Declare NQC coefficients */
    EOBNonQCCoeffs *nqcCoeffs = NULL;

    /* Set up pointers for GSL */ 
    params.values  = cartValues;
    params.params  = (SpinEOBParams *)funcParams;
    nqcCoeffs = params.params->nqcCoeffs;

    s1Vec = params.params->s1Vec;
    s2Vec = params.params->s2Vec;
    sKerr = params.params->sigmaKerr;
    sStar = params.params->sigmaStar;

    F.function = &GSLSpinAlignedHamiltonianWrapper;
    F.params   = &params;

    mass1 = params.params->m1;
    mass2 = params.params->m2;
    eta   = params.params->eta;

    // SpinAlignedEOBversion = 4;

    r = values[0];

    /* Since this is spin aligned, I make the assumption */
    /* that the spin vector is along the z-axis.         */
    a  = sKerr->data[2];

    /* Calculate the potential functions and the tortoise coordinate factor csi,
        given by Eq. 28 of Pan et al. PRD 81, 084041 (2010) */
    DeltaT = XLALSimIMRSpinEOBHamiltonianDeltaT( params.params->seobCoeffs, r, eta, a );
    DeltaR = XLALSimIMRSpinEOBHamiltonianDeltaR( params.params->seobCoeffs, r, eta, a );
    csi    = sqrt( DeltaT * DeltaR ) / (r*r + a*a);
    //printf("DeltaT = %.16e, DeltaR = %.16e, a = %.16e\n",DeltaT,DeltaR,a);
    // print_debug( "csi in derivatives function = %.16e\n", 1./csi );

    /* Populate the Cartesian values vector, using polar coordinate values */
    /* We can assume phi is zero wlog */
    memset( cartValues, 0, sizeof( cartValues ) );
    cartValues[0] = values[0];
    cartValues[3] = values[2];
    cartValues[4] = values[3] / values[0];

    /* Now calculate derivatives w.r.t. each Cartesian variable */
    for ( i = 0; i < 6; i++ )
    {
        params.varyParam = i;
        gslStatus = gsl_deriv_central( &F, cartValues[i], 
                        STEP_SIZE, &tmpDValues[i], &absErr );

        if ( gslStatus != GSL_SUCCESS )
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "Failure in GSL function" );
            return CEV_FAILURE;
        }
    }

    /* Calculate the Cartesian vectors rVec and pVec */
    polarDynamics.length = 4;
    polarDynamics.data   = polData;
    cartDynamics.length = 6;
    cartDynamics.data = cartValues;

    memcpy( polData, values, sizeof( polData ) );
    rVec.length = pVec.length = 3;
    rVec.data   = rData;
    pVec.data   = pData;

    memset( rData, 0, sizeof(rData) );
    memset( pData, 0, sizeof(pData) );

    rData[0] = values[0];
    pData[0] = values[2];
    pData[1] = values[3] / values[0];
// t_p1 = clock();
    /* Calculate Hamiltonian using Cartesian vectors rVec and pVec */
    H =  EOBHamiltonian( eta, &rVec, &pVec, s1Vec, s2Vec, sKerr, sStar, params.params->tortoise, params.params->seobCoeffs );
// t_p2 = clock();
    params.params->cache[0] = H;
    // print_debug( "invcsi = %.16e, csi = %.16e, ham = %.16e ( tortoise = %d)\n", 1./csi, csi, H, params.params->tortoise );
    //exit(1);
    //if ( values[0] > 1.3 && values[0] < 3.9 ) printf( "r = %e\n", values[0] );
    //if ( values[0] > 1.3 && values[0] < 3.9 ) printf( "Hamiltonian = %e\n", H );
    H = H * (mass1 + mass2);


    /*if ( values[0] > 1.3 && values[0] < 3.9 ) printf( "Cartesian derivatives:\n%f %f %f %f %f %f\n",
        tmpDValues[3], tmpDValues[4], tmpDValues[5], -tmpDValues[0], -tmpDValues[1], -tmpDValues[2] );*/

    /* Now calculate omega, and hence the flux */
    omega = tmpDValues[4] / r;
    flux  = XLALInspiralSpinFactorizedFlux( &polarDynamics, &cartDynamics, nqcCoeffs, omega, csi * tmpDValues[3], omega*r, params.params, H/(mass1+mass2), lMax );
// t_p3 = clock();
    /* Looking at the non-spinning model, I think we need to divide the flux by eta */
    flux = flux / eta;

    // print_debug( "Flux in derivatives function = %.16e\n", flux );

    /* Now we can calculate the final (spherical) derivatives */
    /* csi is needed because we use the tortoise co-ordinate */
    /* Right hand side of Eqs. 10a - 10d of Pan et al. PRD 84, 124052 (2011) */
    dvalues[0] = csi * tmpDValues[3];
    dvalues[1] = omega;
    /* Note: in this special coordinate setting, namely y = z = 0, dpr/dt = dpx/dt + dy/dt * py/r, where py = pphi/r */ 
    dvalues[2] = - tmpDValues[0] + tmpDValues[4] * values[3] / (r*r);
    // CORRECTIONFR
// print_debug("d0 = %.16e, d1 = %.16e, d2 = %.16e\n", dvalues[0], dvalues[1], dvalues[2]);
// print_debug("values = (%.16e, %.16e, %.16e, %.16e)\n", values[0], values[1], values[2], values[3]);
// print_debug("cartvalues = (%.16e, %.16e, %.16e, %.16e, %.16e, %.16e)\n", 
//         cartValues[0], cartValues[1], cartValues[2], cartValues[3], cartValues[4], cartValues[5]);

#if 1
    REAL8 cFr, cFf, prDot;
    // CalculateEccCorrectionToFlux(eta, params.params->chi1, params.params->chi1, r, dvalues[0], dvalues[2], &cFr, &cFf);
    prDot = dvalues[2] - ( values[2] / values[3] / csi ) * flux / omega;
    CalculateEccCorrectionToFluxV3(eta, params.params->chi1, params.params->chi2, r, dvalues[0], prDot, &cFr, &cFf, params.params->e0);
    dvalues[2] = dvalues[2] * csi - ( values[2] / values[3] ) * flux * cFr / omega;
    dvalues[3] = - flux * cFf / omega;
    // print_debug("r = %.16e, dr = %.16e, prDot = %.16e\n", r, dvalues[0], prDot);
    // print_debug("Old:cFr = %.16e, cFf = %.16e, flux = %.16e\n", cFr, cFf, flux);
#else
    dvalues[2] = dvalues[2] * csi - ( values[2] / values[3] ) * flux / omega;
    dvalues[3] = - flux / omega;
#endif
    // print_debug("cFr = %f, cFf = %f, vr = %f, prDot = %f\n", cFr, cFf, dvalues[0], dvalues[2]);
    //if ( values[0] > 1.3 && values[0] < 3.9 ) printf("Values:\n%f %f %f %f\n", values[0], values[1], values[2], values[3] );

    //if ( values[0] > 1.3 && values[0] < 3.9 ) printf("Derivatives:\n%f %f %f %f\n", dvalues[0], r*dvalues[1], dvalues[2], dvalues[3] );

    if ( isnan( dvalues[0] ) || isnan( dvalues[1] ) || isnan( dvalues[2] ) || isnan( dvalues[3] ) )
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Deriv is nan: %e %e %e %e\n", dvalues[0], dvalues[1], dvalues[2], dvalues[3] );
        return 1;
    }
// te = clock();
// print_debug("Time Cost Old: t_p1 = %.16e, t_p12 = %.16e, t_p23 = %.16e, t_p3e = %.16e\n", 
//         ((REAL8)(t_p1-t0)/CLOCKS_PER_SEC),
//         ((REAL8)(t_p2-t_p1)/CLOCKS_PER_SEC),
//         ((REAL8)(t_p3-t_p2)/CLOCKS_PER_SEC),
//         ((REAL8)(te-t_p3)/CLOCKS_PER_SEC));
    return CEV_SUCCESS;
}

int XLALSpinAlignedHcapDerivative_SA(
                  double t,          /**< UNUSED */
                  const REAL8   values[],   /**< dynamical varables */
                  REAL8         dvalues[],  /**< time derivative of dynamical variables */
                  void         *funcParams  /**< EOB parameters */
                  )
{
// clock_t t0, te, t_p1, t_p2, t_p3;
// t0 = clock();
    static const REAL8 STEP_SIZE = 1.0e-4;

    static const INT lMax = 8;

    HcapDerivParams params;

    /* Since we take numerical derivatives wrt dynamical variables */
    /* but we want them wrt time, we use this temporary vector in  */
    /* the conversion */
    REAL8           tmpDValues[6];

    /* Cartesian values for calculating the Hamiltonian */
    REAL8           cartValues[6];

    REAL8           H; //Hamiltonian
    REAL8           flux;

    gsl_function F;
    INT4         gslStatus;
    // UINT SpinAlignedEOBversion;
    UINT i;

    REAL8Vector rVec, pVec;
    REAL8 rData[3], pData[3];

    /* We need r, phi, pr, pPhi to calculate the flux */
    REAL8       r;
    REAL8Vector polarDynamics;
    REAL8Vector cartDynamics;
    REAL8       polData[4];

    REAL8 mass1, mass2, eta;

    /* Spins */
    REAL8Vector *s1Vec = NULL;
    REAL8Vector *s2Vec = NULL;
    REAL8Vector *sKerr = NULL;
    REAL8Vector *sStar = NULL;

    REAL8 a;

    REAL8 omega;

    /* EOB potential functions */
    REAL8 DeltaT, DeltaR;
    REAL8 csi;

    /* The error in a derivative as measured by GSL */
    REAL8 absErr;

    /* Declare NQC coefficients */
    EOBNonQCCoeffs *nqcCoeffs = NULL;

    /* Set up pointers for GSL */ 
    params.values  = cartValues;
    params.params  = (SpinEOBParams *)funcParams;
    nqcCoeffs = params.params->nqcCoeffs;

    s1Vec = params.params->s1Vec;
    s2Vec = params.params->s2Vec;
    sKerr = params.params->sigmaKerr;
    sStar = params.params->sigmaStar;

    F.function = &GSLSpinAlignedHamiltonianWrapper_SA;
    F.params   = &params;

    mass1 = params.params->m1;
    mass2 = params.params->m2;
    eta   = params.params->eta;

    // SpinAlignedEOBversion = 4;

    r = values[0];

    /* Since this is spin aligned, I make the assumption */
    /* that the spin vector is along the z-axis.         */
    a  = sKerr->data[2];

    /* Calculate the potential functions and the tortoise coordinate factor csi,
        given by Eq. 28 of Pan et al. PRD 81, 084041 (2010) */
    // DeltaT = XLALSimIMRSpinEOBHamiltonianDeltaT( params.params->seobCoeffs, r, eta, a );
    // DeltaR = XLALSimIMRSpinEOBHamiltonianDeltaR( params.params->seobCoeffs, r, eta, a );
    // csi    = sqrt( DeltaT * DeltaR ) / (r*r + a*a);
    //printf("DeltaT = %.16e, DeltaR = %.16e, a = %.16e\n",DeltaT,DeltaR,a);
    // print_debug( "csi in derivatives function = %.16e\n", 1./csi );

    /* Populate the Cartesian values vector, using polar coordinate values */
    /* We can assume phi is zero wlog */
    memset( cartValues, 0, sizeof( cartValues ) );
    cartValues[0] = values[0];
    cartValues[3] = values[2];
    cartValues[4] = values[3] / values[0];

    /* Now calculate derivatives w.r.t. each Cartesian variable */
    for ( i = 0; i < 6; i++ )
    {
        params.varyParam = i;
        gslStatus = gsl_deriv_central( &F, cartValues[i], 
                        STEP_SIZE, &tmpDValues[i], &absErr );

        if ( gslStatus != GSL_SUCCESS )
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "Failure in GSL function" );
            return CEV_FAILURE;
        }
    }

    /* Calculate the Cartesian vectors rVec and pVec */
    polarDynamics.length = 4;
    polarDynamics.data   = polData;
    cartDynamics.length = 6;
    cartDynamics.data = cartValues;

    memcpy( polData, values, sizeof( polData ) );
    rVec.length = pVec.length = 3;
    rVec.data   = rData;
    pVec.data   = pData;

    memset( rData, 0, sizeof(rData) );
    memset( pData, 0, sizeof(pData) );

    rData[0] = values[0];
    pData[0] = values[2];
    pData[1] = values[3] / values[0];
// t_p1 = clock();
    /* Calculate Hamiltonian using Cartesian vectors rVec and pVec */
    // H =  EOBHamiltonian( eta, &rVec, &pVec, s1Vec, s2Vec, sKerr, sStar, params.params->tortoise, params.params->seobCoeffs );
    H = EOBSAHamiltonian(values[0], values[2], values[3], params.params->saCoeffs, &csi);
// t_p2 = clock();
    params.params->cache[0] = H;
    // print_debug( "invcsi = %.16e, csi = %.16e, ham = %.16e ( tortoise = %d)\n", 1./csi, csi, H, params.params->tortoise );
    //exit(1);
    //if ( values[0] > 1.3 && values[0] < 3.9 ) printf( "r = %e\n", values[0] );
    //if ( values[0] > 1.3 && values[0] < 3.9 ) printf( "Hamiltonian = %e\n", H );
    H = H * (mass1 + mass2);


    /*if ( values[0] > 1.3 && values[0] < 3.9 ) printf( "Cartesian derivatives:\n%f %f %f %f %f %f\n",
        tmpDValues[3], tmpDValues[4], tmpDValues[5], -tmpDValues[0], -tmpDValues[1], -tmpDValues[2] );*/

    /* Now calculate omega, and hence the flux */
    omega = tmpDValues[4] / r;
    flux  = XLALInspiralSpinFactorizedFlux( &polarDynamics, &cartDynamics, nqcCoeffs, omega, csi * tmpDValues[3], omega*r, params.params, H/(mass1+mass2), lMax );
// t_p3 = clock();
    /* Looking at the non-spinning model, I think we need to divide the flux by eta */
    flux = flux / eta;

    // print_debug( "Flux in derivatives function = %.16e\n", flux );

    /* Now we can calculate the final (spherical) derivatives */
    /* csi is needed because we use the tortoise co-ordinate */
    /* Right hand side of Eqs. 10a - 10d of Pan et al. PRD 84, 124052 (2011) */
    dvalues[0] = csi * tmpDValues[3];
    dvalues[1] = omega;
    /* Note: in this special coordinate setting, namely y = z = 0, dpr/dt = dpx/dt + dy/dt * py/r, where py = pphi/r */ 
    dvalues[2] = - tmpDValues[0] + tmpDValues[4] * values[3] / (r*r);
    // CORRECTIONFR
// print_debug("d0 = %.16e, d1 = %.16e, d2 = %.16e\n", dvalues[0], dvalues[1], dvalues[2]);
// print_debug("values = (%.16e, %.16e, %.16e, %.16e)\n", values[0], values[1], values[2], values[3]);
// print_debug("cartvalues = (%.16e, %.16e, %.16e, %.16e, %.16e, %.16e)\n", 
//         cartValues[0], cartValues[1], cartValues[2], cartValues[3], cartValues[4], cartValues[5]);
    // CORRECTIONFR
#if 1
    REAL8 cFr, cFf, prDot;
    // CalculateEccCorrectionToFlux(eta, params.params->chi1, params.params->chi1, r, dvalues[0], dvalues[2], &cFr, &cFf);
    prDot = dvalues[2] - ( values[2] / values[3] / csi ) * flux / omega;
    CalculateEccCorrectionToFluxV3(eta, params.params->chi1, params.params->chi2, r, dvalues[0], prDot, &cFr, &cFf, params.params->e0);
    // CalculateEccCorrectionToFluxV3X(r, dvalues[0], prDot, &cFr, &cFf, params.params->e0, params.params->eccCoeffs);
    dvalues[2] = dvalues[2] * csi - ( values[2] / values[3] ) * flux * cFr / omega;
    dvalues[3] = - flux * cFf / omega;
    // print_debug("r = %.16e, dr = %.16e, prDot = %.16e\n", r, dr, prDot);
    // print_debug("SA:cFr = %.16e, cFf = %.16e, flux = %.16e\n", cFr, cFf, flux);
#else
    dvalues[2] = dvalues[2] * csi - ( values[2] / values[3] ) * flux / omega;
    dvalues[3] = - flux / omega;
#endif
    // print_debug("cFr = %f, cFf = %f, vr = %f, prDot = %f\n", cFr, cFf, dvalues[0], dvalues[2]);
    //if ( values[0] > 1.3 && values[0] < 3.9 ) printf("Values:\n%f %f %f %f\n", values[0], values[1], values[2], values[3] );

    //if ( values[0] > 1.3 && values[0] < 3.9 ) printf("Derivatives:\n%f %f %f %f\n", dvalues[0], r*dvalues[1], dvalues[2], dvalues[3] );

    if ( isnan( dvalues[0] ) || isnan( dvalues[1] ) || isnan( dvalues[2] ) || isnan( dvalues[3] ) )
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Deriv is nan: %e %e %e %e\n", dvalues[0], dvalues[1], dvalues[2], dvalues[3] );
        return CEV_FAILURE;
    }
// te = clock();
// print_debug("Time Cost Old: t_p1 = %.16e, t_p12 = %.16e, t_p23 = %.16e, t_p3e = %.16e\n", 
//         ((REAL8)(t_p1-t0)/CLOCKS_PER_SEC),
//         ((REAL8)(t_p2-t_p1)/CLOCKS_PER_SEC),
//         ((REAL8)(t_p3-t_p2)/CLOCKS_PER_SEC),
//         ((REAL8)(te-t_p3)/CLOCKS_PER_SEC));

    return CEV_SUCCESS;
}

int XLALSpinAlignedHcapDerivative_SAConserve(
                  double t,          /**< UNUSED */
                  const REAL8   values[],   /**< dynamical varables */
                  REAL8         dvalues[],  /**< time derivative of dynamical variables */
                  void         *funcParams  /**< EOB parameters */
                  )
{
// clock_t t0, te, t_p1, t_p2, t_p3;
// t0 = clock();
    static const REAL8 STEP_SIZE = 1.0e-4;

    static const INT lMax = 8;

    HcapDerivParams params;

    /* Since we take numerical derivatives wrt dynamical variables */
    /* but we want them wrt time, we use this temporary vector in  */
    /* the conversion */
    REAL8           tmpDValues[6];

    /* Cartesian values for calculating the Hamiltonian */
    REAL8           cartValues[6];

    REAL8           H; //Hamiltonian
    REAL8           flux;

    gsl_function F;
    INT4         gslStatus;
    // UINT SpinAlignedEOBversion;
    UINT i;

    REAL8Vector rVec, pVec;
    REAL8 rData[3], pData[3];

    /* We need r, phi, pr, pPhi to calculate the flux */
    REAL8       r;
    REAL8Vector polarDynamics;
    REAL8Vector cartDynamics;
    REAL8       polData[4];

    REAL8 mass1, mass2, eta;

    /* Spins */
    REAL8Vector *s1Vec = NULL;
    REAL8Vector *s2Vec = NULL;
    REAL8Vector *sKerr = NULL;
    REAL8Vector *sStar = NULL;

    REAL8 a;

    REAL8 omega;

    /* EOB potential functions */
    REAL8 DeltaT, DeltaR;
    REAL8 csi;

    /* The error in a derivative as measured by GSL */
    REAL8 absErr;

    /* Declare NQC coefficients */
    EOBNonQCCoeffs *nqcCoeffs = NULL;

    /* Set up pointers for GSL */ 
    params.values  = cartValues;
    params.params  = (SpinEOBParams *)funcParams;
    nqcCoeffs = params.params->nqcCoeffs;

    s1Vec = params.params->s1Vec;
    s2Vec = params.params->s2Vec;
    sKerr = params.params->sigmaKerr;
    sStar = params.params->sigmaStar;

    F.function = &GSLSpinAlignedHamiltonianWrapper_SA;
    F.params   = &params;

    mass1 = params.params->m1;
    mass2 = params.params->m2;
    eta   = params.params->eta;

    // SpinAlignedEOBversion = 4;

    r = values[0];

    /* Since this is spin aligned, I make the assumption */
    /* that the spin vector is along the z-axis.         */
    a  = sKerr->data[2];

    /* Calculate the potential functions and the tortoise coordinate factor csi,
        given by Eq. 28 of Pan et al. PRD 81, 084041 (2010) */
    // DeltaT = XLALSimIMRSpinEOBHamiltonianDeltaT( params.params->seobCoeffs, r, eta, a );
    // DeltaR = XLALSimIMRSpinEOBHamiltonianDeltaR( params.params->seobCoeffs, r, eta, a );
    // csi    = sqrt( DeltaT * DeltaR ) / (r*r + a*a);
    //printf("DeltaT = %.16e, DeltaR = %.16e, a = %.16e\n",DeltaT,DeltaR,a);
    // print_debug( "csi in derivatives function = %.16e\n", 1./csi );

    /* Populate the Cartesian values vector, using polar coordinate values */
    /* We can assume phi is zero wlog */
    memset( cartValues, 0, sizeof( cartValues ) );
    cartValues[0] = values[0];
    cartValues[3] = values[2];
    cartValues[4] = values[3] / values[0];

    /* Now calculate derivatives w.r.t. each Cartesian variable */
    for ( i = 0; i < 6; i++ )
    {
        params.varyParam = i;
        gslStatus = gsl_deriv_central( &F, cartValues[i], 
                        STEP_SIZE, &tmpDValues[i], &absErr );

        if ( gslStatus != GSL_SUCCESS )
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "Failure in GSL function" );
            return CEV_FAILURE;
        }
    }

    /* Calculate the Cartesian vectors rVec and pVec */
    polarDynamics.length = 4;
    polarDynamics.data   = polData;
    cartDynamics.length = 6;
    cartDynamics.data = cartValues;

    memcpy( polData, values, sizeof( polData ) );
    rVec.length = pVec.length = 3;
    rVec.data   = rData;
    pVec.data   = pData;

    memset( rData, 0, sizeof(rData) );
    memset( pData, 0, sizeof(pData) );

    rData[0] = values[0];
    pData[0] = values[2];
    pData[1] = values[3] / values[0];
// t_p1 = clock();
    /* Calculate Hamiltonian using Cartesian vectors rVec and pVec */
    // H =  EOBHamiltonian( eta, &rVec, &pVec, s1Vec, s2Vec, sKerr, sStar, params.params->tortoise, params.params->seobCoeffs );
    H = EOBSAHamiltonian(values[0], values[2], values[3], params.params->saCoeffs, &csi);
// t_p2 = clock();
    params.params->cache[0] = H;
    // print_debug( "invcsi = %.16e, csi = %.16e, ham = %.16e ( tortoise = %d)\n", 1./csi, csi, H, params.params->tortoise );
    //exit(1);
    //if ( values[0] > 1.3 && values[0] < 3.9 ) printf( "r = %e\n", values[0] );
    //if ( values[0] > 1.3 && values[0] < 3.9 ) printf( "Hamiltonian = %e\n", H );
    // H = H * (mass1 + mass2);


    /*if ( values[0] > 1.3 && values[0] < 3.9 ) printf( "Cartesian derivatives:\n%f %f %f %f %f %f\n",
        tmpDValues[3], tmpDValues[4], tmpDValues[5], -tmpDValues[0], -tmpDValues[1], -tmpDValues[2] );*/

    /* Now calculate omega, and hence the flux */
    omega = tmpDValues[4] / r;
    // flux  = XLALInspiralSpinFactorizedFlux( &polarDynamics, &cartDynamics, nqcCoeffs, omega, csi * tmpDValues[3], omega*r, params.params, H/(mass1+mass2), lMax );
// t_p3 = clock();
    /* Looking at the non-spinning model, I think we need to divide the flux by eta */
    // flux = flux / eta;

    // print_debug( "Flux in derivatives function = %.16e\n", flux );

    /* Now we can calculate the final (spherical) derivatives */
    /* csi is needed because we use the tortoise co-ordinate */
    /* Right hand side of Eqs. 10a - 10d of Pan et al. PRD 84, 124052 (2011) */
    dvalues[0] = csi * tmpDValues[3];
    dvalues[1] = omega;
    /* Note: in this special coordinate setting, namely y = z = 0, dpr/dt = dpx/dt + dy/dt * py/r, where py = pphi/r */ 
    dvalues[2] = - tmpDValues[0] + tmpDValues[4] * values[3] / (r*r);
    // CORRECTIONFR
// print_debug("d0 = %.16e, d1 = %.16e, d2 = %.16e\n", dvalues[0], dvalues[1], dvalues[2]);
// print_debug("values = (%.16e, %.16e, %.16e, %.16e)\n", values[0], values[1], values[2], values[3]);
// print_debug("cartvalues = (%.16e, %.16e, %.16e, %.16e, %.16e, %.16e)\n", 
//         cartValues[0], cartValues[1], cartValues[2], cartValues[3], cartValues[4], cartValues[5]);
    // CORRECTIONFR
#if 0
    REAL8 cFr, cFf, prDot;
    // CalculateEccCorrectionToFlux(eta, params.params->chi1, params.params->chi1, r, dvalues[0], dvalues[2], &cFr, &cFf);
    prDot = dvalues[2] - ( values[2] / values[3] / csi ) * flux / omega;
    CalculateEccCorrectionToFluxV3(eta, params.params->chi1, params.params->chi2, r, dvalues[0], prDot, &cFr, &cFf, params.params->e0);
    // CalculateEccCorrectionToFluxV3X(r, dvalues[0], prDot, &cFr, &cFf, params.params->e0, params.params->eccCoeffs);
    dvalues[2] = dvalues[2] * csi - ( values[2] / values[3] ) * flux * cFr / omega;
    dvalues[3] = - flux * cFf / omega;
    // print_debug("r = %.16e, dr = %.16e, prDot = %.16e\n", r, dr, prDot);
    // print_debug("SA:cFr = %.16e, cFf = %.16e, flux = %.16e\n", cFr, cFf, flux);
#else
    dvalues[2] = dvalues[2] * csi;
    dvalues[3] = 0.0;
#endif
    // print_debug("cFr = %f, cFf = %f, vr = %f, prDot = %f\n", cFr, cFf, dvalues[0], dvalues[2]);
    //if ( values[0] > 1.3 && values[0] < 3.9 ) printf("Values:\n%f %f %f %f\n", values[0], values[1], values[2], values[3] );

    //if ( values[0] > 1.3 && values[0] < 3.9 ) printf("Derivatives:\n%f %f %f %f\n", dvalues[0], r*dvalues[1], dvalues[2], dvalues[3] );

    if ( isnan( dvalues[0] ) || isnan( dvalues[1] ) || isnan( dvalues[2] ) || isnan( dvalues[3] ) )
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Deriv is nan: %e %e %e %e\n", dvalues[0], dvalues[1], dvalues[2], dvalues[3] );
        return CEV_FAILURE;
    }
// te = clock();
// print_debug("Time Cost Old: t_p1 = %.16e, t_p12 = %.16e, t_p23 = %.16e, t_p3e = %.16e\n", 
//         ((REAL8)(t_p1-t0)/CLOCKS_PER_SEC),
//         ((REAL8)(t_p2-t_p1)/CLOCKS_PER_SEC),
//         ((REAL8)(t_p3-t_p2)/CLOCKS_PER_SEC),
//         ((REAL8)(te-t_p3)/CLOCKS_PER_SEC));

    return CEV_SUCCESS;
}


int XLALSpinAlignedHcapDerivative_Conserve(
                  double t,          /**< UNUSED */
                  const REAL8   values[],   /**< dynamical varables */
                  REAL8         dvalues[],  /**< time derivative of dynamical variables */
                  void         *funcParams  /**< EOB parameters */
                  )
{
    static const REAL8 STEP_SIZE = 1.0e-4;

    static const INT lMax = 8;

    HcapDerivParams params;

    /* Since we take numerical derivatives wrt dynamical variables */
    /* but we want them wrt time, we use this temporary vector in  */
    /* the conversion */
    REAL8           tmpDValues[6];

    /* Cartesian values for calculating the Hamiltonian */
    REAL8           cartValues[6];

    REAL8           H; //Hamiltonian
    REAL8           flux;

    gsl_function F;
    INT4         gslStatus;
    // UINT SpinAlignedEOBversion;
    UINT i;

    REAL8Vector rVec, pVec;
    REAL8 rData[3], pData[3];

    /* We need r, phi, pr, pPhi to calculate the flux */
    REAL8       r;
    REAL8Vector polarDynamics;
    REAL8Vector cartDynamics;
    REAL8       polData[4];

    REAL8 mass1, mass2, eta;

    /* Spins */
    REAL8Vector *s1Vec = NULL;
    REAL8Vector *s2Vec = NULL;
    REAL8Vector *sKerr = NULL;
    REAL8Vector *sStar = NULL;

    REAL8 a;

    REAL8 omega;

    /* EOB potential functions */
    REAL8 DeltaT, DeltaR;
    REAL8 csi;

    /* The error in a derivative as measured by GSL */
    REAL8 absErr;

    /* Declare NQC coefficients */
    EOBNonQCCoeffs *nqcCoeffs = NULL;

    /* Set up pointers for GSL */ 
    params.values  = cartValues;
    params.params  = (SpinEOBParams *)funcParams;
    nqcCoeffs = params.params->nqcCoeffs;

    s1Vec = params.params->s1Vec;
    s2Vec = params.params->s2Vec;
    sKerr = params.params->sigmaKerr;
    sStar = params.params->sigmaStar;

    F.function = &GSLSpinAlignedHamiltonianWrapper;
    F.params   = &params;

    mass1 = params.params->m1;
    mass2 = params.params->m2;
    eta   = params.params->eta;

    // SpinAlignedEOBversion = 4;

    r = values[0];

    /* Since this is spin aligned, I make the assumption */
    /* that the spin vector is along the z-axis.         */
    a  = sKerr->data[2];

    /* Calculate the potential functions and the tortoise coordinate factor csi,
        given by Eq. 28 of Pan et al. PRD 81, 084041 (2010) */
    DeltaT = XLALSimIMRSpinEOBHamiltonianDeltaT( params.params->seobCoeffs, r, eta, a );
    DeltaR = XLALSimIMRSpinEOBHamiltonianDeltaR( params.params->seobCoeffs, r, eta, a );
    csi    = sqrt( DeltaT * DeltaR ) / (r*r + a*a);
    //printf("DeltaT = %.16e, DeltaR = %.16e, a = %.16e\n",DeltaT,DeltaR,a);
    //printf( "csi in derivatives function = %.16e\n", csi );

    /* Populate the Cartesian values vector, using polar coordinate values */
    /* We can assume phi is zero wlog */
    memset( cartValues, 0, sizeof( cartValues ) );
    cartValues[0] = values[0];
    cartValues[3] = values[2];
    cartValues[4] = values[3] / values[0];

    /* Now calculate derivatives w.r.t. each Cartesian variable */
    for ( i = 0; i < 6; i++ )
    {
        params.varyParam = i;
        gslStatus = gsl_deriv_central( &F, cartValues[i], 
                        STEP_SIZE, &tmpDValues[i], &absErr );

        if ( gslStatus != GSL_SUCCESS )
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "Failure in GSL function" );
            return CEV_FAILURE;
        }
    }

    /* Calculate the Cartesian vectors rVec and pVec */
    polarDynamics.length = 4;
    polarDynamics.data   = polData;
    cartDynamics.length = 6;
    cartDynamics.data = cartValues;

    memcpy( polData, values, sizeof( polData ) );
    rVec.length = pVec.length = 3;
    rVec.data   = rData;
    pVec.data   = pData;

    memset( rData, 0, sizeof(rData) );
    memset( pData, 0, sizeof(pData) );

    rData[0] = values[0];
    pData[0] = values[2];
    pData[1] = values[3] / values[0];
    /* Calculate Hamiltonian using Cartesian vectors rVec and pVec */
    H =  EOBHamiltonian( eta, &rVec, &pVec, s1Vec, s2Vec, sKerr, sStar, params.params->tortoise, params.params->seobCoeffs );

    //printf( "csi = %.16e, ham = %.16e ( tortoise = %d)\n", csi, H, params.params->tortoise );
    //exit(1);
    //if ( values[0] > 1.3 && values[0] < 3.9 ) printf( "r = %e\n", values[0] );
    //if ( values[0] > 1.3 && values[0] < 3.9 ) printf( "Hamiltonian = %e\n", H );
    H = H * (mass1 + mass2);


    /*if ( values[0] > 1.3 && values[0] < 3.9 ) printf( "Cartesian derivatives:\n%f %f %f %f %f %f\n",
        tmpDValues[3], tmpDValues[4], tmpDValues[5], -tmpDValues[0], -tmpDValues[1], -tmpDValues[2] );*/

    /* Now calculate omega, and hence the flux */
    omega = tmpDValues[4] / r;
    // flux  = XLALInspiralSpinFactorizedFlux( &polarDynamics, &cartDynamics, nqcCoeffs, omega, csi * tmpDValues[3], omega*r, params.params, H/(mass1+mass2), lMax );

    /* Looking at the non-spinning model, I think we need to divide the flux by eta */
    // flux = flux / eta;

    //printf( "Flux in derivatives function = %.16e\n", flux );

    /* Now we can calculate the final (spherical) derivatives */
    /* csi is needed because we use the tortoise co-ordinate */
    /* Right hand side of Eqs. 10a - 10d of Pan et al. PRD 84, 124052 (2011) */
    dvalues[0] = csi * tmpDValues[3];
    dvalues[1] = omega;
    /* Note: in this special coordinate setting, namely y = z = 0, dpr/dt = dpx/dt + dy/dt * py/r, where py = pphi/r */ 
    dvalues[2] = - tmpDValues[0] + tmpDValues[4] * values[3] / (r*r);
    dvalues[2] = dvalues[2] * csi;
    // CORRECTIONFR
#if 0
    REAL8 cFr, cFf, prDot;
    // CalculateEccCorrectionToFlux(eta, params.params->chi1, params.params->chi1, r, dvalues[0], dvalues[2], &cFr, &cFf);
    prDot = dvalues[2] - ( values[2] / values[3] / csi ) * flux / omega;
    CalculateEccCorrectionToFluxV3(eta, params.params->chi1, params.params->chi1, r, dvalues[0], prDot, &cFr, &cFf, params.params->e0);
    dvalues[2] = dvalues[2] * csi - ( values[2] / values[3] ) * flux * cFr / omega;
    dvalues[3] = - flux * cFf / omega;
    // print_debug("%e, %e, %e\n", cFr, cFf, prDot);
#endif
    if ( isnan( dvalues[0] ) || isnan( dvalues[1] ) || isnan( dvalues[2] ) || isnan( dvalues[3] ) )
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Deriv is nan: %e %e %e %e\n", dvalues[0], dvalues[1], dvalues[2], dvalues[3] );
        print_debug("values = (%.16e, %.16e, %.16e, %.16e)\n", values[0], values[1], values[2], values[3]);
        print_debug("DValues = (%.16e, %.16e, %.16e)\n\t(%.16e, %.16e, %.16e)\n", 
            tmpDValues[0], tmpDValues[1], tmpDValues[2], 
            tmpDValues[3], tmpDValues[4], tmpDValues[5]);
        return 1;
    }
    return CEV_SUCCESS;
}

int XLALSpinAlignedHcapDerivative_invConserve(
                  double t,          /**< UNUSED */
                  const REAL8   values[],   /**< dynamical varables */
                  REAL8         dvalues[],  /**< time derivative of dynamical variables */
                  void         *funcParams  /**< EOB parameters */
                  )
{
    static const REAL8 STEP_SIZE = 1.0e-4;

    static const INT lMax = 8;

    HcapDerivParams params;

    /* Since we take numerical derivatives wrt dynamical variables */
    /* but we want them wrt time, we use this temporary vector in  */
    /* the conversion */
    REAL8           tmpDValues[6];

    /* Cartesian values for calculating the Hamiltonian */
    REAL8           cartValues[6];

    REAL8           H; //Hamiltonian
    REAL8           flux;

    gsl_function F;
    INT4         gslStatus;
    // UINT SpinAlignedEOBversion;
    UINT i;

    REAL8Vector rVec, pVec;
    REAL8 rData[3], pData[3];

    /* We need r, phi, pr, pPhi to calculate the flux */
    REAL8       r;
    REAL8Vector polarDynamics;
    REAL8Vector cartDynamics;
    REAL8       polData[4];

    REAL8 mass1, mass2, eta;

    /* Spins */
    REAL8Vector *s1Vec = NULL;
    REAL8Vector *s2Vec = NULL;
    REAL8Vector *sKerr = NULL;
    REAL8Vector *sStar = NULL;

    REAL8 a;

    REAL8 omega;

    /* EOB potential functions */
    REAL8 DeltaT, DeltaR;
    REAL8 csi;

    /* The error in a derivative as measured by GSL */
    REAL8 absErr;

    /* Declare NQC coefficients */
    EOBNonQCCoeffs *nqcCoeffs = NULL;

    /* Set up pointers for GSL */ 
    params.values  = cartValues;
    params.params  = (SpinEOBParams *)funcParams;
    nqcCoeffs = params.params->nqcCoeffs;

    s1Vec = params.params->s1Vec;
    s2Vec = params.params->s2Vec;
    sKerr = params.params->sigmaKerr;
    sStar = params.params->sigmaStar;

    F.function = &GSLSpinAlignedHamiltonianWrapper;
    F.params   = &params;

    mass1 = params.params->m1;
    mass2 = params.params->m2;
    eta   = params.params->eta;

    // SpinAlignedEOBversion = 4;

    r = values[0];

    /* Since this is spin aligned, I make the assumption */
    /* that the spin vector is along the z-axis.         */
    a  = sKerr->data[2];

    /* Calculate the potential functions and the tortoise coordinate factor csi,
        given by Eq. 28 of Pan et al. PRD 81, 084041 (2010) */
    DeltaT = XLALSimIMRSpinEOBHamiltonianDeltaT( params.params->seobCoeffs, r, eta, a );
    DeltaR = XLALSimIMRSpinEOBHamiltonianDeltaR( params.params->seobCoeffs, r, eta, a );
    csi    = sqrt( DeltaT * DeltaR ) / (r*r + a*a);
    //printf("DeltaT = %.16e, DeltaR = %.16e, a = %.16e\n",DeltaT,DeltaR,a);
    //printf( "csi in derivatives function = %.16e\n", csi );

    /* Populate the Cartesian values vector, using polar coordinate values */
    /* We can assume phi is zero wlog */
    memset( cartValues, 0, sizeof( cartValues ) );
    cartValues[0] = values[0];
    cartValues[3] = values[2];
    cartValues[4] = values[3] / values[0];

    /* Now calculate derivatives w.r.t. each Cartesian variable */
    for ( i = 0; i < 6; i++ )
    {
        params.varyParam = i;
        gslStatus = gsl_deriv_central( &F, cartValues[i], 
                        STEP_SIZE, &tmpDValues[i], &absErr );

        if ( gslStatus != GSL_SUCCESS )
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "Failure in GSL function" );
            return CEV_FAILURE;
        }
    }

    /* Calculate the Cartesian vectors rVec and pVec */
    polarDynamics.length = 4;
    polarDynamics.data   = polData;
    cartDynamics.length = 6;
    cartDynamics.data = cartValues;

    memcpy( polData, values, sizeof( polData ) );
    rVec.length = pVec.length = 3;
    rVec.data   = rData;
    pVec.data   = pData;

    memset( rData, 0, sizeof(rData) );
    memset( pData, 0, sizeof(pData) );

    rData[0] = values[0];
    pData[0] = values[2];
    pData[1] = values[3] / values[0];
    /* Calculate Hamiltonian using Cartesian vectors rVec and pVec */
    H =  EOBHamiltonian( eta, &rVec, &pVec, s1Vec, s2Vec, sKerr, sStar, params.params->tortoise, params.params->seobCoeffs );

    //printf( "csi = %.16e, ham = %.16e ( tortoise = %d)\n", csi, H, params.params->tortoise );
    //exit(1);
    //if ( values[0] > 1.3 && values[0] < 3.9 ) printf( "r = %e\n", values[0] );
    //if ( values[0] > 1.3 && values[0] < 3.9 ) printf( "Hamiltonian = %e\n", H );
    H = H * (mass1 + mass2);


    /*if ( values[0] > 1.3 && values[0] < 3.9 ) printf( "Cartesian derivatives:\n%f %f %f %f %f %f\n",
        tmpDValues[3], tmpDValues[4], tmpDValues[5], -tmpDValues[0], -tmpDValues[1], -tmpDValues[2] );*/

    /* Now calculate omega, and hence the flux */
    omega = tmpDValues[4] / r;
    // flux  = XLALInspiralSpinFactorizedFlux( &polarDynamics, &cartDynamics, nqcCoeffs, omega, csi * tmpDValues[3], omega*r, params.params, H/(mass1+mass2), lMax );

    /* Looking at the non-spinning model, I think we need to divide the flux by eta */
    // flux = flux / eta;

    //printf( "Flux in derivatives function = %.16e\n", flux );

    /* Now we can calculate the final (spherical) derivatives */
    /* csi is needed because we use the tortoise co-ordinate */
    /* Right hand side of Eqs. 10a - 10d of Pan et al. PRD 84, 124052 (2011) */
    dvalues[0] = -csi * tmpDValues[3];
    dvalues[1] = -omega;
    /* Note: in this special coordinate setting, namely y = z = 0, dpr/dt = dpx/dt + dy/dt * py/r, where py = pphi/r */ 
    dvalues[2] = tmpDValues[0] - tmpDValues[4] * values[3] / (r*r);
    // CORRECTIONFR
#if 0
    REAL8 cFr, cFf, prDot;
    // CalculateEccCorrectionToFlux(eta, params.params->chi1, params.params->chi1, r, dvalues[0], dvalues[2], &cFr, &cFf);
    prDot = dvalues[2] - ( values[2] / values[3] / csi ) * flux / omega;
    CalculateEccCorrectionToFluxV3(eta, params.params->chi1, params.params->chi1, r, dvalues[0], prDot, &cFr, &cFf, params.params->e0);
    dvalues[2] = dvalues[2] * csi - ( values[2] / values[3] ) * flux * cFr / omega;
    dvalues[3] = - flux * cFf / omega;
    // print_debug("%e, %e, %e\n", cFr, cFf, prDot);
#endif
    if ( isnan( dvalues[0] ) || isnan( dvalues[1] ) || isnan( dvalues[2] ) || isnan( dvalues[3] ) )
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Deriv is nan: %e %e %e %e\n", dvalues[0], dvalues[1], dvalues[2], dvalues[3] );
        print_debug("values = (%.16e, %.16e, %.16e, %.16e)\n", values[0], values[1], values[2], values[3]);
        print_debug("DValues = (%.16e, %.16e, %.16e)\n\t(%.16e, %.16e, %.16e)\n", 
            tmpDValues[0], tmpDValues[1], tmpDValues[2], 
            tmpDValues[3], tmpDValues[4], tmpDValues[5]);
        return 1;
    }
    return CEV_SUCCESS;
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
#define KRONECKER_DELTA(i,j) ((REAL8)((i)==(j)))

REAL8 XLALSimIMRSpinPrecEOBHamiltonian(
               const REAL8    eta,                  /**<< Symmetric mass ratio */
               REAL8Vector    * x,         /**<< Position vector */
               REAL8Vector    * p,	    /**<< Momentum vector (tortoise radial component pr*) */
               REAL8Vector    * s1Vec,     /**<< Spin vector 1 */
               REAL8Vector    * s2Vec,     /**<< Spin vector 2 */
               REAL8Vector    * sigmaKerr, /**<< Spin vector sigma_kerr */
               REAL8Vector    * sigmaStar, /**<< Spin vector sigma_star */
               INT4                      tortoise,  /**<< flag to state whether the momentum is the tortoise co-ord */
	            SpinEOBHCoeffs *coeffs,               /**<< Structure containing various coefficients */
                HyperParams *hParams
               )
{
    /* Flag for debug output */
    // int debugPK = 0;

    /* Dump out inputs when debug flag is set */
    // if(IS_DEBUG){
    //     PRINT_LOG_INFO(LOG_DEBUG, "In Hamiltonian: tortoise flag = %d\n", (int) tortoise );
    //     PRINT_LOG_INFO(LOG_DEBUG, "x = %.16e\t%.16e\t%.16e\n", x->data[0], x->data[1], x->data[2] );
    //     PRINT_LOG_INFO(LOG_DEBUG, "p = %.16e\t%.16e\t%.16e\n", p->data[0], p->data[1], p->data[2] );
    //     PRINT_LOG_INFO(LOG_DEBUG, "sStar = %.16e\t%.16e\t%.16e\n", sigmaStar->data[0],
    //         sigmaStar->data[1], sigmaStar->data[2] );
    //     PRINT_LOG_INFO(LOG_DEBUG, "sKerr = %.16e\t%.16e\t%.16e\n", sigmaKerr->data[0],
    //         sigmaKerr->data[1], sigmaKerr->data[2] );
    // }

    /* Update the Hamiltonian coefficients, if spins are evolving. Right
        now, this code path is always executed. In the future, v3 and v2
        code may be merged, and we want to skip this step in the
        non-precessing limit. */
    int UsePrecH = 1, k;
    UINT jj;
    SpinEOBHCoeffs tmpCoeffs;
    REAL8 L[3] = {0, 0, 0};
    REAL8 Lhat[3] = {0, 0, 0.0};
    cross_product3d(x->data, p->data, L); // Note that L = r x p is invariant under tortoise transform
    REAL8 L_mag = sqrt(inner_product3d(L, L));
    for ( jj = 0; jj < 3; jj++)
    {
        Lhat[jj] = L[jj] / L_mag;
    }

    REAL8 tempS1_p = inner_product3d(s1Vec->data, Lhat);
    REAL8 tempS2_p = inner_product3d(s2Vec->data, Lhat);
    REAL8 S1_perp[3] = {0, 0, 0};
    REAL8 S2_perp[3] = {0, 0, 0};
    REAL8 S_perp[3] = {0,0,0};
    for (jj = 0; jj < 3; jj++)
    {
        S1_perp[jj] = s1Vec->data[jj] - tempS1_p * Lhat[jj];
        S2_perp[jj] = s2Vec->data[jj] - tempS2_p * Lhat[jj];
        S_perp[jj] = S1_perp[jj]+S2_perp[jj];
    }
    REAL8 sKerr_norm = sqrt(inner_product3d(sigmaKerr->data, sigmaKerr->data));
    REAL8 S_con = 0.0;
    if (sKerr_norm>1e-6)
    {
        S_con = sigmaKerr->data[0] * Lhat[0] + sigmaKerr->data[1] * Lhat[1] + sigmaKerr->data[2] * Lhat[2];
        S_con /= (1 - 2 * eta);
        // Last division by 2 is to ensure the spin oebys the Kerr bound.
        S_con += inner_product3d(S_perp, sigmaKerr->data) / sKerr_norm / (1 - 2 * eta) / 2.;
    }

    REAL8 chi = S_con;
    if ( UsePrecH && coeffs->updateHCoeffs )
    {

        REAL8 tmpa; // = magnitude of S_1 + S_2
        tmpa = sqrt(sigmaKerr->data[0]*sigmaKerr->data[0]
                    + sigmaKerr->data[1]*sigmaKerr->data[1]
                    + sigmaKerr->data[2]*sigmaKerr->data[2]);

        // Update coefficients, checking for errors
        // if (coeffs->SpinAlignedEOBversion ==4)
        // {
        if ( XLALSimIMRCalculateSpinPrecEOBHCoeffs_v2( &tmpCoeffs, eta,
            tmpa, chi, 4 ,hParams) == CEV_FAILURE )
        {
            return REAL8_FAIL_NAN;
        }
        // }
        // else 
        // {
        //     if ( XLALSimIMRCalculateSpinPrecEOBHCoeffs( &tmpCoeffs, eta,
        //         tmpa, coeffs->SpinAlignedEOBversion ) == CEV_FAILURE )
        //     {
        //         return REAL8_FAIL_NAN;
        //     }
        // }


        // Copy over underlying model version number
        // tmpCoeffs.SpinAlignedEOBversion = coeffs->SpinAlignedEOBversion;
        tmpCoeffs.updateHCoeffs = coeffs->updateHCoeffs;

        coeffs = &tmpCoeffs;
    }

    REAL8 r, r2, nx, ny, nz;
    REAL8 sKerr_x, sKerr_y, sKerr_z, a, a2;
    REAL8 sStar_x, sStar_y, sStar_z;
    REAL8 e3_x, e3_y, e3_z;
    REAL8 costheta; /* Cosine of angle between Skerr and r */
    REAL8 xi2, xi_x, xi_y, xi_z; /* Cross product of unit vectors in direction of Skerr and r */
    REAL8 vx, vy, vz, pxir, pvr, pn, prT, pr, pf, ptheta2; /*prT is the tortoise pr */
    REAL8 w2, rho2;
    REAL8 u, u2, u3, u4, u5;
    REAL8 bulk, deltaT, deltaR, Lambda;
    REAL8 D, qq, ww, B, w, BR, wr, nur, mur;
    REAL8 wcos, nucos, mucos, ww_r, Lambda_r;
    REAL8 logTerms, deltaU, deltaU_u, Q, deltaT_r, pn2, pp;
    REAL8 deltaSigmaStar_x, deltaSigmaStar_y, deltaSigmaStar_z;
    REAL8 sx, sy, sz, sxi, sv, sn, s3;
    REAL8 H, Hns, Hs, Hss, Hreal, Hwcos, Hwr, HSOL, HSONL;

    /* Terms which come into the 3.5PN mapping of the spins */
    REAL8 sMultiplier1, sMultiplier2;

    /*Temporary p vector which we will make non-tortoise */
    REAL8 tmpP[3] = {0.};

    REAL8 csi;
    REAL8 logu;

    r2 = x->data[0]*x->data[0] + x->data[1]*x->data[1] + x->data[2]*x->data[2];
    r  = sqrt(r2);
    u  = 1./r;
    u2 = u*u;
    u3 = u2*u;
    u4 = u2*u2;
    u5 = u4*u;

    nx = x->data[0] *u;
    ny = x->data[1] *u;
    nz = x->data[2] *u;

    sKerr_x = sigmaKerr->data[0];
    sKerr_y = sigmaKerr->data[1];
    sKerr_z = sigmaKerr->data[2];

    sStar_x = sigmaStar->data[0];
    sStar_y = sigmaStar->data[1];
    sStar_z = sigmaStar->data[2];

    a2 = sKerr_x*sKerr_x + sKerr_y*sKerr_y + sKerr_z*sKerr_z;
    a  = sqrt( a2 );

    if(a !=0.)
    {
        const REAL8 inva = 1./a;
        e3_x = sKerr_x * inva;
        e3_y = sKerr_y * inva;
        e3_z = sKerr_z * inva;
    }
    else
    {
        e3_x = 1./sqrt(3.);
        e3_y = 1./sqrt(3.);
        e3_z = 1./sqrt(3.);
    }
    REAL8 result[3] = {0.0, 0.0, 0.0};
    REAL8 e3[3] = {e3_x, e3_y, e3_z};
    REAL8 nhat[3] = {nx,ny,nz};
    REAL8 lambda_hat[3]={0.0,0.0,0.0};
    cross_product3d(Lhat,nhat,lambda_hat);
    REAL8 nrm = sqrt(inner_product3d(lambda_hat,lambda_hat));
    for (k=0;k<3;k++)
    {
        lambda_hat[k]/=nrm;
    }
    // Check if e_3 is aligned with n

    if (1. - fabs(e3_x * nx + e3_y * ny + e3_z * nz) <= 1.e-8)
    {
        // if (coeffs->SpinAlignedEOBversion == 4)
        // {
        REAL8 angle = 1.8e-3; // This is ~0.1 degrees
        rotate_vector(e3, lambda_hat, angle, result);
        e3_x = result[0];
        e3_y = result[1];
        e3_z = result[2];
        // }
        // else
        // {
        //     e3_x = e3_x+0.1;
        //     e3_y = e3_y+0.1;
        //     const REAL8 invnorm = 1./sqrt(e3_x*e3_x + e3_y*e3_y + e3_z*e3_z);
        //     e3_x = e3_x*invnorm;
        //     e3_y = e3_y*invnorm;
        //     e3_z = e3_z*invnorm;
        // }
    }
    costheta = e3_x*nx + e3_y*ny + e3_z*nz;

    xi2=1. - costheta*costheta;

    xi_x = -e3_z*ny + e3_y*nz;
    xi_y =  e3_z*nx - e3_x*nz;
    xi_z = -e3_y*nx + e3_x*ny;

    vx = -nz*xi_y + ny*xi_z;
    vy =  nz*xi_x - nx*xi_z;
    vz = -ny*xi_x + nx*xi_y;

    w2 = r2 + a2;
    rho2 = r2 + a2*costheta*costheta;

    // if(debugPK)XLAL_PRINT_INFO( "KK = %.16e\n", coeffs->KK );
    const REAL8 invm1PlusetaKK = 1./(-1. + eta * coeffs->KK);
    /* Eq. 5.75 of BB1 */
    bulk = invm1PlusetaKK*(invm1PlusetaKK + (2.*u)) + a2*u2;
    /* Eq. 5.73 of BB1 */
    // use ln(u) = log_2(u)/log_2(e) and the fact that log2 is faster than ln
    // this relies on the compiler evaluating the expression at compile time.
    // which apparently not all do so in stead of 1./log2(exp(1.)) I use the
    // result returned by Maple.
    const REAL8 invlog_2e = 0.69314718055994530941723212145817656807550013436026;
    logu = log2(u)*invlog_2e;
    const REAL8 logarg = coeffs->k1*u + coeffs->k2*u2 + coeffs->k3*u3 + coeffs->k4*u4
                                                + coeffs->k5*u5 + coeffs->k5l*u5*logu;
    logTerms = 1. + eta*coeffs->k0 + eta*log1p(fabs(1. + logarg) - 1.);
    //if(debugPK)XLAL_PRINT_INFO( "bulk = %.16e, logTerms = %.16e\n", bulk, logTerms );
    /* Eq. 5.73 of BB1 */
    deltaU = fabs(bulk*logTerms);
    /* Eq. 5.71 of BB1 */
    deltaT = r2*deltaU;
    /* ddeltaU/du */
    deltaU_u = 2.*(invm1PlusetaKK + a2*u)*logTerms +
        bulk * (eta*(coeffs->k1 + u*(2.*coeffs->k2 + u*(3.*coeffs->k3 + u*(4.*coeffs->k4 + 5.*(coeffs->k5+coeffs->k5l*logu)*u)))))
            / (1. + logarg);
    /* ddeltaT/dr */
    deltaT_r = 2.*r*deltaU - deltaU_u;
    /* Eq. 5.39 of BB1 */
    Lambda = fabs(w2*w2 - a2*deltaT*xi2);
    // RH: this is horrible, but faster than 3 divisions
    const REAL8 invrho2xi2Lambda = 1./(rho2*xi2*Lambda);
    const REAL8 invrho2 = xi2 * (Lambda*invrho2xi2Lambda);
    const REAL8 invxi2 = rho2 * (Lambda*invrho2xi2Lambda);
    const REAL8 invLambda = xi2*rho2*invrho2xi2Lambda;
    /* Eq. 5.83 of BB1, inverse */
    D = 1. + log1p(6.*eta*u2 + 2.*(26. - 3.*eta)*eta*u3);
    /* Eq. 5.38 of BB1 */
    deltaR = deltaT*D;
    /* See Hns below, Eq. 4.34 of Damour et al. PRD 62, 084011 (2000) */
    qq = 2.*eta*(4. - 3.*eta);
    /* See Hns below. In Sec. II D of BB2 b3 and bb3 coeffs are chosen to be zero. */
    ww=2.*a*r + coeffs->b3*eta*a2*a*u + coeffs->bb3*eta*a*u;

    /* We need to transform the momentum to get the tortoise co-ord */
    /* Eq. 28 of Pan et al. PRD 81, 084041 (2010) */
    // RH: this assumes that tortoise can be 0 or 1 or 2.
    csi = sqrt( fabs(deltaT * deltaR) )/ w2;
    // non-unity only for tortoise==1
    const REAL8 csi1 = 1.0 + (1.-fabs(1.-tortoise)) * (csi - 1.0);
    // non-unity only for tortoise==2
    const REAL8 csi2 = 1.0 + (0.5-copysign(0.5, 1.5-tortoise)) * (csi - 1.0);

    // if(debugPK){
    //     XLAL_PRINT_INFO( "csi1(miami) = %.16e\n", csi1 );
    //     XLAL_PRINT_INFO( "csi2(miami) = %.16e\n", csi2 );}

    prT = (p->data[0]*nx + p->data[1]*ny + p->data[2]*nz)*csi2;
    /* p->data is BL momentum vector; tmpP is tortoise momentum vector */
    tmpP[0] = p->data[0] - nx * prT * (1. - 1./csi1);
    tmpP[1] = p->data[1] - ny * prT * (1. - 1./csi1);
    tmpP[2] = p->data[2] - nz * prT * (1. - 1./csi1);

    pxir = (tmpP[0]*xi_x + tmpP[1]*xi_y + tmpP[2]*xi_z) * r;
    pvr  = (tmpP[0]*vx + tmpP[1]*vy + tmpP[2]*vz) * r;
    pn   = tmpP[0]*nx + tmpP[1]*ny + tmpP[2]*nz;

    pr = pn;
    pf = pxir;
    ptheta2 = pvr * pvr *invxi2;

    // if(debugPK) {XLAL_PRINT_INFO( "pr = %.16e, prT = %.16e\n", pr, prT );
    // XLAL_PRINT_INFO( " a = %.16e, r = %.16e\n", a, r );
    // XLAL_PRINT_INFO( "D = %.16e, ww = %.16e, rho = %.16e, Lambda = %.16e, xi = %.16e\npr = %.16e, pf = %.16e, deltaR = %.16e, deltaT = %.16e\n",
    //     D, ww, sqrt(rho2), Lambda, sqrt(xi2), pr, pf, deltaR, deltaT );}

    /* Eqs. 5.36 - 5.46 of BB1 */
    /* Note that the tortoise prT appears only in the quartic term, explained in Eqs. 14 and 15 of Tarrachini et al. */
    Hns = sqrt((1. + ((prT*prT)*(prT*prT))*qq*u2 + ptheta2*invrho2 + pf*pf*rho2*invLambda*invxi2 + pr*pr*deltaR*invrho2)
                * (rho2*deltaT) * invLambda) + pf*ww*invLambda;
    // if(IS_DEBUG){
    // print_debug( "term 1 in Hns: %.16e\n",  prT*prT*prT*prT*qq*u2 );
    // print_debug( "term 2 in Hns: %.16e\n", ptheta2/rho2 );
    // print_debug( "term 3 in Hns = %.16e\n", pf*pf*rho2/(Lambda*xi2) );
    // print_debug( "term 4 in Hns = %.16e\n", pr*pr*deltaR/rho2 );
    // print_debug( "term 5 in Hns = %.16e\n", Lambda/(rho2*deltaT) );
    // print_debug( "term 6 in Hns = %.16e\n", pf*ww/Lambda );
    // print_debug( "term 6 in Hns = %.16e\n", pf );
    // }

    /* Eqs. 5.30 - 5.33 of BB1 */
    B = sqrt(deltaT);
    // RH: this is horrible but faster than 3 divisions
    const REAL8 sqrtdeltaT = B;
    const REAL8 sqrtdeltaR = sqrt(deltaR);
    const REAL8 invdeltaTsqrtdeltaTsqrtdeltaR = 1./(sqrtdeltaT*deltaT*sqrtdeltaR);
    const REAL8 invdeltaT = sqrtdeltaT*(sqrtdeltaR*invdeltaTsqrtdeltaTsqrtdeltaR);
    const REAL8 invsqrtdeltaT = deltaT*(sqrtdeltaR*invdeltaTsqrtdeltaTsqrtdeltaR);
    const REAL8 invsqrtdeltaR = deltaT*sqrtdeltaT*invdeltaTsqrtdeltaTsqrtdeltaR;
    w = ww*invLambda;
    //nu = 0.5 * log(deltaT*rho2/Lambda);
    //MU = 0.5 * log(rho2);
    const REAL8 expnu = sqrt(deltaT*rho2*invLambda);
    const REAL8 expMU = sqrt(rho2);
    // RH: this is horrible but faster than 2 divisions
    const REAL8 invexpnuexpMU = 1./(expnu*expMU);
    const REAL8 invexpnu = expMU*invexpnuexpMU;
    const REAL8 invexpMU = expnu*invexpnuexpMU;
    /* dLambda/dr */
    Lambda_r = 4.*r*w2 - a2*deltaT_r*xi2;

    ww_r=2.*a - (a2*a*coeffs->b3*eta)*u2 - coeffs->bb3*eta*a*u2;
    /* Eqs. 5.47a - 5.47d of BB1 */
    BR = (-deltaT*invsqrtdeltaR + deltaT_r*0.5)*invsqrtdeltaT;
    wr = (-Lambda_r*ww + Lambda*ww_r)*(invLambda*invLambda);
    nur = (r*invrho2 + (w2 * (-4.*r*deltaT + w2*deltaT_r) ) * 0.5*invdeltaT*invLambda );
    mur = (r*invrho2 - invsqrtdeltaR);
    /* Eqs. 5.47f - 5.47h of BB1 */
    wcos  = -2.*(a2*costheta)*deltaT*ww*(invLambda*invLambda);
    nucos = (a2*costheta)*w2*(w2-deltaT)*(invrho2*invLambda);
    mucos = (a2*costheta)*invrho2;
    /* Eq. 5.52 of BB1, (YP) simplified */
    Q = 1. + pvr*pvr*invrho2*invxi2 + pxir*pxir*rho2*invLambda*invxi2 + pn*pn*deltaR*invrho2;
    // if(debugPK){
    //     XLAL_PRINT_INFO( "Q = %.16e, pvr = %.16e, xi2 = %.16e , deltaT = %.16e, rho2 = %.16e, Lambda = %.16e, pxir = %.16e, B = %.16e\n", Q, pvr, xi2, deltaT, rho2, Lambda, pxir, B );
    // }
    pn2 = pr * pr * deltaR * invrho2;
    pp  = Q - 1.;

    // if(debugPK){
    //     XLAL_PRINT_INFO( "pn2 = %.16e, pp = %.16e\n", pn2, pp );
    //     XLAL_PRINT_INFO( "sigmaKerr = %.16e, sigmaStar = %.16e\n", sKerr_z, sStar_z );}

    /* Eq. 5.68 of BB1, (YP) simplified for aa=bb=0. */
    deltaSigmaStar_x=eta*((-8. - 3.*r*(12.*pn2 - pp))*sKerr_x + (14. + (- 30.*pn2 + 4.*pp)*r)*sStar_x)*(1./12.)*u;

    deltaSigmaStar_y=eta*((-8. - 3.*r*(12.*pn2 - pp))*sKerr_y + (14. + (- 30.*pn2 + 4.*pp)*r)*sStar_y)*(1./12.)*u;

    deltaSigmaStar_z=eta*((-8. - 3.*r*(12.*pn2 - pp))*sKerr_z + (14. + (- 30.*pn2 + 4.*pp)*r)*sStar_z)*(1./12.)*u;


    /* Now compute the additional 3.5PN terms. */
    /* Eq. 52 of BB2, (YP) simplified for zero gauge parameters */
    // RH: below is horner(%, [eta,r])
    // sMultiplier1 = -(2.*eta*(-353. + 27.*eta) + 2.*(103.*eta - 60.*eta*eta)*pp*r
    //             + (120.*(-3.))*(eta*eta)*(pn2*pn2)*(r*r) + (eta*(23. + 3.*eta))*(pp*pp)*(r*r )
    //             + 6.*pn2*r*(- 47.*eta + 54.*(eta*eta) + (- 16.*eta + 21.*(eta*eta))*pp*r))
    //             * (1./72.) * u2;
    sMultiplier1 = (-706.0+(206.0*pp-282.0*pn2+(-96.0*pn2*pp+23.0*pp*pp)*r)*r
                    +(54.0+( -120.0*pp+324.0*pn2+(-360.0*pn2*pn2+126.0*pn2*pp
                                                    +3.0*pp*pp)*r)*r)*eta)*eta*u2
                    *(-1./72.0);
    /* Eq. 52 of BB2, (YP) simplified for zero gauge parameters */
    //RH: below is horner(expand(%), [eta,r])
    // sMultiplier2 = (-16.*(7.*eta*(8. + 3.*eta)) + 4.*(- 109.*eta + 51.*eta*eta)*pp*r
    //             + 810.*(eta*eta)*(pn2*pn2)*(r*r) - 45.*eta*(pp*pp)*(r*r)
    //             - 6.*pn2*r*(16.*eta + 147.*eta*eta + (- 6.*eta + 39.*(eta*eta))*pp*r))
    //             * (1./144.) * u2;
    sMultiplier2 = (-56.0/9.0*u2+(-2.0/3.0*pn2*u2-109.0/36.0*pp*u2
                                    +(pn2*pp*u2/4.0-5.0/16.0*pp*pp*u2)*r)*r
                                +(-7.0/3.0*u2+(-49.0/8.0*pn2*u2+17.0/12.0*pp*u2
                                                +(45.0/8.0* pn2*pn2*u2
                                                -13.0/8.0*pn2*pp*u2)*r)*r)*eta)
                    *eta;

// if(IS_DEBUG)
// {
//     print_debug("deltaSigmaStar0 = (%.16e, %.16e, %.16e), pnBar2 = %.16e\n", deltaSigmaStar_x, deltaSigmaStar_y, deltaSigmaStar_z, pn2);
//     print_debug("sMultiplier1 = %.16e, sMultiplier2 = %.16e\n", sMultiplier1, sMultiplier2);
//     print_debug("coeffs->d1 = %.16e, coeffs->d1v2 = %.16e\n", coeffs->d1, coeffs->d1v2);
// }

    /* Eq. 52 of BB2 */
    deltaSigmaStar_x += sMultiplier1*sigmaStar->data[0] + sMultiplier2*sigmaKerr->data[0];
    deltaSigmaStar_y += sMultiplier1*sigmaStar->data[1] + sMultiplier2*sigmaKerr->data[1];
    deltaSigmaStar_z += sMultiplier1*sigmaStar->data[2] + sMultiplier2*sigmaKerr->data[2];

    /* And now the (calibrated) 4.5PN term */
    deltaSigmaStar_x += coeffs->d1 * eta * sigmaStar->data[0] * u3;
    deltaSigmaStar_y += coeffs->d1 * eta * sigmaStar->data[1] * u3;
    deltaSigmaStar_z += coeffs->d1 * eta * sigmaStar->data[2] * u3;
    deltaSigmaStar_x += coeffs->d1v2 * eta * sigmaKerr->data[0] * u3;
    deltaSigmaStar_y += coeffs->d1v2 * eta * sigmaKerr->data[1] * u3;
    deltaSigmaStar_z += coeffs->d1v2 * eta * sigmaKerr->data[2] * u3;


    sx = sStar_x + deltaSigmaStar_x;
    sy = sStar_y + deltaSigmaStar_y;
    sz = sStar_z + deltaSigmaStar_z;


    sxi = sx*xi_x + sy*xi_y + sz*xi_z;
    sv  = sx*vx + sy*vy + sz*vz;
    sn  = sx*nx + sy*ny + sz*nz;

    s3 = sx*e3_x + sy*e3_y + sz*e3_z;
    /* Eq. 3.45 of BB1, second term */
    const REAL8 sqrtQ = sqrt(Q);
    const REAL8 inv2B1psqrtQsqrtQ = 1./(2.*B*(1. + sqrtQ)*sqrtQ);
    Hwr = ((invexpMU*invexpMU*invexpMU*invexpnu)*sqrtdeltaR*((expMU*expMU)*(expnu*expnu)*(pxir*pxir)*sv - B*(expMU*expnu)*pvr*pxir*sxi +
                                                        B*B*xi2*((expMU*expMU)*(sqrtQ + Q)*sv + pn*pvr*sn*sqrtdeltaR - pn*pn*sv*deltaR)))*inv2B1psqrtQsqrtQ*invxi2;
    /* Eq. 3.45 of BB1, third term */
    Hwcos = ((invexpMU*invexpMU*invexpMU*invexpnu)*(sn*(-((expMU*expMU)*(expnu*expnu)*(pxir*pxir)) + B*B*(pvr*pvr - (expMU*expMU)*(sqrtQ + Q)*xi2)) -
                                                B*pn*(B*pvr*sv - (expMU*expnu)*pxir*sxi)*sqrtdeltaR))*inv2B1psqrtQsqrtQ;
    /* Eq. 3.44 of BB1, leading term */
    HSOL = ((expnu*expnu*invexpMU)*(-B + (expMU*expnu))*pxir*s3)/(deltaT*sqrtQ)*invxi2;
    /* Eq. 3.44 of BB1, next-to-leading term */
    HSONL = ((expnu*(invexpMU*invexpMU))*(-(B*expMU*expnu*nucos*pxir*(1. + 2.*sqrtQ)*sn*xi2) +
            (-(BR*(expMU*expnu)*pxir*(1. + sqrtQ)*sv) + B*((expMU*expnu)*nur*pxir*(1. + 2.*sqrtQ)*sv + B*mur*pvr*sxi +
            B*sxi*(-(mucos*pn*xi2) + sqrtQ*(mur*pvr - nur*pvr + (-mucos + nucos)*pn*xi2))))*sqrtdeltaR))*invxi2/(deltaT*(sqrtQ + Q));
    /* Eq. 3.43 and 3.45 of BB1 */
    Hs = w*s3 + Hwr*wr + Hwcos*wcos + HSOL + HSONL;
    /* Eq. 5.70 of BB1, last term */
    Hss = -0.5*u3 * (sx*sx + sy*sy + sz*sz - 3.*sn*sn);
    /* Eq. 5.70 of BB1 */
    H = Hns + Hs + Hss;

    /* Add the additional calibrated term */
    H += coeffs->dheffSS * eta * (sKerr_x*sStar_x + sKerr_y*sStar_y + sKerr_z*sStar_z) *u4;
    /* One more calibrated term proportional to S1^2+S2^2. Note that we use symmetric expressions of m1,m2 and S1,S2 */
    H += coeffs->dheffSSv2 * eta * u4
                            * (s1Vec->data[0]*s1Vec->data[0] + s1Vec->data[1]*s1Vec->data[1] + s1Vec->data[2]*s1Vec->data[2]
                            +s2Vec->data[0]*s2Vec->data[0] + s2Vec->data[1]*s2Vec->data[1] + s2Vec->data[2]*s2Vec->data[2]);
    // if(debugPK){
    //     XLAL_PRINT_INFO( "Hns = %.16e, Hs = %.16e, Hss = %.16e\n", Hns, Hs, Hss );
    //     XLAL_PRINT_INFO( "H = %.16e\n", H );}
    /* Real Hamiltonian given by Eq. 2, ignoring the constant -1. */
    Hreal = sqrt(1. + 2.*eta *(fabs(H) - 1.));

    // if(debugPK)
    //     XLAL_PRINT_INFO( "Hreal = %.16e\n", Hreal );

    if(isnan(Hreal) || IS_DEBUG) 
    {
// print_debug("tmp1 = %.16e, tmp2 = %.16e, tmp3 = %.16e, tmp4 = %.16e\n",(invexpMU*invexpMU*invexpMU*invexpnu)*sqrtdeltaR*inv2B1psqrtQsqrtQ*invxi2,
//     ((expMU*expMU)*(expnu*expnu)*(pxir*pxir)*sv), (- B*(expMU*expnu)*pvr*pxir*sxi),(B*B*xi2*((expMU*expMU)*(sqrtQ + Q)*sv + pn*pvr*sn*sqrtdeltaR - pn*pn*sv*deltaR)));
        PRINT_LOG_INFO(
        LOG_DEBUG,"Inside Hamiltonian: Hreal is a NAN. Printing its components below:");        

        PRINT_LOG_INFO(LOG_DEBUG, "In Hamiltonian: tortoise flag = %d", (int) tortoise );
        PRINT_LOG_INFO(LOG_DEBUG, "x = %.16e\t%.16e\t%.16e", x->data[0], x->data[1], x->data[2] );
        PRINT_LOG_INFO(LOG_DEBUG, "p = %.16e\t%.16e\t%.16e", p->data[0], p->data[1], p->data[2] );
        PRINT_LOG_INFO(LOG_DEBUG, "sStar = %.16e\t%.16e\t%.16e", sigmaStar->data[0],
            sigmaStar->data[1], sigmaStar->data[2] );
        PRINT_LOG_INFO(LOG_DEBUG, "sKerr = %.16e\t%.16e\t%.16e", sigmaKerr->data[0],
            sigmaKerr->data[1], sigmaKerr->data[2] );
        PRINT_LOG_INFO(LOG_DEBUG,"csi = %.16e, Q = %.16e, pvr = %.16e, xi2 = %.16e , deltaT = %.16e, rho2 = %.16e, Lambda = %.16e, pxir = %.16e, B = %.16e", csi,Q, pvr, xi2, deltaT, rho2, Lambda, pxir, B );
        PRINT_LOG_INFO(LOG_DEBUG, "alpha = %.16e, beta_f = %.16e, gamma_rr = %.16e, gamma_thth = %.16e, gamma_ff = %.16e",
            sqrt((rho2*deltaT) * invLambda), ww*invLambda, deltaR*invrho2, invrho2, rho2*invLambda*invxi2);
        PRINT_LOG_INFO(LOG_DEBUG, "w = %.16e, wr = %.16e, wcos = %.16e", w, wr, wcos);

        PRINT_LOG_INFO(LOG_DEBUG, "ehat = (%.16e, %.16e, %.16e)", e3_x, e3_y, e3_z);
        PRINT_LOG_INFO(LOG_DEBUG, "xihat = (%.16e, %.16e, %.16e)", xi_x, xi_y, xi_z);
        PRINT_LOG_INFO(LOG_DEBUG, "vhat = (%.16e, %.16e, %.16e)", vx, vy, vz);
        PRINT_LOG_INFO(LOG_DEBUG, "Lhat = (%.16e, %.16e, %.16e)", Lhat[0], Lhat[1], Lhat[2]);

        PRINT_LOG_INFO(LOG_DEBUG, "KK = %.16e, costheta = %.16e", coeffs->KK, costheta);
        PRINT_LOG_INFO(LOG_DEBUG, "bulk = %.16e, logTerms = %.16e", bulk, logTerms );
        PRINT_LOG_INFO(LOG_DEBUG,"csi(miami) = %.16e", csi);
        PRINT_LOG_INFO(LOG_DEBUG, "(deltaU, bulk, logTerms, log arg) = (%.16e, %.16e, %.16e, %.16e)", deltaU, bulk, logTerms, 1. + coeffs->k1*u + coeffs->k2*u2 + coeffs->k3*u3 + coeffs->k4*u4
                + coeffs->k5*u5 + coeffs->k5l*u5*logu);
        // PRINT_LOG_INFO(LOG_DEBUG, "tmpBR = (%.16e, %.16e, %.16e)\n",  
        //     (-deltaT*invsqrtdeltaR + deltaT_r*0.5)*invsqrtdeltaT,
        //      (-deltaT*invsqrtdeltaR)*invsqrtdeltaT,
        //      (deltaT_r*0.5)*invsqrtdeltaT);

        PRINT_LOG_INFO(LOG_DEBUG, "a = %.16e, r = %.16e, S_con = %.16e", a, r, S_con);
        PRINT_LOG_INFO(LOG_DEBUG, "D = %.16e, ww = %.16e, rho = %.16e, Lambda = %.16e, xi = %.16e pr = %.16e, pf = %.16e, deltaR = %.16e, deltaT = %.16e",
            D, ww, sqrt(rho2), Lambda, sqrt(xi2), pr, pf, deltaR, deltaT );
        PRINT_LOG_INFO(LOG_DEBUG, "pr = %.16e, prT = %.16e", pr, prT );
        PRINT_LOG_INFO(LOG_DEBUG, "pxir = %.16e, pvr = %.16e", pxir, pvr );
        PRINT_LOG_INFO(LOG_DEBUG, "pn2 = %.16e, pp = %.16e, ptheta2 = %.16e", pn2, pp, ptheta2);
        PRINT_LOG_INFO(LOG_DEBUG, "deltaSigmaStar_x = %.16e, deltaSigmaStar_y = %.16e, deltaSigmaStar_z = %.16e",
            deltaSigmaStar_x, deltaSigmaStar_y, deltaSigmaStar_z );

        PRINT_LOG_INFO(LOG_DEBUG, "term 1 in Hns: %.16e",  prT*prT*prT*prT*qq*u2 );
        PRINT_LOG_INFO(LOG_DEBUG, "term 2 in Hns: %.16e", ptheta2/rho2 );
        PRINT_LOG_INFO(LOG_DEBUG, "term 3 in Hns = %.16e", pf*pf*rho2/(Lambda*xi2) );
        PRINT_LOG_INFO(LOG_DEBUG, "term 4 in Hns = %.16e", pr*pr*deltaR/rho2 );
        PRINT_LOG_INFO(LOG_DEBUG, "term 5 in Hns = %.16e", Lambda/(rho2*deltaT) );
        PRINT_LOG_INFO(LOG_DEBUG, "term 6 in Hns = %.16e", pf*ww/Lambda );

        PRINT_LOG_INFO(LOG_DEBUG, "expMu = %.16e, expNu = %.16e\n", expMU, expnu);
        PRINT_LOG_INFO(LOG_DEBUG, "SStar = (%.16e, %.16e, %.16e)\n", sx, sy, sz);

        PRINT_LOG_INFO(LOG_DEBUG, "Hns = %.16e, Hs = %.16e, Hss = %.16e\n", Hns, Hs, Hss );
        PRINT_LOG_INFO(LOG_DEBUG, "BR = %.16e, sv = %.16e, sn = %.16e, sxi = %.16e, nur = %.16e, mur = %.16e, nucos = %.16e, mucos = %.16e\n", 
            BR, sv, sn, sxi, nur, mur, nucos, mucos );
        PRINT_LOG_INFO(LOG_DEBUG, "HSOL = %.16e, HSONL = %.16e, Hwr = %.16e, Hwcos = %.16e\n", HSOL, HSONL, Hwr, Hwcos );
        PRINT_LOG_INFO(LOG_DEBUG, "H = %.16e", H );

        PRINT_LOG_INFO(LOG_DEBUG ,"Done printing components.");
        PRINT_LOG_INFO( LOG_CRITICAL, "Hreal = %.16e in Hamiltonian", Hreal);
        return REAL8_FAIL_NAN;
    }

    return Hreal;
}

REAL8 XLALSimIMRSpinPrecEOBHamiltonianDeltaT(
        SpinEOBHCoeffs *coeffs, /**<< Pre-computed coefficients which appear in the function */
        const REAL8    r,       /**<< Current orbital radius (in units of total mass) */
        const REAL8    eta,     /**<< Symmetric mass ratio */
        const REAL8    a        /**<< Normalized deformed Kerr spin */
        )
{

    REAL8 a2;
    REAL8 u, u2, u3, u4, u5;
    REAL8 m1PlusetaKK;

    REAL8 bulk;
    REAL8 logTerms;
    REAL8 deltaU;
    REAL8 deltaT;

    u  = 1./r;
    u2 = u*u;
    u3 = u2*u;
    u4 = u2*u2;
    u5 = u4*u;

    a2 = a*a;

    m1PlusetaKK = -1. + eta * coeffs->KK;

    bulk = 1./(m1PlusetaKK*m1PlusetaKK) + (2.*u)/m1PlusetaKK + a2*u2;

    logTerms = 1. + eta*coeffs->k0 + eta*log(fabs(1. + coeffs->k1*u + coeffs->k2*u2 + coeffs->k3*u3 + coeffs->k4*u4
                                                + coeffs->k5*u5 + coeffs->k5l*u5*log(u)));
    deltaU = bulk*logTerms;
    deltaU = fabs(deltaU);
    deltaT = r*r*deltaU;

    return deltaT;
}

REAL8 XLALSimIMRSpinPrecEOBHamiltonianDeltaR(
        SpinEOBHCoeffs *coeffs, /**<< Pre-computed coefficients which appear in the function */
        const REAL8    r,       /**<< Current orbital radius (in units of total mass) */
        const REAL8    eta,     /**<< Symmetric mass ratio */
        const REAL8    a        /**<< Normalized deformed Kerr spin */
        )
{
    REAL8 u2, u3;
    REAL8 D;
    REAL8 deltaT; /* The potential function, not a time interval... */
    REAL8 deltaR;

    u2 = 1./(r*r);
    u3 = u2 / r;

    D = 1. + log(1. + 6.*eta*u2 + 2.*(26. - 3.*eta)*eta*u3);

    deltaT = XLALSimIMRSpinPrecEOBHamiltonianDeltaT( coeffs, r, eta, a );

    deltaR = deltaT*D;
    return deltaR;
}
/**
 * Wrapper for GSL to call the Hamiltonian function
 */
REAL8
GSLSpinPrecHamiltonianWrapper(double x, void *params)
{
	HcapDerivParams *dParams = (HcapDerivParams *) params;

	SpinEOBParams      *seobParams = (SpinEOBParams *) dParams->params;
	SpinEOBHCoeffs *coeffs = (SpinEOBHCoeffs *) dParams->params->seobCoeffs;

	REAL8		tmpVec   [12];
	REAL8		s1normData[3], s2normData[3], sKerrData[3], sStarData[3];

	/*
	 * These are the vectors which will be used in the call to the
	 * Hamiltonian
	 */
	REAL8Vector	r  , p, spin1, spin2, spin1norm, spin2norm;
	REAL8Vector	sigmaKerr, sigmaStar;

	INT4		i;
	REAL8		a;
	REAL8		m1 = seobParams->m1;
	REAL8		m2 = seobParams->m2;
	REAL8 m_total = m1 + m2;
	m1 /=m_total;
	m2 /=m_total;
	REAL8 mT2 = (m1 + m2) * (m1 + m2);
	REAL8 eta = m1 * m2 / mT2;

	INT4	oldTortoise = dParams->params->tortoise;
	/* Use a temporary vector to avoid corrupting the main function */
	memcpy(tmpVec, dParams->values, sizeof(tmpVec));

	/* Set the relevant entry in the vector to the correct value */
	tmpVec[dParams->varyParam] = x;

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

	memcpy(spin1norm.data, tmpVec + 6, 3 * sizeof(REAL8));
	memcpy(spin2norm.data, tmpVec + 9, 3 * sizeof(REAL8));

	/*
	 * To compute the SigmaKerr and SigmaStar, we need the non-normalized
	 * spin values, i.e. S_i. The spins being evolved are S_i/M^2.
	 */
	for (i = 0; i < 3; i++) 
    {
		spin1.data[i] *= mT2;
		spin2.data[i] *= mT2;
	}

	/* Calculate various spin parameters */
    /* Note that XLALSimIMRSpinEOBCalculateSigmaKerr and XLALSimIMRSpinEOBCalculateSigmaStar return a unitless quantity */
	EOBCalculateSigmaKerr(&sigmaKerr, &spin1norm, &spin2norm);
	EOBCalculateSigmaStar(&sigmaStar, m1,
					    m2, &spin1norm, &spin2norm);
	a = sqrt(sigmaKerr.data[0] * sigmaKerr.data[0]
		 + sigmaKerr.data[1] * sigmaKerr.data[1]
		 + sigmaKerr.data[2] * sigmaKerr.data[2]);

	if (isnan(a) ) 
    {
          PRINT_LOG_INFO(LOG_CRITICAL, "a = nan  \n");
          return REAL8_FAIL_NAN;
	}

	REAL8 SpinEOBH=0.0;
	SpinEOBH = XLALSimIMRSpinPrecEOBHamiltonian(seobParams->eta, &r, &p, 
        &spin1norm, &spin2norm, &sigmaKerr, &sigmaStar, 
        dParams->params->tortoise, dParams->params->seobCoeffs, dParams->params->hParams) / seobParams->eta;
	if (dParams->varyParam < 3)
		dParams->params->tortoise = oldTortoise;
	return SpinEOBH;
}


REAL8
XLALSpinPrecHcapNumDerivWRTParam(
			     const INT paramIdx,	/**<< Index of the parameters */
			     const REAL8 values[],	/**<< Dynamical variables */
			     SpinEOBParams * funcParams	/**<< EOB Parameters */
)
{
	static const REAL8 STEP_SIZE = 2.0e-3;

	HcapDerivParams	params;

	REAL8		result;

	gsl_function	F;
	INT		gslStatus;
    UINT i;
	REAL8		mass1   , mass2;

	/* The error in a derivative as measured by GSL */
	REAL8		absErr;

	/* Set up pointers for GSL */
	params.params = funcParams;


	F.function = &GSLSpinPrecHamiltonianWrapper;
	F.params = &params;
	params.varyParam = paramIdx;

	mass1 = params.params->m1;
	mass2 = params.params->m2;

    REAL8 mT2 = (mass1 + mass2) * (mass1 + mass2);
    REAL8 tmpValues[14];
    for ( i = 0; i < 3; i++) {
        tmpValues[i] = values[i];
        tmpValues[i + 3] = values[i + 3];
        tmpValues[i + 6] = values[i + 6]/mT2;
        tmpValues[i + 9] = values[i + 9]/mT2;
    }
    tmpValues[12] = values[12];
    tmpValues[13] = values[13];
    params.values = tmpValues;
    REAL8 m_total = mass1 + mass2;
    mass1 /=m_total;
    mass2 /=m_total;
    params.params->seobCoeffs->updateHCoeffs = 1;
	/* Now calculate derivatives w.r.t. the required parameter */
    if (paramIdx >=0 && paramIdx < 6) {
        gslStatus = gsl_deriv_central(&F, values[paramIdx], STEP_SIZE, &result, &absErr);
    }
    else if (paramIdx >= 6 && paramIdx < 9) {
      params.params->seobCoeffs->updateHCoeffs = 1;
		gslStatus = gsl_deriv_central(&F, values[paramIdx],
			      STEP_SIZE * mass1 * mass1, &result, &absErr);
	}
    else if (paramIdx >= 9 && paramIdx < 12) {
      params.params->seobCoeffs->updateHCoeffs = 1;
		gslStatus = gsl_deriv_central(&F, values[paramIdx],
			      STEP_SIZE * mass2 * mass2, &result, &absErr);
	}
    else {
        PRINT_LOG_INFO(LOG_CRITICAL, "Requested partial derivative is not availabble");
        return REAL8_FAIL_NAN;
    }
	if (gslStatus != GSL_SUCCESS) {
		PRINT_LOG_INFO(LOG_CRITICAL, "Failure in GSL function");
		return REAL8_FAIL_NAN;
	}
	//XLAL_PRINT_INFO("Abserr = %e\n", absErr);
	return result;
}

int	XLALSpinPrecHcapNumericalDerivative(
                    double	t,	/**<< UNUSED */
                    const	REAL8	values[],	/**<< Dynamical variables */
                    REAL8	dvalues[],	/**<< Time derivatives of variables (returned) */
                    void    *funcParams	/**<< EOB parameters */
)
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
	REAL8		rData    [3], pData[3];

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
    REAL8  	eobD_r, deltaU_u, deltaU_r;
	REAL8		dcsi, csi;

	REAL8		tmpValues[12];
	REAL8		Tmatrix  [3][3], invTmatrix[3][3], dTijdXk[3][3][3];
	REAL8		tmpPdotT1[3], tmpPdotT2[3], tmpPdotT3[3];
	//3 terms of Eq.A5

	/* NQC coefficients container */
    EOBNonQCCoeffs * nqcCoeffs = NULL;

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
	if(SpinAlignedEOBversion==4)
    {
	  STEP_SIZE = 2.0e-3; //Allow a different step size for v4P 
	}
	else
    {
	  STEP_SIZE = 2.0e-4;
	}
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
	rMag = sqrt(rVec.data[0] * rVec.data[0] + rVec.data[1] * rVec.data[1] + rVec.data[2] * rVec.data[2]);
    /* This is p^*.r/|r| */
	prT = pVec.data[0] * (rVec.data[0] / rMag) + pVec.data[1] * (rVec.data[1] / rMag)
		+ pVec.data[2] * (rVec.data[2] / rMag);

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

#if 0
	if (debugPK) {
		XLAL_PRINT_INFO("csi = %.12e\n", csi);
		for (i = 0; i < 3; i++)
			XLAL_PRINT_INFO("p,p*: %.12e\t%.12e\n", pData[i], tmpP[i]);
	}
#endif
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

	//Print out the T - matrix for comparison
#if 0
		if (debugPK) {
			XLAL_PRINT_INFO("\nT-Matrix:\n");
			for (i = 0; i < 3; i++)
				XLAL_PRINT_INFO("%le\t%le\t%le\n", Tmatrix[i][0], Tmatrix[i][1], Tmatrix[i][2]);

			for (i = 0; i < 3; i++) {
				XLAL_PRINT_INFO("dT[%d][j]dX[k]:\n", i);
				for (j = 0; j < 3; j++)
					XLAL_PRINT_INFO("%.12e\t%.12e\t%.12e\n", dTijdXk[i][j][0],
					dTijdXk[i][j][1], dTijdXk[i][j][2]);
				XLAL_PRINT_INFO("\n");
			}
		}
#endif
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
	Lx = values[1] * values[5] - values[2] * values[4];
	Ly = values[2] * values[3] - values[0] * values[5];
	Lz = values[0] * values[4] - values[1] * values[3];

	magL = sqrt(Lx * Lx + Ly * Ly + Lz * Lz);

	Lhatx = Lx / magL;
	Lhaty = Ly / magL;
	Lhatz = Lz / magL;

	/* Calculate the polar data */
	polarDynamics.length = 4;
	polarDynamics.data = polData;

	r = polarDynamics.data[0] = sqrt(values[0] * values[0] + values[1] * values[1]
			      + values[2] * values[2]);
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
    // if (r>24.5)
    //     print_debug("get chi1 = (%g, %g, %g), chi2 = (%g, %g, %g)\n", 
    //         c1x, c1y, c1z, c2x, c2y, c2z);

	//
	// Compute \vec{L_N} = \vec{r} \times \.{\vec{r}}, \vec{S_i} \dot
	// \vec{L_N} and chiS and chiA
	//

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
#if 0
	if (debugPK) {
		XLAL_PRINT_INFO("chiS = %.12e, chiA = %.12e\n", chiS, chiA);
		fflush(NULL);
	}
	if (debugPK) {
		XLAL_PRINT_INFO("Computing derivatives for values\n%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n\n",
		    (double)values[0], (double)values[1], (double)values[2],
		    (double)values[3], (double)values[4], (double)values[5],
            (double)values[6], (double)values[7], (double)values[8],
            (double)values[9], (double)values[10], (double)values[11]);
		XLAL_PRINT_INFO("tmpDvalues\n%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t\n",
		       (double)tmpDValues[0], (double)tmpDValues[1], (double)tmpDValues[2],
		       (double)tmpDValues[3], (double)tmpDValues[4], (double)tmpDValues[5],
		       (double)tmpDValues[6], (double)tmpDValues[7], (double)tmpDValues[8],
		       (double)tmpDValues[9], (double)tmpDValues[10], (double)tmpDValues[11]);
	}
#endif
	/*
	 * Compute the test-particle limit spin of the deformed-Kerr
	 * background
	 */
	switch (SpinAlignedEOBversion) {
	case 1:
        /* See below Eq. 17 of PRD 86, 041011 (2012) */
		tplspin = 0.0;
		break;
	case 2:
	case 3:
    case 4:
        /* See below Eq. 4 of PRD 89, 061502(R) (2014)*/
		tplspin = (1. - 2. * eta) * chiS + (mass1 - mass2) / (mass1 + mass2) * chiA;
        break;
	default:
		PRINT_LOG_INFO(LOG_INFO, "Unknown SEOBNR version! At present only v1 and v2 are available.\n");
		return CEV_FAILURE;
		break;
	}


    memcpy(params.params->s1Vec->data, s1norm.data, 3*sizeof(*params.params->s1Vec->data));
    memcpy(params.params->s2Vec->data, s2norm.data, 3*sizeof(*params.params->s2Vec->data));
    memcpy(params.params->sigmaStar->data, sStar.data, 3*sizeof(*params.params->sigmaStar->data));
    memcpy(params.params->sigmaKerr->data, sKerr.data, 3*sizeof(*params.params->sigmaKerr->data));


	params.params->a = a;
    if (params.params->alignedSpins==1) 
    {
        if (CODE_VERSION == 3)
            EccPrec_CalcSpinPrecFacWaveformCoefficients(params.params->hCoeffs, mass1, mass2, eta, tplspin,
                            chiS, chiA, csx, csy, csz, cax, cay, caz, 451);
        else
            XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients(
                params.params->hCoeffs, mass1, mass2, eta, tplspin,
                            chiS, chiA, SpinAlignedEOBversion);
    }
    else 
    {
			/* This is needed because SpinAlignedEOBversion is set to 2 for v3 */
			/* while XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients requires 3 ... */
        if ( SpinAlignedEOBversion == 2 ) SpinAlignedEOBversionForWaveformCoefficients = 3;
        else SpinAlignedEOBversionForWaveformCoefficients = SpinAlignedEOBversion;

        // if ( params.params->use_hm ) SpinAlignedEOBversionForWaveformCoefficients = 451;
        // else SpinAlignedEOBversionForWaveformCoefficients = SpinAlignedEOBversion;
        if (CODE_VERSION == 3)
            EccPrec_CalcSpinPrecFacWaveformCoefficients(params.params->hCoeffs, mass1, mass2, eta, tplspin,
                            chiS, chiA, csx, csy, csz, cax, cay, caz, 451);
        else
            XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients(params.params->hCoeffs, mass1, mass2, eta, tplspin,
                chiS, chiA, SpinAlignedEOBversionForWaveformCoefficients);
    }
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

	H = H * (mass1 + mass2);

	/* Now we have the ingredients to compute the flux */
	memcpy(tmpValues, values, 12 * sizeof(REAL8));
	cartDynamics.data = tmpValues;
#if 0
	if (debugPK) {
		XLAL_PRINT_INFO("params.params->a = %.12e, %.12e\n", a, params.params->a);
		fflush(NULL);
	}
#endif
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
	/*
	 * Looking at consistency with the non-spinning model, we have to divide the
	 * flux by eta
	 */
	flux = flux / eta;

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
#if 0
	if (debugPK) {
		XLAL_PRINT_INFO("omega = %.12e \n flux = %.12e \n Lmag = %.12e\n", omega, flux, magL);
		XLAL_PRINT_INFO("rrForce = %.12e %.12e %.12e\n", -flux * values[3] / (omega * magL), -flux * values[4] / (omega * magL), -flux * values[5] / (omega * magL));
	}
#endif
	/* Now pDot */
	/* Compute the first and second terms in eq. A5 of PRD 81, 084041 (2010) */
# if 1
	for (i = 0; i < 3; i++) 
    {
		for (j = 0, tmpPdotT1[i] = 0.; j < 3; j++)
			tmpPdotT1[i] += -tmpDValues[j] * Tmatrix[i][j];
		tmpPdotT2[i] = -flux * values[i + 3] / (omega * magL);
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

#if 0
	if (debugPK) {
		for (i = 0; i < 3; i++)
			for (j = 0; j < 1; j++)
				for (l = 0; l < 3; l++) {
					double		sum = 0;
					for (k = 0; k < 3; k++)
						sum += dTijdXk[i][k][j] * invTmatrix[k][l];
					XLAL_PRINT_INFO("\n sum[%d][%d][%d] = %.12e", i, j, l, sum);

				}
		XLAL_PRINT_INFO("\n\n Printing dTdX * Tinverse:\n");
		for (l = 0; l < 3; l++) {
			for (i = 0; i < 3; i++)
				for (j = 0; j < 3; j++) {
					double		sum = 0;
					for (k = 0; k < 3; k++) {
						sum += dTijdXk[i][k][l] * invTmatrix[k][j];
						XLAL_PRINT_INFO("\n sum[%d][%d][%d] = %.12e", l, i, j, sum);
					}
				}
		}
	}
	if (debugPK)
		XLAL_PRINT_INFO("\npData: {%.12e, %.12e, %.12e}\n", pData[0], pData[1], pData[2]);
#endif
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

	/* Add them to obtain pDot */
	for (i = 0; i < 3; i++)
		dvalues[i + 3] = tmpPdotT1[i] + tmpPdotT2[i] + tmpPdotT3[i];

#if 0
	if (debugPK) {
		XLAL_PRINT_INFO("\ntmpPdotT3 = ");
		for (i = 0; i < 3; i++)
			XLAL_PRINT_INFO("%.12e ", tmpPdotT3[i]);
		XLAL_PRINT_INFO("\n");
	}
#endif

    /* Eqs. 11c-11d of PRD 89, 084006 (2014) */
	/* spin1 */
#if 0
	if (debugPK) {
		XLAL_PRINT_INFO("Raw spin1 derivatives = %.12e %.12e %.12e\n", tmpDValues[6], tmpDValues[7], tmpDValues[8]);
		XLAL_PRINT_INFO("Raw spin2 derivatives = %.12e %.12e %.12e\n", tmpDValues[9], tmpDValues[10], tmpDValues[11]);
	}
#endif
     /* The factor eta is there becasue Hreal is normalized by eta */
	dvalues[6] = eta * (tmpDValues[7] * values[8] - tmpDValues[8] * values[7]);
	dvalues[7] = eta * (tmpDValues[8] * values[6] - tmpDValues[6] * values[8]);
	dvalues[8] = eta * (tmpDValues[6] * values[7] - tmpDValues[7] * values[6]);

	/* spin2 */
    /* The factor eta is there becasue Hreal is normalized by eta */
	dvalues[9] = eta * (tmpDValues[10] * values[11] - tmpDValues[11] * values[10]);
	dvalues[10] = eta * (tmpDValues[11] * values[9] - tmpDValues[9] * values[11]);
	dvalues[11] = eta * (tmpDValues[9] * values[10] - tmpDValues[10] * values[9]);

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

#if 0
	if (debugPK) {
    XLAL_PRINT_INFO("\nIn XLALSpinPrecHcapNumericalDerivative:\n");
		/* Print out all mass parameters */
		XLAL_PRINT_INFO("m1 = %12.12lf, m2 = %12.12lf, eta = %12.12lf\n",
          (double)mass1, (double)mass2, (double)eta);
		/* Print out all spin parameters */
		XLAL_PRINT_INFO("spin1 = {%12.12lf,%12.12lf,%12.12lf}, spin2 = {%12.12lf,%12.12lf,%12.12lf}\n",
            (double)s1.data[0], (double)s1.data[1], (double)s1.data[2],
            (double)s2.data[0], (double)s2.data[1], (double)s2.data[2]);
		XLAL_PRINT_INFO("sigmaStar = {%12.12lf,%12.12lf,%12.12lf}, sigmaKerr = {%12.12lf,%12.12lf,%12.12lf}\n",
		       (double)sStar.data[0], (double)sStar.data[1],
		       (double)sStar.data[2], (double)sKerr.data[0],
		       (double)sKerr.data[1], (double)sKerr.data[2]);
		XLAL_PRINT_INFO("L = {%12.12lf,%12.12lf,%12.12lf}, |L| = %12.12lf\n",
        (double)Lx, (double)Ly, (double)Lz, (double)magL);
		XLAL_PRINT_INFO("dLdt = {%12.12lf,%12.12lf,%12.12lf}, d|L|dt = %12.12lf\n",
        (double)dLx, (double)dLy, (double)dLz, (double)dMagL);
		XLAL_PRINT_INFO("Polar coordinates = {%12.12lf, %12.12lf, %12.12lf, %12.12lf}\n",
		       (double)polData[0], (double)polData[1], (double)polData[2],
           (double)polData[3]);

		XLAL_PRINT_INFO("Cartesian coordinates: {%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf}\n",
		       (double)values[0], (double)values[1], (double)values[2],
           (double)values[3], (double)values[4], (double)values[5],
           (double)values[6], (double)values[7], (double)values[8],
           (double)values[9], (double)values[10], (double)values[11],
		       (double)values[12], (double)values[13]);
		XLAL_PRINT_INFO("Cartesian derivatives: {%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf}\n",
		       (double)dvalues[0], (double)dvalues[1], (double)dvalues[2],
           (double)dvalues[3], (double)dvalues[4], (double)dvalues[5],
           (double)dvalues[6], (double)dvalues[7], (double)dvalues[8],
           (double)dvalues[9], (double)dvalues[10], (double)dvalues[11],
           (double)dvalues[12], (double)dvalues[13]);

     XLAL_PRINT_INFO("Hamiltonian = %12.12lf, Flux = %12.12lf, Omega = %12.12lf\n",
              H/ (mass1 + mass2), eta*flux, omega);
		fflush(NULL);
	}

    if(debugPK){
    for(i=0; i<14; i++)
    if(dvalues[i] > 1e3)
    {
        XLAL_PRINT_INFO("\nIn XLALSpinPrecHcapNumericalDerivative:\n");
        XLAL_PRINT_INFO("Derivatives have blown up!\n");
        for(j=0; j<14; XLAL_PRINT_INFO("dvalues[%d] = %3.12f\n", j, dvalues[j]), j++);
        XLAL_PRINT_INFO("Flux = %3.12f\n\n", flux);
        break;
        }
    }
#endif
    return CEV_SUCCESS;
}


int	XLALSpinPrecHcapNumericalDerivative_Conserve(
                    double	t,	/**<< UNUSED */
                    const	REAL8	values[],	/**<< Dynamical variables */
                    REAL8	dvalues[],	/**<< Time derivatives of variables (returned) */
                    void    *funcParams	/**<< EOB parameters */
)
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
	REAL8		rData    [3], pData[3];

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
    REAL8  	eobD_r, deltaU_u, deltaU_r;
	REAL8		dcsi, csi;

	REAL8		tmpValues[12];
	REAL8		Tmatrix  [3][3], invTmatrix[3][3], dTijdXk[3][3][3];
	REAL8		tmpPdotT1[3], tmpPdotT2[3], tmpPdotT3[3];
	//3 terms of Eq.A5

	/* NQC coefficients container */
    EOBNonQCCoeffs * nqcCoeffs = NULL;

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
	if(SpinAlignedEOBversion==4)
    {
	  STEP_SIZE = 2.0e-3; //Allow a different step size for v4P 
	}
	else
    {
	  STEP_SIZE = 2.0e-4;
	}
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
	rMag = sqrt(rVec.data[0] * rVec.data[0] + rVec.data[1] * rVec.data[1] + rVec.data[2] * rVec.data[2]);
    /* This is p^*.r/|r| */
	prT = pVec.data[0] * (rVec.data[0] / rMag) + pVec.data[1] * (rVec.data[1] / rMag)
		+ pVec.data[2] * (rVec.data[2] / rMag);

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

#if 0
	if (debugPK) {
		XLAL_PRINT_INFO("csi = %.12e\n", csi);
		for (i = 0; i < 3; i++)
			XLAL_PRINT_INFO("p,p*: %.12e\t%.12e\n", pData[i], tmpP[i]);
	}
#endif
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

	//Print out the T - matrix for comparison
#if 0
		if (debugPK) {
			XLAL_PRINT_INFO("\nT-Matrix:\n");
			for (i = 0; i < 3; i++)
				XLAL_PRINT_INFO("%le\t%le\t%le\n", Tmatrix[i][0], Tmatrix[i][1], Tmatrix[i][2]);

			for (i = 0; i < 3; i++) {
				XLAL_PRINT_INFO("dT[%d][j]dX[k]:\n", i);
				for (j = 0; j < 3; j++)
					XLAL_PRINT_INFO("%.12e\t%.12e\t%.12e\n", dTijdXk[i][j][0],
					dTijdXk[i][j][1], dTijdXk[i][j][2]);
				XLAL_PRINT_INFO("\n");
			}
		}
#endif
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
	Lx = values[1] * values[5] - values[2] * values[4];
	Ly = values[2] * values[3] - values[0] * values[5];
	Lz = values[0] * values[4] - values[1] * values[3];

	magL = sqrt(Lx * Lx + Ly * Ly + Lz * Lz);

	Lhatx = Lx / magL;
	Lhaty = Ly / magL;
	Lhatz = Lz / magL;

	/* Calculate the polar data */
	polarDynamics.length = 4;
	polarDynamics.data = polData;

	r = polarDynamics.data[0] = sqrt(values[0] * values[0] + values[1] * values[1]
			      + values[2] * values[2]);
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
    // if (r>24.5)
    //     print_debug("get chi1 = (%g, %g, %g), chi2 = (%g, %g, %g)\n", 
    //         c1x, c1y, c1z, c2x, c2y, c2z);

	//
	// Compute \vec{L_N} = \vec{r} \times \.{\vec{r}}, \vec{S_i} \dot
	// \vec{L_N} and chiS and chiA
	//

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
#if 0
	if (debugPK) {
		XLAL_PRINT_INFO("chiS = %.12e, chiA = %.12e\n", chiS, chiA);
		fflush(NULL);
	}
	if (debugPK) {
		XLAL_PRINT_INFO("Computing derivatives for values\n%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n\n",
		    (double)values[0], (double)values[1], (double)values[2],
		    (double)values[3], (double)values[4], (double)values[5],
            (double)values[6], (double)values[7], (double)values[8],
            (double)values[9], (double)values[10], (double)values[11]);
		XLAL_PRINT_INFO("tmpDvalues\n%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t\n",
		       (double)tmpDValues[0], (double)tmpDValues[1], (double)tmpDValues[2],
		       (double)tmpDValues[3], (double)tmpDValues[4], (double)tmpDValues[5],
		       (double)tmpDValues[6], (double)tmpDValues[7], (double)tmpDValues[8],
		       (double)tmpDValues[9], (double)tmpDValues[10], (double)tmpDValues[11]);
	}
#endif
	/*
	 * Compute the test-particle limit spin of the deformed-Kerr
	 * background
	 */
	switch (SpinAlignedEOBversion) {
	case 1:
        /* See below Eq. 17 of PRD 86, 041011 (2012) */
		tplspin = 0.0;
		break;
	case 2:
	case 3:
    case 4:
        /* See below Eq. 4 of PRD 89, 061502(R) (2014)*/
		tplspin = (1. - 2. * eta) * chiS + (mass1 - mass2) / (mass1 + mass2) * chiA;
        break;
	default:
		PRINT_LOG_INFO(LOG_INFO, "Unknown SEOBNR version! At present only v1 and v2 are available.\n");
		return CEV_FAILURE;
		break;
	}


    memcpy(params.params->s1Vec->data, s1norm.data, 3*sizeof(*params.params->s1Vec->data));
    memcpy(params.params->s2Vec->data, s2norm.data, 3*sizeof(*params.params->s2Vec->data));
    memcpy(params.params->sigmaStar->data, sStar.data, 3*sizeof(*params.params->sigmaStar->data));
    memcpy(params.params->sigmaKerr->data, sKerr.data, 3*sizeof(*params.params->sigmaKerr->data));


	params.params->a = a;
    if (params.params->alignedSpins==1) 
    {
        if (CODE_VERSION == 3)
            EccPrec_CalcSpinPrecFacWaveformCoefficients(params.params->hCoeffs, mass1, mass2, eta, tplspin,
                            chiS, chiA, csx, csy, csz, cax, cay, caz, 451);
        else
            XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients(
                params.params->hCoeffs, mass1, mass2, eta, tplspin,
                            chiS, chiA, SpinAlignedEOBversion);
    }
    else 
    {
			/* This is needed because SpinAlignedEOBversion is set to 2 for v3 */
			/* while XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients requires 3 ... */
        if ( SpinAlignedEOBversion == 2 ) SpinAlignedEOBversionForWaveformCoefficients = 3;
        else SpinAlignedEOBversionForWaveformCoefficients = SpinAlignedEOBversion;

        // if ( params.params->use_hm ) SpinAlignedEOBversionForWaveformCoefficients = 451;
        // else SpinAlignedEOBversionForWaveformCoefficients = SpinAlignedEOBversion;
        if (CODE_VERSION == 3)
            EccPrec_CalcSpinPrecFacWaveformCoefficients(params.params->hCoeffs, mass1, mass2, eta, tplspin,
                            chiS, chiA, csx, csy, csz, cax, cay, caz, 451);
        else
            XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients(params.params->hCoeffs, mass1, mass2, eta, tplspin,
                chiS, chiA, SpinAlignedEOBversionForWaveformCoefficients);
    }
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

	H = H * (mass1 + mass2);

	/* Now we have the ingredients to compute the flux */
	memcpy(tmpValues, values, 12 * sizeof(REAL8));
	cartDynamics.data = tmpValues;
#if 0
	if (debugPK) {
		XLAL_PRINT_INFO("params.params->a = %.12e, %.12e\n", a, params.params->a);
		fflush(NULL);
	}
#endif
    /* Eq. 13 of PRD 86, 024011 (2012) */
    if ( 0 ) 
    {
        flux = XLALInspiralPrecSpinFactorizedFlux(&polarDynamics, &cartDynamics,
						  nqcCoeffs, omega, dr, ncrv, params.params, H / (mass1 + mass2), lMax, SpinAlignedEOBversion);
    }
    else if (1) 
    {
            flux = 0.;
    }
    else 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Wrong ignorflux option in XLALSpinPrecHcapNumericalDerivative!");
        return CEV_FAILURE;
    }
	/*
	 * Looking at consistency with the non-spinning model, we have to divide the
	 * flux by eta
	 */
	flux = flux / eta;

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
#if 0
	if (debugPK) {
		XLAL_PRINT_INFO("omega = %.12e \n flux = %.12e \n Lmag = %.12e\n", omega, flux, magL);
		XLAL_PRINT_INFO("rrForce = %.12e %.12e %.12e\n", -flux * values[3] / (omega * magL), -flux * values[4] / (omega * magL), -flux * values[5] / (omega * magL));
	}
#endif
	/* Now pDot */
	/* Compute the first and second terms in eq. A5 of PRD 81, 084041 (2010) */
# if 1
	for (i = 0; i < 3; i++) 
    {
		for (j = 0, tmpPdotT1[i] = 0.; j < 3; j++)
			tmpPdotT1[i] += -tmpDValues[j] * Tmatrix[i][j];
		tmpPdotT2[i] = -flux * values[i + 3] / (omega * magL);
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

#if 0
	if (debugPK) {
		for (i = 0; i < 3; i++)
			for (j = 0; j < 1; j++)
				for (l = 0; l < 3; l++) {
					double		sum = 0;
					for (k = 0; k < 3; k++)
						sum += dTijdXk[i][k][j] * invTmatrix[k][l];
					XLAL_PRINT_INFO("\n sum[%d][%d][%d] = %.12e", i, j, l, sum);

				}
		XLAL_PRINT_INFO("\n\n Printing dTdX * Tinverse:\n");
		for (l = 0; l < 3; l++) {
			for (i = 0; i < 3; i++)
				for (j = 0; j < 3; j++) {
					double		sum = 0;
					for (k = 0; k < 3; k++) {
						sum += dTijdXk[i][k][l] * invTmatrix[k][j];
						XLAL_PRINT_INFO("\n sum[%d][%d][%d] = %.12e", l, i, j, sum);
					}
				}
		}
	}
	if (debugPK)
		XLAL_PRINT_INFO("\npData: {%.12e, %.12e, %.12e}\n", pData[0], pData[1], pData[2]);
#endif
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

	/* Add them to obtain pDot */
	for (i = 0; i < 3; i++)
		dvalues[i + 3] = tmpPdotT1[i] + tmpPdotT2[i] + tmpPdotT3[i];

#if 0
	if (debugPK) {
		XLAL_PRINT_INFO("\ntmpPdotT3 = ");
		for (i = 0; i < 3; i++)
			XLAL_PRINT_INFO("%.12e ", tmpPdotT3[i]);
		XLAL_PRINT_INFO("\n");
	}
#endif

    /* Eqs. 11c-11d of PRD 89, 084006 (2014) */
	/* spin1 */
#if 0
	if (debugPK) {
		XLAL_PRINT_INFO("Raw spin1 derivatives = %.12e %.12e %.12e\n", tmpDValues[6], tmpDValues[7], tmpDValues[8]);
		XLAL_PRINT_INFO("Raw spin2 derivatives = %.12e %.12e %.12e\n", tmpDValues[9], tmpDValues[10], tmpDValues[11]);
	}
#endif
     /* The factor eta is there becasue Hreal is normalized by eta */
	dvalues[6] = eta * (tmpDValues[7] * values[8] - tmpDValues[8] * values[7]);
	dvalues[7] = eta * (tmpDValues[8] * values[6] - tmpDValues[6] * values[8]);
	dvalues[8] = eta * (tmpDValues[6] * values[7] - tmpDValues[7] * values[6]);

	/* spin2 */
    /* The factor eta is there becasue Hreal is normalized by eta */
	dvalues[9] = eta * (tmpDValues[10] * values[11] - tmpDValues[11] * values[10]);
	dvalues[10] = eta * (tmpDValues[11] * values[9] - tmpDValues[9] * values[11]);
	dvalues[11] = eta * (tmpDValues[9] * values[10] - tmpDValues[10] * values[9]);

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

#if 0
	if (debugPK) {
    XLAL_PRINT_INFO("\nIn XLALSpinPrecHcapNumericalDerivative:\n");
		/* Print out all mass parameters */
		XLAL_PRINT_INFO("m1 = %12.12lf, m2 = %12.12lf, eta = %12.12lf\n",
          (double)mass1, (double)mass2, (double)eta);
		/* Print out all spin parameters */
		XLAL_PRINT_INFO("spin1 = {%12.12lf,%12.12lf,%12.12lf}, spin2 = {%12.12lf,%12.12lf,%12.12lf}\n",
            (double)s1.data[0], (double)s1.data[1], (double)s1.data[2],
            (double)s2.data[0], (double)s2.data[1], (double)s2.data[2]);
		XLAL_PRINT_INFO("sigmaStar = {%12.12lf,%12.12lf,%12.12lf}, sigmaKerr = {%12.12lf,%12.12lf,%12.12lf}\n",
		       (double)sStar.data[0], (double)sStar.data[1],
		       (double)sStar.data[2], (double)sKerr.data[0],
		       (double)sKerr.data[1], (double)sKerr.data[2]);
		XLAL_PRINT_INFO("L = {%12.12lf,%12.12lf,%12.12lf}, |L| = %12.12lf\n",
        (double)Lx, (double)Ly, (double)Lz, (double)magL);
		XLAL_PRINT_INFO("dLdt = {%12.12lf,%12.12lf,%12.12lf}, d|L|dt = %12.12lf\n",
        (double)dLx, (double)dLy, (double)dLz, (double)dMagL);
		XLAL_PRINT_INFO("Polar coordinates = {%12.12lf, %12.12lf, %12.12lf, %12.12lf}\n",
		       (double)polData[0], (double)polData[1], (double)polData[2],
           (double)polData[3]);

		XLAL_PRINT_INFO("Cartesian coordinates: {%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf}\n",
		       (double)values[0], (double)values[1], (double)values[2],
           (double)values[3], (double)values[4], (double)values[5],
           (double)values[6], (double)values[7], (double)values[8],
           (double)values[9], (double)values[10], (double)values[11],
		       (double)values[12], (double)values[13]);
		XLAL_PRINT_INFO("Cartesian derivatives: {%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf}\n",
		       (double)dvalues[0], (double)dvalues[1], (double)dvalues[2],
           (double)dvalues[3], (double)dvalues[4], (double)dvalues[5],
           (double)dvalues[6], (double)dvalues[7], (double)dvalues[8],
           (double)dvalues[9], (double)dvalues[10], (double)dvalues[11],
           (double)dvalues[12], (double)dvalues[13]);

     XLAL_PRINT_INFO("Hamiltonian = %12.12lf, Flux = %12.12lf, Omega = %12.12lf\n",
              H/ (mass1 + mass2), eta*flux, omega);
		fflush(NULL);
	}

    if(debugPK){
    for(i=0; i<14; i++)
    if(dvalues[i] > 1e3)
    {
        XLAL_PRINT_INFO("\nIn XLALSpinPrecHcapNumericalDerivative:\n");
        XLAL_PRINT_INFO("Derivatives have blown up!\n");
        for(j=0; j<14; XLAL_PRINT_INFO("dvalues[%d] = %3.12f\n", j, dvalues[j]), j++);
        XLAL_PRINT_INFO("Flux = %3.12f\n\n", flux);
        break;
        }
    }
#endif
    return CEV_SUCCESS;
}

/**
 * Wrapper for GSL to call the Hamiltonian function. This is simply the function
 * GSLSpinPrecHamiltonianWrapper copied over. The alternative was to make it non-static
 * which increases runtime as static functions can be better optimized.
 */
static double GSLSpinPrecHamiltonianWrapperFordHdpphi( double x, void *params )
{
  HcapDerivParams *dParams = (HcapDerivParams *)params;

  SpinEOBParams *seobParams = (SpinEOBParams*) dParams->params;
  SpinEOBHCoeffs *coeffs = (SpinEOBHCoeffs*) dParams->params->seobCoeffs;

  REAL8 tmpVec[12] = {0.};
  REAL8 rpolar[3] = {0.}, rcart[3] = {0.}, ppolar[3] = {0.}, pcart[3] = {0.};
  REAL8 s1normData[3] = {0.}, s2normData[3] = {0.}, sKerrData[3] = {0.}, sStarData[3] = {0.};

  /* These are the vectors which will be used in the call to the Hamiltonian */
  REAL8Vector r, p, spin1, spin2, spin1norm, spin2norm;
  REAL8Vector sigmaKerr, sigmaStar;

  INT4 i;
  REAL8 a;
  REAL8 m1 = seobParams->m1;
  REAL8 m2 = seobParams->m2;
  REAL8 mT2 = (m1+m2)*(m1+m2);
  REAL8 eta = m1*m2/mT2;

  /* Use a temporary vector to avoid corrupting the main function */
  memcpy( tmpVec, dParams->values, sizeof(tmpVec) );

  /* Set the relevant entry in the vector to the correct value */
  tmpVec[dParams->varyParam] = x;

  /* Set the LAL-style vectors to point to the appropriate things */
  r.length = p.length = spin1.length = spin2.length = spin1norm.length = spin2norm.length = 3;
  sigmaKerr.length = sigmaStar.length = 3;

  /* Now rotate the R and P vectors from polar coordinates to Cartesian */
  memcpy( rpolar, tmpVec, 3*sizeof(REAL8));
  memcpy( ppolar, tmpVec+3, 3*sizeof(REAL8));

  rcart[0] = rpolar[0] * cos(rpolar[1]);
  rcart[1] =-rpolar[0] * sin(rpolar[1])*sin(rpolar[2]);
  rcart[2] = rpolar[0] * sin(rpolar[1])*cos(rpolar[2]);

  if( rpolar[1]==0. || rpolar[1]==CST_PI )
  {
    rpolar[1] = CST_PI/2.;

    if( rpolar[1]==0.)
      rpolar[2] = 0.;
    else
      rpolar[2] = CST_PI;

    pcart[0] = ppolar[0]*sin(rpolar[1])*cos(rpolar[2])
					+ ppolar[1]/rpolar[0]*cos(rpolar[1])*cos(rpolar[2])
					- ppolar[2]/rpolar[0]/sin(rpolar[1])*sin(rpolar[2]);
    pcart[1] = ppolar[0]*sin(rpolar[1])*sin(rpolar[2])
					+ ppolar[1]/rpolar[0]*cos(rpolar[1])*sin(rpolar[2])
					+ ppolar[2]/rpolar[0]/sin(rpolar[1])*cos(rpolar[2]);
    pcart[2] = ppolar[0]*cos(rpolar[1])- ppolar[1]/rpolar[0]*sin(rpolar[1]);
  }
  else
  {
    pcart[0] = ppolar[0]*cos(rpolar[1]) -ppolar[1]/rpolar[0]*sin(rpolar[1]);
    pcart[1] =-ppolar[0]*sin(rpolar[1])*sin(rpolar[2])
					-ppolar[1]/rpolar[0]*cos(rpolar[1])*sin(rpolar[2])
					-ppolar[2]/rpolar[0]/sin(rpolar[1])*cos(rpolar[2]);
    pcart[2] = ppolar[0]*sin(rpolar[1])*cos(rpolar[2])
					+ppolar[1]/rpolar[0]*cos(rpolar[1])*cos(rpolar[2])
					-ppolar[2]/rpolar[0]/sin(rpolar[1])*sin(rpolar[2]);
  }

  r.data     = rcart;
  p.data     = pcart;

  spin1.data = tmpVec+6;
  spin2.data = tmpVec+9;
  spin1norm.data = s1normData;
  spin2norm.data = s2normData;
  sigmaKerr.data = sKerrData;
  sigmaStar.data = sStarData;

  memcpy( s1normData, tmpVec+6, 3*sizeof(REAL8) );
  memcpy( s2normData, tmpVec+9, 3*sizeof(REAL8) );

  /* To compute the SigmaKerr and SigmaStar, we need the non-normalized
   * spin values, i.e. S_i. The spins being evolved are S_i/M^2. */
  for ( i = 0; i < 3; i++ )
  {
	 spin1.data[i]  *= mT2;
	 spin2.data[i]  *= mT2;
  }

  /* Calculate various spin parameters */
  EOBCalculateSigmaKerr( &sigmaKerr, &spin1norm, &spin2norm );
  EOBCalculateSigmaStar( &sigmaStar, seobParams->m1,
				seobParams->m2, &spin1norm, &spin2norm );
  a = sqrt( sigmaKerr.data[0]*sigmaKerr.data[0]
			+ sigmaKerr.data[1]*sigmaKerr.data[1]
            + sigmaKerr.data[2]*sigmaKerr.data[2] );
  if ( isnan( a ) )
  {
      PRINT_LOG_INFO(LOG_CRITICAL, "a is nan in GSLSpinPrecHamiltonianWrapperFordHdpphi !!");
      PRINT_LOG_INFO(LOG_CRITICAL, "rpolar, ppolar = %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f", rpolar[0], rpolar[1], rpolar[2], ppolar[0], ppolar[1], ppolar[2]);
      PRINT_LOG_INFO(LOG_CRITICAL, "rcart, pcart = %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f", rcart[0], rcart[1], rcart[2], pcart[0], pcart[1], pcart[2]);
      PRINT_LOG_INFO(LOG_CRITICAL, "a = nan");
      return CEV_FAILURE;;
  }
  REAL8 SpinEOBH = XLALSimIMRSpinPrecEOBHamiltonian( seobParams->eta, &r, &p, &spin1norm, 
    &spin2norm, &sigmaKerr, &sigmaStar, 
    dParams->params->tortoise, dParams->params->seobCoeffs , dParams->params->hParams) / seobParams->eta;

  return SpinEOBH;
}

/**
 * Wrapper for GSL to call the Hamiltonian function. This is simply the function
 * GSLSpinPrecHamiltonianWrapper copied over. The alternative was to make it non-static
 * which increases runtime as static functions can be better optimized.
 */
static double GSLSpinPrecHamiltonianWrapperForRvecDerivs( double x, void *params )
{
    // int debugPK = 1;
    HcapDerivParams *dParams = (HcapDerivParams *)params;

    SpinEOBParams *seobParams = (SpinEOBParams*) dParams->params;
    SpinEOBHCoeffs *coeffs = (SpinEOBHCoeffs*) dParams->params->seobCoeffs;

    REAL8 tmpVec[12]= {0.};
    REAL8 s1normData[3]= {0.}, s2normData[3]= {0.}, sKerrData[3]= {0.}, sStarData[3]= {0.};

    /* These are the vectors which will be used in the call to the Hamiltonian */
    REAL8Vector r, p, spin1, spin2, spin1norm, spin2norm;
    REAL8Vector sigmaKerr, sigmaStar;

    INT4 i;
    REAL8 a;
    REAL8 m1 = seobParams->m1;
    REAL8 m2 = seobParams->m2;
    REAL8 mT2 = (m1+m2)*(m1+m2);
    REAL8 eta = m1*m2/mT2;

    INT4 oldTortoise = dParams->params->tortoise;
    /* Use a temporary vector to avoid corrupting the main function */
    memcpy( tmpVec, dParams->values, sizeof(tmpVec) );

    // if (debugPK){
    //     for( i =0; i < 12; i++)
    //     if( isnan(tmpVec[i]) ) {
    //         XLAL_PRINT_INFO("GSLSpinPrecHamiltonianWrapperForRvecDerivs (from input)::tmpVec %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f\n", tmpVec[0], tmpVec[1], tmpVec[2], tmpVec[3], tmpVec[4], tmpVec[5], tmpVec[6], tmpVec[7], tmpVec[8], tmpVec[9], tmpVec[10], tmpVec[11]);
    //         }
    //     }

    /* Set the relevant entry in the vector to the correct value */
    tmpVec[dParams->varyParam] = x;

    /* Set the LAL-style vectors to point to the appropriate things */
    r.length = p.length = spin1.length = spin2.length = spin1norm.length = spin2norm.length = 3;
    sigmaKerr.length = sigmaStar.length = 3;
    r.data     = tmpVec;
    p.data     = tmpVec+3;

    spin1.data = tmpVec+6;
    spin2.data = tmpVec+9;
    spin1norm.data = s1normData;
    spin2norm.data = s2normData;
    sigmaKerr.data = sKerrData;
    sigmaStar.data = sStarData;

    memcpy( s1normData, tmpVec+6, 3*sizeof(REAL8) );
    memcpy( s2normData, tmpVec+9, 3*sizeof(REAL8) );

    /* To compute the SigmaKerr and SigmaStar, we need the non-normalized
    * spin values, i.e. S_i. The spins being evolved are S_i/M^2. */
    for ( i = 0; i < 3; i++ )
    {
        spin1.data[i]  *= mT2;
        spin2.data[i]  *= mT2;
    }

    /* Calculate various spin parameters */
    XLALSimIMRSpinEOBCalculateSigmaKerr( &sigmaKerr, seobParams->m1,
                    seobParams->m2, &spin1, &spin2 );
    XLALSimIMRSpinEOBCalculateSigmaStar( &sigmaStar, seobParams->m1,
                    seobParams->m2, &spin1, &spin2 );
    a = sqrt( sigmaKerr.data[0]*sigmaKerr.data[0]
                + sigmaKerr.data[1]*sigmaKerr.data[1]
                + sigmaKerr.data[2]*sigmaKerr.data[2] );
    if ( isnan( a ) )
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "a is nan in GSLSpinPrecHamiltonianWrapperForRvecDerivs!!");
        PRINT_LOG_INFO(LOG_CRITICAL, "a = nan");
        return CEV_FAILURE;
    }

    double magR = r.data[0]*r.data[0] + r.data[1]*r.data[1] + r.data[2]*r.data[2];

    // if(debugPK) {
    //     if(0 && magR < 1.96 * 1.96) {
    //     XLAL_PRINT_INFO("GSLSpinPrecHamiltonianWrapperForRvecDerivs (JUST inputs)::tmpVec %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f\n", tmpVec[0], tmpVec[1], tmpVec[2], tmpVec[3], tmpVec[4], tmpVec[5], tmpVec[6], tmpVec[7], tmpVec[8], tmpVec[9], tmpVec[10], tmpVec[11]);

    //     XLAL_PRINT_INFO(" R = %3.10f\n\n", sqrt(magR));
    //     }
    // }

    REAL8 SpinEOBH = XLALSimIMRSpinPrecEOBHamiltonian( seobParams->eta, 
        &r, &p, &spin1norm, &spin2norm, &sigmaKerr, &sigmaStar, 
        dParams->params->tortoise, dParams->params->seobCoeffs, dParams->params->hParams) / seobParams->eta;

    if( isnan(SpinEOBH) )
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "H is nan in GSLSpinPrecHamiltonianWrapperForRvecDerivs.");
        PRINT_LOG_INFO(LOG_CRITICAL, "GSLSpinPrecHamiltonianWrapperForRvecDerivs (JUST inputs)::tmpVec %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f", tmpVec[0], tmpVec[1], tmpVec[2], tmpVec[3], tmpVec[4], tmpVec[5], tmpVec[6], tmpVec[7], tmpVec[8], tmpVec[9], tmpVec[10], tmpVec[11]);
        PRINT_LOG_INFO(LOG_CRITICAL, " R = %3.10f", sqrt(magR));
        PRINT_LOG_INFO(LOG_CRITICAL, "H = nan");
        return CEV_FAILURE;
    }

    if ( dParams->varyParam < 3 )dParams->params->tortoise = oldTortoise;
    return SpinEOBH;
}

/**
 * Function to calculate numerical derivatives of the spin EOB Hamiltonian,
 * to get \f$dr/dt\f$, as decribed in Eqs. A4 of PRD 81, 084041 (2010)
 * This function is not used by the spin-aligned SEOBNRv1 model.
 */
int XLALSpinPrecHcapRvecDerivative(
            double     t,         /**<< UNUSED */
            const  REAL8      values[],  /**<< Dynamical variables */
            REAL8             dvalues[], /**<< Time derivatives of variables (returned) */
            void             *funcParams /**<< EOB parameters */
                               )
{
    // UNUSED int debugPK = 1;
    //if (debugPK){
    // for(int i =0; i < 12; i++){
    //     if( isnan(values[i]) ) {
    //     XLAL_PRINT_INFO("XLALSpinPrecHcapRvecDerivative::values %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f\n", values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9], values[10], values[11]);
    //         XLALPrintError( "XLAL Error - %s: nan in input values \n", __func__);
    //         XLAL_ERROR( XLAL_EINVAL );
    //     }

    //     if( isnan(dvalues[i]) ) {
    //     XLAL_PRINT_INFO("XLALSpinPrecHcapRvecDerivative::dvalues %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f\n", dvalues[0], dvalues[1], dvalues[2], dvalues[3], dvalues[4], dvalues[5], dvalues[6], dvalues[7], dvalues[8], dvalues[9], dvalues[10], dvalues[11]);
    //         XLALPrintError( "XLAL Error - %s: nan in the input dvalues \n", __func__);
    //         XLAL_ERROR( XLAL_EINVAL );
    //     }
    // }
    //}

    static const REAL8 STEP_SIZE = 1.0e-4;

    static const INT4 lMax = 8;

    HcapDerivParams params;

    /* Since we take numerical derivatives wrt dynamical variables */
    /* but we want them wrt time, we use this temporary vector in  */
    /* the conversion */
    REAL8           tmpDValues[14] = {0.};

    REAL8           H; //Hamiltonian
    //REAL8           flux;

    gsl_function F;
    INT4         gslStatus;
    UINT SpinAlignedEOBversion;
    UINT SpinAlignedEOBversionForWaveformCoefficients;

    UINT i, j, k, jj;//, l;

    REAL8Vector rVec, pVec;
    REAL8 rData[3] = {0.}, pData[3] = {0.};

    /* We need r, phi, pr, pPhi to calculate the flux */
    REAL8        r;
    REAL8Vector  polarDynamics;
    REAL8       polData[4] = {0.};

    REAL8 mass1, mass2, eta;
    REAL8 rrTerm2, pDotS1, pDotS2;
    REAL8Vector s1, s2, s1norm, s2norm, sKerr, sStar;
    REAL8       s1Data[3]= {0.}, s2Data[3]= {0.}, s1DataNorm[3]= {0.}, s2DataNorm[3]= {0.};
    REAL8       sKerrData[3]= {0.}, sStarData[3]= {0.};
    REAL8 /*magS1, magS2,*/ chiS, chiA, a, tplspin;
    REAL8 s1dotL, s2dotL;
    REAL8  rcrossrDot[3]= {0.}, rcrossrDotMag, s1dotLN, s2dotLN;


    /* Orbital angular momentum */
    REAL8 Lx, Ly, Lz, magL;
    REAL8 Lhatx, Lhaty, Lhatz;
    //REAL8 dLx, dLy, dLz;
    //REAL8 dLhatx, dLhaty, dMagL;

    //REAL8 alphadotcosi;

    //REAL8 rCrossV_x, rCrossV_y, rCrossV_z, omega;

    /* The error in a derivative as measured by GSL */
    REAL8 absErr;

    REAL8 tmpP[3]= {0.}, rMag, rMag2, prT;
    REAL8 u, u2, u3, u4, u5, w2, a2;
    REAL8 D, m1PlusetaKK, bulk, logTerms, deltaU, deltaT, deltaR;
    REAL8 eobD_r, deltaU_u, deltaU_r, deltaT_r;
    REAL8 dcsi, csi;

    REAL8 tmpValues[12]= {0.};
    REAL8 Tmatrix[3][3]= {{0.}}, invTmatrix[3][3]= {{0.}}, dTijdXk[3][3][3]= {{{0.}}};
    //REAL8 tmpPdotT1[3], tmpPdotT2[3], tmpPdotT3[3]; // 3 terms of Eq. A5

    /* Set up pointers for GSL */
    params.values  = values;
    params.params  = (SpinEOBParams *)funcParams;

    F.function = &GSLSpinPrecHamiltonianWrapperForRvecDerivs;
    F.params   = &params;

    mass1 = params.params->m1;
    mass2 = params.params->m2;
    eta   = params.params->eta;
    SpinAlignedEOBversion = 4;
    SpinEOBHCoeffs *coeffs = (SpinEOBHCoeffs*) params.params->seobCoeffs;

    /* For precessing binaries, the effective spin of the Kerr
    * background evolves with time. The coefficients used to compute
    * the Hamiltonian depend on the Kerr spin, and hence need to
    * be updated for the current spin values */
#if 0
    if ( 0 )
    {
        /*{{{*/
        /* Set up structures and calculate necessary (spin-only) PN parameters */
        /* Due to precession, these need to get calculated in every step */
        //memset( params.params->seobCoeffs, 0, sizeof(SpinEOBHCoeffs) );

        REAL8 tmps1Data[3]= {0.}, tmps2Data[3]= {0.}; 
        REAL8Vector tmps1Vec, tmps2Vec;
        memcpy( tmps1Data, values+6, 3*sizeof(REAL8) );
        memcpy( tmps2Data, values+9, 3*sizeof(REAL8) );
        tmps1Vec.data   = tmps1Data; tmps2Vec.data   = tmps2Data;
        tmps1Vec.length = tmps2Vec.length = 3;

        REAL8Vector *tmpsigmaKerr = NULL;
        REAL8Vector *tmpsigmaStar = NULL;
        if ( !(tmpsigmaKerr = CreateREAL8Vector( 3 )) )
        {
            return CEV_FAILURE;
        }

        if ( !(tmpsigmaStar = CreateREAL8Vector( 3 )) )
        {
            return CEV_FAILURE;
        }

        if ( XLALSimIMRSpinEOBCalculateSigmaKerr( tmpsigmaKerr, mass1, mass2,
                                    &tmps1Vec, &tmps2Vec ) == CEV_FAILURE )
        {
            DestroyREAL8Vector( tmpsigmaKerr );
            return CEV_FAILURE;
        }

        if ( XLALSimIMRSpinEOBCalculateSigmaStar( tmpsigmaStar, mass1, mass2,
                                    &tmps1Vec, &tmps2Vec ) == CEV_FAILURE )
        {
            DestroyREAL8Vector( tmpsigmaKerr );
            DestroyREAL8Vector( tmpsigmaStar );
            return CEV_FAILURE;
        }

        /* Update a with the Kerr background spin
            * Pre-compute the Hamiltonian coefficients            */
        //REAL8Vector *delsigmaKerr 	= params.params->sigmaKerr;
        params.params->sigmaKerr 	= tmpsigmaKerr;
        params.params->sigmaStar 	= tmpsigmaStar;
        params.params->a 		= sqrt( tmpsigmaKerr->data[0]*tmpsigmaKerr->data[0]
                    + tmpsigmaKerr->data[1]*tmpsigmaKerr->data[1]
                    + tmpsigmaKerr->data[2]*tmpsigmaKerr->data[2] );
        //tmpsigmaKerr->data[2];
        if ( XLALSimIMRCalculateSpinPrecEOBHCoeffs( params.params->seobCoeffs, eta,
                params.params->a, SpinAlignedEOBversion ) == CEV_FAILURE )
        {
            DestroyREAL8Vector( tmpsigmaKerr );
            DestroyREAL8Vector( tmpsigmaStar );
            return CEV_FAILURE;
        }

        // params.params->seobCoeffs->SpinAlignedEOBversion = SpinAlignedEOBversion;
        /* Release the old memory */
        //if(0)DestroyREAL8Vector( delsigmaKerr );
        /*}}}*/}
#endif
        /* Set the position/momenta vectors to point to the appropriate things */
        rVec.length = pVec.length = 3;
        rVec.data   = rData;
        pVec.data   = pData;
        memcpy( rData, values, sizeof(rData) );
        memcpy( pData, values+3, sizeof(pData) );

        /* We need to re-calculate the parameters at each step as precessing
        * spins will not be constant */

        /* We cannot point to the values vector directly as it leads to a warning */
        s1.length = s2.length = s1norm.length = s2norm.length = 3;
        s1.data = s1Data;
        s2.data = s2Data;
        s1norm.data = s1DataNorm;
        s2norm.data = s2DataNorm;

        memcpy( s1Data, values+6, 3*sizeof(REAL8) );
        memcpy( s2Data, values+9, 3*sizeof(REAL8) );
        memcpy( s1DataNorm, values+6, 3*sizeof(REAL8) );
        memcpy( s2DataNorm, values+9, 3*sizeof(REAL8) );

        for ( i = 0; i < 3; i++ )
        {
            s1Data[i] *= (mass1+mass2)*(mass1+mass2);
            s2Data[i] *= (mass1+mass2)*(mass1+mass2);
        }

        sKerr.length = 3;
        sKerr.data   = sKerrData;
        XLALSimIMRSpinEOBCalculateSigmaKerr( &sKerr, mass1, mass2, &s1, &s2 );

        sStar.length = 3;
        sStar.data   = sStarData;
        XLALSimIMRSpinEOBCalculateSigmaStar( &sStar, mass1, mass2, &s1, &s2 );

        a = sqrt(sKerr.data[0]*sKerr.data[0] + sKerr.data[1]*sKerr.data[1]
            + sKerr.data[2]*sKerr.data[2]);

        if (isnan(a))
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "a = nan");
            return CEV_FAILURE;
        }
        // if(debugPK && isnan(a))
        // XLAL_PRINT_INFO("a is nan in XLALSpinPrecHcapRvecDerivative \n");

        ///* set the tortoise flag to 2 */
        //INT4 oldTortoise = params.params->tortoise;
        //params.params->tortoise = 2;

        /* Convert momenta to p */
        rMag = sqrt(rData[0]*rData[0] + rData[1]*rData[1] + rData[2]*rData[2]);
        prT = pData[0]*(rData[0]/rMag) + pData[1]*(rData[1]/rMag)
                        + pData[2]*(rData[2]/rMag);

        rMag2 = rMag * rMag;
        u  = 1./rMag;
        u2 = u*u;
        u3 = u2*u;
        u4 = u2*u2;
        u5 = u4*u;
        a2 = a*a;
        w2 = rMag2 + a2;
        /* Eq. 5.83 of BB1, inverse */
        D = 1. + log(1. + 6.*eta*u2 + 2.*(26. - 3.*eta)*eta*u3);
        eobD_r =  (u2/(D*D))*(12.*eta*u + 6.*(26. - 3.*eta)*eta*u2)/(1.
            + 6.*eta*u2 + 2.*(26. - 3.*eta)*eta*u3);
        m1PlusetaKK = -1. + eta * coeffs->KK;
        /* Eq. 5.75 of BB1 */
        bulk = 1./(m1PlusetaKK*m1PlusetaKK) + (2.*u)/m1PlusetaKK + a2*u2;
        /* Eq. 5.73 of BB1 */
        logTerms = 1. + eta*coeffs->k0 + eta*log(fabs(1. + coeffs->k1*u
        + coeffs->k2*u2 + coeffs->k3*u3 + coeffs->k4*u4
        + coeffs->k5*u5 + coeffs->k5l*u5*log(u)));
        /* Eq. 5.73 of BB1 */
        deltaU = bulk*logTerms;
        deltaU = fabs(deltaU);

        /* Eq. 5.71 of BB1 */
        deltaT = rMag2*deltaU;
        /* ddeltaU/du */
        deltaU_u = 2.*(1./m1PlusetaKK + a2*u)*logTerms +
        bulk * (eta*(coeffs->k1 + u*(2.*coeffs->k2 + u*(3.*coeffs->k3
        + u*(4.*coeffs->k4 + 5.*(coeffs->k5+coeffs->k5l*log(u))*u)))))
        / (1. + coeffs->k1*u + coeffs->k2*u2 + coeffs->k3*u3
        + coeffs->k4*u4 + (coeffs->k5+coeffs->k5l*log(u))*u5);
        deltaU_r = -u2 * deltaU_u;
        /* Eq. 5.38 of BB1 */
        deltaR = deltaT*D;
        if ( params.params->tortoise )
            csi = sqrt( fabs(deltaT * deltaR) )/ w2; /* Eq. 28 of Pan et al. PRD 81, 084041 (2010) */
        else
            csi = 1.0;

        for( i = 0; i < 3; i++ )
        {
            tmpP[i] = pData[i] - (rData[i]/rMag) * prT * (csi-1.)/csi;
        }


        /* Calculate the T-matrix, required to convert P from tortoise to
        * non-tortoise coordinates, and/or vice-versa. This is given explicitly
        * in Eq. A3 of 0912.3466 */
        for( i = 0; i < 3; i++ )
            for( j = 0; j <= i; j++ )
            {
                Tmatrix[i][j] = Tmatrix[j][i] = (rData[i]*rData[j]/rMag2)
                                            * (csi - 1.);

                invTmatrix[i][j] = invTmatrix[j][i] =
                        - (csi - 1)/csi * (rData[i]*rData[j]/rMag2);

                if( i==j ){
                    Tmatrix[i][j]++;
                    invTmatrix[i][j]++;  }

            }

        dcsi = csi * (2./rMag + deltaU_r/deltaU) + csi*csi*csi
            / (2.*rMag2*rMag2 * deltaU*deltaU) * ( rMag*(-4.*w2)/D - eobD_r*(w2*w2));

        for( i = 0; i < 3; i++ )
            for( j = 0; j < 3; j++ )
                for( k = 0; k < 3; k++ )
                {
                    dTijdXk[i][j][k]  =
                        (rData[i]*KRONECKER_DELTA(j,k) + KRONECKER_DELTA(i,k)*rData[j])
                        *(csi - 1.)/rMag2
                        + rData[i]*rData[j]*rData[k]/rMag2/rMag*(-2./rMag*(csi - 1.) + dcsi);
                }

        //if (debugPK){
        for(i =0; i < 12; i++)
        {
            if( isnan(values[i]) ) 
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "values %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f", 
                    values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9], values[10], values[11]);
                PRINT_LOG_INFO(LOG_CRITICAL, "values = nan");
                return CEV_FAILURE;
            }

            if( isnan(dvalues[i]) ) 
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "dvalues %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f\n", dvalues[0], dvalues[1], dvalues[2], dvalues[3], dvalues[4], dvalues[5], dvalues[6], dvalues[7], dvalues[8], dvalues[9], dvalues[10], dvalues[11]);
                PRINT_LOG_INFO(LOG_CRITICAL, "dvalues = nan");
                return CEV_FAILURE;
            }
        }
    //}

    /* Now calculate derivatives w.r.t. each parameter */
    // RH: Check if components i=0..5 (position and momenta) do not change the a parameter used for spin
    // RH: this breaks the loop below for i>=6
    SpinEOBHCoeffs tmpCoeffs;
    {
        // RH: taken from GSLSpinHamiltonianWrapperForRvecDerivs
        /* These are the vectors which will be used in the call to the Hamiltonian */
        REAL8Vector spin1, spin2;
        REAL8Vector sigmaKerr;
        REAL8 tmpVec[12]= {0.};
        REAL8 tmpsKerrData[3]= {0.};
        REAL8 mT2 = (mass1+mass2)*(mass1+mass2);

        /* Use a temporary vector to avoid corrupting the main function */
        memcpy( tmpVec, values, sizeof(tmpVec) );

        /* Set the LAL-style vectors to point to the appropriate things */
        sigmaKerr.length = 3;
        spin1.length = 3;
        spin2.length = 3;

        spin1.data = tmpVec+6;
        spin2.data = tmpVec+9;
        sigmaKerr.data = tmpsKerrData;

        /* To compute the SigmaKerr and SigmaStar, we need the non-normalized
            * spin values, i.e. S_i. The spins being evolved are S_i/M^2. */
        for ( i = 0; i < 3; i++ )
        {
            spin1.data[i]  *= mT2;
            spin2.data[i]  *= mT2;
        }

        /* Calculate various spin parameters */
        XLALSimIMRSpinEOBCalculateSigmaKerr( &sigmaKerr, mass1, mass2,
                                                &spin1, &spin2 );

        REAL8 tmpa;
        /* Calculate the orbital angular momentum */
        Lx = values[1] * values[5] - values[2] * values[4];
        Ly = values[2] * values[3] - values[0] * values[5];
        Lz = values[0] * values[4] - values[1] * values[3];

        magL = sqrt(Lx * Lx + Ly * Ly + Lz * Lz);
        Lhatx = Lx / magL;
        Lhaty = Ly / magL;
        Lhatz = Lz / magL;
        REAL8 Lhat[3] = {Lhatx, Lhaty, Lhatz};
        REAL8 tempS1_p = inner_product3d(s1Data, Lhat);
        REAL8 tempS2_p = inner_product3d(s2Data, Lhat);
        REAL8 S1_perp[3] = {0, 0, 0};
        REAL8 S2_perp[3] = {0, 0, 0};
        for(jj=0; jj<3; jj++)
        {
            S1_perp[jj] = 1/mT2*(s1Data[jj]-tempS1_p*Lhat[jj]);
            S2_perp[jj] = 1/mT2*(s2Data[jj]-tempS2_p*Lhat[jj]);
        }
        REAL8 sKerr_norm = sqrt(inner_product3d(sigmaKerr.data, sigmaKerr.data));
        REAL8 S_con = 0.0;
        if (sKerr_norm>1e-6)
        {
            S_con = sigmaKerr.data[0] * Lhat[0] + sigmaKerr.data[1] * Lhat[1] + sigmaKerr.data[2] * Lhat[2];
            S_con /= (1 - 2 * eta);
            S_con += (inner_product3d(S1_perp, sigmaKerr.data) + inner_product3d(S2_perp, sigmaKerr.data)) / sKerr_norm / (1 - 2 * eta) / 2.;
        }

        REAL8 chi = S_con;
        tmpa = sqrt(sigmaKerr.data[0]*sigmaKerr.data[0]
                    + sigmaKerr.data[1]*sigmaKerr.data[1]
                    + sigmaKerr.data[2]*sigmaKerr.data[2]);
        // if (SpinAlignedEOBversion == 4)
        // {
        if ( XLALSimIMRCalculateSpinPrecEOBHCoeffs_v2( &tmpCoeffs, eta,
            tmpa, chi, SpinAlignedEOBversion , params.params->hParams) == CEV_FAILURE )
        {
            return CEV_FAILURE;
        }
        // }
        // else
        // {
        //     if ( XLALSimIMRCalculateSpinPrecEOBHCoeffs( &tmpCoeffs, eta,
        //         tmpa, coeffs->SpinAlignedEOBversion ) == CEV_FAILURE )
        //     {
        //     return CEV_FAILURE;
        //     }
        // }


        // tmpCoeffs.SpinAlignedEOBversion = SpinAlignedEOBversion;
        tmpCoeffs.updateHCoeffs = 0;
    }
    SpinEOBHCoeffs *oldCoeffs = params.params->seobCoeffs;
    params.params->seobCoeffs = &tmpCoeffs;
    for ( i = 0; i < 6; i++ )
    {
        params.varyParam = i;
        if ( i >=6 && i < 9 )
        {
            return CEV_FAILURE; // this should never happen
            params.params->seobCoeffs->updateHCoeffs = 1;
            gslStatus = gsl_deriv_central( &F, values[i],
                            STEP_SIZE*mass1*mass1, &tmpDValues[i], &absErr );
        }
        else if ( i >= 9 )
        {
            return CEV_FAILURE; // this should never happen
            params.params->seobCoeffs->updateHCoeffs = 1;
            gslStatus = gsl_deriv_central( &F, values[i],
                            STEP_SIZE*mass2*mass2, &tmpDValues[i], &absErr );
        }
        else if ( i < 3 )
        {
            // return CEV_FAILURE; // this should never happen
            params.params->tortoise = 2;
            memcpy( tmpValues, params.values, sizeof(tmpValues) );
            tmpValues[3] = tmpP[0]; tmpValues[4] = tmpP[1]; tmpValues[5] = tmpP[2];
            params.values = tmpValues;

            gslStatus = gsl_deriv_central( &F, values[i],
                            STEP_SIZE, &tmpDValues[i], &absErr );

            params.values = values;
            params.params->tortoise = 1;
        }
        else
        {
            params.params->seobCoeffs->updateHCoeffs = 1;

            gslStatus = gsl_deriv_central( &F, values[i],
                            STEP_SIZE, &tmpDValues[i], &absErr );
        }
        if ( gslStatus != GSL_SUCCESS )
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "Failure in GSL function");
            return CEV_FAILURE;
        }
    }
    params.params->seobCoeffs = oldCoeffs;
    // if (debugPK){
    // for( i =0; i < 12; i++)
    //     if( isnan(tmpDValues[i]) ) 
    //     {
    //         XLAL_PRINT_INFO("XLALSpinPrecHcapRvecDerivative (just after diff)::tmpDValues %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f\n", tmpDValues[0], tmpDValues[1], tmpDValues[2], tmpDValues[3], tmpDValues[4], tmpDValues[5], tmpDValues[6], tmpDValues[7], tmpDValues[8], tmpDValues[9], tmpDValues[10], tmpDValues[11]);
    //     }
    // }

    /* Calculate the orbital angular momentum */
    Lx = values[1]*values[5] - values[2]*values[4];
    Ly = values[2]*values[3] - values[0]*values[5];
    Lz = values[0]*values[4] - values[1]*values[3];

    magL = sqrt( Lx*Lx + Ly*Ly + Lz*Lz );

    Lhatx = Lx/magL;
    Lhaty = Ly/magL;
    Lhatz = Lz/magL;

    /* Calculate the polar data */
    polarDynamics.length = 4;
    polarDynamics.data   = polData;

    r = polData[0] = sqrt( values[0]*values[0] + values[1]*values[1]
                        + values[2]*values[2] );
    polData[1] = 0;
    polData[2] = (values[0]*values[3] + values[1]*values[4]
                + values[2]*values[5]) / polData[0];
    polData[3] = magL;


    // Compute \vec{S_i} \dot \vec{L}
    s1dotL = (s1Data[0]*Lhatx + s1Data[1]*Lhaty + s1Data[2]*Lhatz)
            / (mass1*mass1);
    s2dotL = (s2Data[0]*Lhatx + s2Data[1]*Lhaty + s2Data[2]*Lhatz)
            / (mass2*mass2);

    // Compute \vec{L_N} = \vec{r} \times \.{\vec{r}},
    // \vec{S_i} \dot \vec{L_N} and chiS and chiA
    rcrossrDot[0] = values[1]*tmpDValues[5] - values[2]*tmpDValues[4];
    rcrossrDot[1] = values[2]*tmpDValues[3] - values[0]*tmpDValues[5];
    rcrossrDot[2] = values[0]*tmpDValues[4] - values[1]*tmpDValues[3];
    rcrossrDotMag = sqrt( rcrossrDot[0]*rcrossrDot[0]
        + rcrossrDot[1]*rcrossrDot[1]	+ rcrossrDot[2]*rcrossrDot[2] );

    rcrossrDot[0] /= rcrossrDotMag;
    rcrossrDot[1] /= rcrossrDotMag;
    rcrossrDot[2] /= rcrossrDotMag;

    s1dotLN = (s1Data[0]*rcrossrDot[0] + s1Data[1]*rcrossrDot[1]
                + s1Data[2]*rcrossrDot[2]) / (mass1*mass1);
    s2dotLN = (s2Data[0]*rcrossrDot[0] + s2Data[1]*rcrossrDot[1]
            + s2Data[2]*rcrossrDot[2]) / (mass2*mass2);

    REAL8 cax, cay, caz, csx, csy, csz;
    REAL8 c1x, c1y, c1z, c2x, c2y, c2z;
    REAL8 m1sq = mass1*mass1;
    REAL8 m2sq = mass2*mass2;
    REAL8 magY;
    REAL8 yhat[3] = {0,0,0};
    // Calculate chiAVec, chiSVec
    yhat[0] = rcrossrDot[1]*values[2] - rcrossrDot[2]*values[1];
    yhat[1] = rcrossrDot[2]*values[0] - rcrossrDot[0]*values[2];
    yhat[2] = rcrossrDot[0]*values[1] - rcrossrDot[1]*values[0];
    magY = sqrt(yhat[0]*yhat[0] + yhat[1]*yhat[1] + yhat[2]*yhat[2]);

    c1x = (s1.data[0] * values[0] + s1.data[1] * values[1] + s1.data[2] * values[2]) / (r * m1sq);
    c1z = (s1.data[0] * rcrossrDot[0] + s1.data[1] * rcrossrDot[1] + s1.data[2] * rcrossrDot[2]) / (rcrossrDotMag * m1sq);
    // c1y = sqrt( (s1.data[0]*s1.data[0] + s1.data[1]*s1.data[1] + s1.data[2]*s1.data[2])/(m1sq*m1sq) - c1x*c1x - c1z*c1z);
    c1y = (s1.data[0] * yhat[0] + s1.data[1] * yhat[1] + s1.data[2] * yhat[2]) / m1sq / magY;

    c2x = (s2.data[0] * values[0] + s2.data[1] * values[1] + s2.data[2] * values[2]) / (r * m2sq);
    c2z = (s2.data[0] * rcrossrDot[0] + s2.data[1] * rcrossrDot[1] + s2.data[2] * rcrossrDot[2]) / (rcrossrDotMag * m2sq);
    // c2y = sqrt( (s2.data[0]*s2.data[0] + s2.data[1]*s2.data[1] + s2.data[2]*s2.data[2])/(m2sq*m2sq) - c2x*c2x - c2z*c2z);
    c2y = (s2.data[0] * yhat[0] + s2.data[1] * yhat[1] + s2.data[2] * yhat[2]) / m2sq / magY;

    cax = 0.5 * (c1x - c2x);
    cay = 0.5 * (c1y - c2y);
    caz = 0.5 * (c1z - c2z);

    csx = 0.5 * (c1x + c2x);
    csy = 0.5 * (c1y + c2y);
    csz = 0.5 * (c1z + c2z);

    REAL8 mT2 = (mass1+mass2)*(mass1+mass2);
    if (SpinAlignedEOBversion==4)
    {
        chiS = 0.5 * (s1dotL + s2dotL);
        chiA = 0.5 * (s1dotL - s2dotL);
    }
    else
    {
        chiS = 0.5 * (s1dotLN + s2dotLN);
        chiA = 0.5 * (s1dotLN - s2dotLN);
    }
    REAL8 Lhat[3] = {Lhatx, Lhaty, Lhatz};
    REAL8 tempS1_p = inner_product3d(s1Data, Lhat);
    REAL8 tempS2_p = inner_product3d(s2Data, Lhat);
    REAL8 S1_perp[3] = {0, 0, 0};
    REAL8 S2_perp[3] = {0, 0, 0};

    for ( jj = 0; jj < 3; jj++)
    {
        S1_perp[jj] = 1 / mT2 * (s1Data[jj] - tempS1_p * Lhat[jj]);
        S2_perp[jj] = 1 / mT2 * (s2Data[jj] - tempS2_p * Lhat[jj]);
    }

    REAL8 sKerr_norm = sqrt(inner_product3d(sKerr.data, sKerr.data));

    REAL8 S_con = 0.0;
    if (sKerr_norm > 1e-6)
    {
        S_con = sKerr.data[0] * Lhat[0] + sKerr.data[1] * Lhat[1] + sKerr.data[2] * Lhat[2];
        S_con /= (1 - 2 * eta);
        S_con += (inner_product3d(S1_perp, sKerr.data) + inner_product3d(S2_perp, sKerr.data)) / sKerr_norm / (1 - 2 * eta) / 2.;
    }

    REAL8 chi = S_con;
    /* Compute the test-particle limit spin of the deformed-Kerr background */
    switch ( SpinAlignedEOBversion )
    {
        case 1:
            tplspin = 0.0;
            break;
        case 2:
        case 4:
            tplspin = (1.-2.*eta) * chiS + (mass1 - mass2)/(mass1 + mass2) * chiA;
            break;
        default:
            PRINT_LOG_INFO(LOG_CRITICAL, "Unknown SEOBNR version! At present only v1 and v2 are available.");
            return CEV_FAILURE;
            break;
    }

    for( i = 0; i< 3; i++ )
    {
        params.params->s1Vec->data[i]     = s1norm.data[i];
        params.params->s2Vec->data[i]     = s2norm.data[i];
        params.params->sigmaStar->data[i] = sStar.data[i];
        params.params->sigmaKerr->data[i] = sKerr.data[i];
    }

    //params.params->s1Vec     = &s1norm;
    //params.params->s2Vec     = &s2norm;
    //params.params->sigmaStar = &sStar;
    //params.params->sigmaKerr = &sKerr;
    params.params->a = a;

    if (SpinAlignedEOBversion == 2) 
    {
        SpinAlignedEOBversionForWaveformCoefficients = 3;
    }
    else
    {
        SpinAlignedEOBversionForWaveformCoefficients = SpinAlignedEOBversion;
    }

    if (CODE_VERSION == 3)
    {
        // if (params.params->alignedSpins==1) 
        //     EccPrec_CalcSpinPrecFacWaveformCoefficients(params.params->hCoeffs, mass1, mass2, eta, tplspin,
        //                     chiS, chiA, csx, csy, csz, cax, cay, caz, 451);
        // else
            EccPrec_CalcSpinPrecFacWaveformCoefficients(params.params->hCoeffs, mass1, mass2, eta, tplspin,
                            chiS, chiA, csx, csy, csz, cax, cay, caz, 451);
    }
    else
    {
        if (params.params->alignedSpins==1) {
            XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients(
                                                            params.params->hCoeffs, mass1, mass2, eta, tplspin,
                                                                                                            chiS, chiA, SpinAlignedEOBversion);
                                                                                                        }
        else {
            XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients(
                                                            params.params->hCoeffs, mass1, mass2, eta, tplspin,
                                                            chiS, chiA, SpinAlignedEOBversionForWaveformCoefficients);
                                                        }
    }
    // if (SpinAlignedEOBversion == 4){

    XLALSimIMRCalculateSpinPrecEOBHCoeffs_v2( params.params->seobCoeffs, eta, a, chi,
        SpinAlignedEOBversion , params.params->hParams);

    //H = XLALSimIMRSpinPrecEOBHamiltonian( eta, &rVec, &pVec, s1proj, s2proj,
    //&sKerr, &sStar, params.params->tortoise, params.params->seobCoeffs );
    // }
    // else{
    // XLALSimIMRCalculateSpinPrecEOBHCoeffs( params.params->seobCoeffs, eta, a,
    //     SpinAlignedEOBversion );
    // }
    H = XLALSimIMRSpinPrecEOBHamiltonian( eta, &rVec, &pVec, &s1norm, &s2norm,
                    &sKerr, &sStar, params.params->tortoise, params.params->seobCoeffs , params.params->hParams);
    H = H * (mass1 + mass2);

    /* Now make the conversion */
    /* rVectorDot */
    for( i = 0; i < 3; i++ )
        for( j = 0, dvalues[i] = 0.; j < 3; j++ )
            dvalues[i] += tmpDValues[j+3]*Tmatrix[i][j];
    
    // DEBUG: non-tortoise pVecDot
    // print_debug("%e\t%e\t%e\n", tmpDValues[0], tmpDValues[1], tmpDValues[2]);
    dvalues[3] = -tmpDValues[0];
    dvalues[4] = -tmpDValues[1];
    dvalues[5] = -tmpDValues[2];
    return CEV_SUCCESS;
}


REAL8 XLALSimIMRSpinPrecEOBCalcOmega(
                      const REAL8           values[],   /**<< Dynamical variables */
                      SpinEOBParams         *funcParams /**<< EOB parameters */
                      )
{
    // int debugPK = 1;
    // if (debugPK){
    //     for(int i =0; i < 12; i++)
    //     if( isnan(values[i]) ) {
    //         XLAL_PRINT_INFO("XLALSimIMRSpinPrecEOBCalcOmega::values %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f\n", values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9], values[10], values[11]);
    //         XLALPrintError( "XLAL Error - %s: nan in input values  \n", __func__);
    //         XLAL_ERROR( XLAL_EINVAL );
    //     }
    // }

    /* ********************************************************************* */
    /* ************ Memory Allocation ************************************** */
    /* ********************************************************************* */
    static const REAL8 STEP_SIZE = 1.0e-4;
    REAL8 tmpvar = 0;

    HcapDerivParams params;

    /* Cartesian values for calculating the Hamiltonian */
    REAL8 cartValues[14] = {0.}, dvalues[14] = {0.};
    REAL8 cartvalues[14] = {0.}, polarvalues[6] = {0.}; /* The rotated cartesian/polar values */
    REAL8 polarRPcartSvalues[14] = {0.};
    memcpy( cartValues, values, 14 * sizeof(REAL8) );

    INT4 i, j;

    REAL8 rvec[3]  = {0.,0,0}, pvec[3]  = {0.,0,0};
    REAL8 s1vec[3] = {0.,0,0}, s2vec[3] = {0.,0,0};

    REAL8 rdotvec[3] = {0.,0,0};
    REAL8 rvecprime[3] = {0.,0,0}, pvecprime[3] = {0.,0,0},
            s1vecprime[3]= {0.,0,0}, s2vecprime[3]= {0.,0,0};
    REAL8 rvectmp[3]   = {0.,0,0}, pvectmp[3] = {0.,0,0},
            s1vectmp[3]  = {0.,0,0}, s2vectmp[3]= {0.,0,0};
    REAL8 LNhatprime[3]= {0.,0,0}, LNhatTmp[3]= {0.,0,0};
    REAL8 rcrossrdot[3] = {0.,0,0};

    REAL8 Rot1[3][3] ={{0.}}; // Rotation matrix for prevention of blowing up
    REAL8 Rot2[3][3] ={{0.}} ;
    REAL8 LNhat[3] = {0.,0,0};

    REAL8        Xhat[3] = {1, 0, 0};
    REAL8 Yhat[3] = {0, 1, 0};
    REAL8 Zhat[3] = {0, 0, 1};

    REAL8 Xprime[3] = {0.,0,0}, Yprime[3] = {0.,0,0}, Zprime[3] = {0.,0,0};

    gsl_function F;
    INT4         gslStatus;

    REAL8 omega;

    /* The error in a derivative as measured by GSL */
    REAL8 absErr;

    /* ********************************************************************* */
    /* ************ Main Logic begins ************************************ */
    /* ********************************************************************* */

    /* Copy over the coordinates and spins */
    memcpy( rvec,  values,   3*sizeof(REAL8) );
    memcpy( pvec,  values+3, 3*sizeof(REAL8) );
    memcpy( s1vec, values+6, 3*sizeof(REAL8) );
    memcpy( s2vec, values+9, 3*sizeof(REAL8) );

    /* Calculate rDot = \f$\partial Hreal / \partial p_r\f$ */
    memset( dvalues, 0, 14 * sizeof(REAL8) );

    if( XLALSpinPrecHcapRvecDerivative( 0, values, dvalues,
                                    (void*) funcParams) == CEV_FAILURE )
    {
        return REAL8_FAIL_NAN;
    }
    memcpy( rdotvec, dvalues, 3*sizeof(REAL8) );
    // if (debugPK){
    //     for(int ii =0; ii < 12; ii++)
    //     if( isnan(dvalues[ii]) ) {
    //         XLAL_PRINT_INFO("XLALSimIMRSpinPrecEOBCalcOmega::dvalues %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f\n", dvalues[0], dvalues[1], dvalues[2], dvalues[3], dvalues[4], dvalues[5], dvalues[6], dvalues[7], dvalues[8], dvalues[9], dvalues[10], dvalues[11]);
    //         XLALPrintError( "XLAL Error - %s: nan in dvalues \n", __func__);
    //         XLAL_ERROR( XLAL_EINVAL );
    //     }
    // }

    /* Calculate LN = r cross rDot */
    cross_product3d( rvec, rdotvec, rcrossrdot );
    REAL8 rcrossrdotNorm = sqrt(inner_product3d( rcrossrdot, rcrossrdot ));
    for( i = 0; i < 3; i++ )
        rcrossrdot[i] /= rcrossrdotNorm;
    memcpy( LNhat, rcrossrdot, 3 * sizeof(REAL8) );


    /* ********************************************************************* */
    /* First, the frame is rotated so that L is along the y-axis. */
    /* this rotation includes the spins. */
    /* ********************************************************************* */

    // For Now , set first rotation matrix to identity
    // Check if LNhat and Xhat are too aligned, in which case rotate LNhat
    if( inner_product3d(LNhat, Xhat) < 0.9 )
    {
        Rot1[0][0] = 1.; Rot1[0][1] = 0; Rot1[0][2] = 0;
        Rot1[1][0] = 0.; Rot1[1][1] = 1; Rot1[1][2] = 0;
        Rot1[2][0] = 0.; Rot1[2][1] = 0; Rot1[2][2] = 1;

        memcpy(Xprime, LNhat, 3 * sizeof(REAL8));
        cross_product3d( Xprime, Xhat, Yprime );
        tmpvar = sqrt(inner_product3d(Yprime, Yprime));

        for( i=0; i<3; i++)
        Yprime[i] /= tmpvar;

        cross_product3d(Xprime, Yprime, Zprime);
        tmpvar = sqrt(inner_product3d(Zprime, Zprime));
        for( i=0; i<3; i++)
        Zprime[i] /= tmpvar;
    }
    else
    {
        Rot1[0][0] = 1./sqrt(2); Rot1[0][1] = -1/sqrt(2); Rot1[0][2] = 0;
        Rot1[1][0] = 1./sqrt(2); Rot1[1][1] = 1./sqrt(2); Rot1[1][2] = 0;
        Rot1[2][0] = 0.;         Rot1[2][1] = 0;          Rot1[2][2] = 1;
        LNhatTmp[0] = LNhatTmp[1] = LNhatTmp[2] = 0.;

        for(i=0; i<3; i++)
        for(j=0; j<3; j++)
            LNhatTmp[i] += Rot1[i][j]*LNhat[j];

        memcpy(Xprime, LNhatTmp, 3*sizeof(REAL8));
        cross_product3d(Xprime, Xhat, Yprime);
        tmpvar = sqrt(inner_product3d(Yprime, Yprime));

        for( i=0; i<3; i++)
        Yprime[i] /= tmpvar;

        cross_product3d(Xprime, Yprime, Zprime);
        tmpvar = sqrt(inner_product3d(Zprime, Zprime));
        for( i=0; i<3; i++)
        Zprime[i] /= tmpvar;
    }

    Rot2[0][0] = Xprime[0]; Rot2[0][1] = Xprime[1]; Rot2[0][2] = Xprime[2];
    Rot2[1][0] = Yprime[0]; Rot2[1][1] = Yprime[1]; Rot2[1][2] = Yprime[2];
    Rot2[2][0] = Zprime[0]; Rot2[2][1] = Zprime[1]; Rot2[2][2] = Zprime[2];

    memset( rvectmp,    0, 3 * sizeof(REAL8) );
    memset( pvectmp,    0, 3 * sizeof(REAL8) );
    memset( s1vectmp,   0, 3 * sizeof(REAL8) );
    memset( s2vectmp,   0, 3 * sizeof(REAL8) );
    memset( rvecprime,  0, 3 * sizeof(REAL8) );
    memset( pvecprime,  0, 3 * sizeof(REAL8) );
    memset( s1vecprime, 0, 3 * sizeof(REAL8) );
    memset( s2vecprime, 0, 3 * sizeof(REAL8) );
    memset( LNhatprime, 0, 3 * sizeof(REAL8) );
    memset( LNhatTmp,   0, 3 * sizeof(REAL8) );

    /* Perform the actual rotation */
    for (i=0; i<3; i++)
        for(j=0; j<3; j++)
        {
            rvectmp[i]  += Rot1[i][j]*rvec[j];
            pvectmp[i]  += Rot1[i][j]*pvec[j];
            s1vectmp[i] += Rot1[i][j]*s1vec[j];
            s2vectmp[i] += Rot1[i][j]*s2vec[j];
            LNhatTmp[i] += Rot1[i][j]*LNhat[j];
        }
    for (i=0; i<3; i++)
        for(j=0; j<3; j++)
        {
            rvecprime[i]  += Rot2[i][j]*rvectmp[j];
            pvecprime[i]  += Rot2[i][j]*pvectmp[j];
            s1vecprime[i] += Rot2[i][j]*s1vectmp[j];
            s2vecprime[i] += Rot2[i][j]*s2vectmp[j];
            LNhatprime[i] += Rot2[i][j]*LNhatTmp[j];
        }

    memcpy(cartvalues,   rvecprime,  3*sizeof(REAL8));
    memcpy(cartvalues+3, pvecprime,  3*sizeof(REAL8));
    memcpy(cartvalues+6, s1vecprime, 3*sizeof(REAL8));
    memcpy(cartvalues+9, s2vecprime, 3*sizeof(REAL8));

    /* ********************************************************************* */
    /* Second, \f$\vec{r}\f$ and \f$\vec{p}\f$ are converted to polar
    * coordinates (and not the spins).
    * As L is along the y-axis, \f$\theta\f$ defined as the angle between
    * L and the y-axis is 0, which is a cyclic coordinate now and that fixes
    * \f$p_\theta = 0\f$. */
    /* ********************************************************************* */

    /** the polarvalues, respectively, are
     * \f${r, \theta, \phi, p_r, p_\theta, p_\phi}\f$ */
    polarvalues[0] = sqrt(inner_product3d(rvecprime,rvecprime));
    polarvalues[1] = acos(rvecprime[0] / polarvalues[0]);
    polarvalues[2] = atan2(-rvecprime[1], rvecprime[2]);
    //polarvalues[3] = inner_product3d(rvecprime, pvecprime) / polarvalues[0];
    /* FIX p_r = 0 */
    polarvalues[3] = 0;

    REAL8 rvecprime_x_xhat[3] = {0.}, rvecprime_x_xhat_x_rvecprime[3] = {0.};
    cross_product3d(rvecprime, Xhat, rvecprime_x_xhat);
    cross_product3d(rvecprime_x_xhat, rvecprime, rvecprime_x_xhat_x_rvecprime);

    polarvalues[4] = -inner_product3d(rvecprime_x_xhat_x_rvecprime, pvecprime)
                                / polarvalues[0] / sin(polarvalues[1]);
    polarvalues[5] = -inner_product3d(rvecprime_x_xhat, pvecprime);


    /* ********************************************************************* */  /* Finally, Differentiate Hamiltonian w.r.t. p_\phi, keeping p_r = 0 */
    /* ********************************************************************* */

    /* Populate the vector specifying the dynamical variables in mixed frames */
    memcpy( polarRPcartSvalues, cartvalues, 12*sizeof(REAL8));
    memcpy( polarRPcartSvalues, polarvalues, 6*sizeof(REAL8));

    /* Set up pointers for GSL */
    params.values  = polarRPcartSvalues;
    params.params  = funcParams;

    F.function = &GSLSpinPrecHamiltonianWrapperFordHdpphi;
    F.params   = &params;

    /* Now calculate omega. In the chosen co-ordinate system, */
    /* we need dH/dpphi to calculate this, i.e. varyParam = 5 */
    params.varyParam = 5;
    gslStatus = gsl_deriv_central( &F, polarvalues[5],
                    STEP_SIZE, &omega, &absErr );

    if ( gslStatus != GSL_SUCCESS )
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Failure in GSL function" );
        return REAL8_FAIL_NAN;
    }

    return omega;
}