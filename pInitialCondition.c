/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "pInitialCondition.h"
#include "pHamiltonian.h"
#include "pEnergyFlux.h"
#include "myLog.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_roots.h>

#define GSL_START \
          gsl_error_handler_t *saveGSLErrorHandler_; \
          saveGSLErrorHandler_ = gsl_set_error_handler_off();

#define GSL_END \
          gsl_set_error_handler( saveGSLErrorHandler_ );

/* Spin EOB root */
typedef
struct tagSEOBRootParams
{
    REAL8          values[12]; /**<< Dynamical variables, x, y, z, px, py, pz, S1x, S1y, S1z, S2x, S2y and S2z */
    SpinEOBParams *params;     /**<< Spin EOB parameters -- physical, pre-computed, etc. */
    REAL8          omega;      /**<< Orbital frequency */
    REAL8          e0;
}
SEOBRootParams;

static inline REAL8
CalculateCrossProduct( const INT4 i, const REAL8 a[], const REAL8 b[] )
{
    return a[(i+1)%3]*b[(i+2)%3] - a[(i+2)%3]*b[(i+1)%3];
}

static inline int NormalizeVector( REAL8 a[] )
{
    REAL8 norm = sqrt( a[0]*a[0] + a[1]*a[1] + a[2]*a[2] );
    
    a[0] /= norm;
    a[1] /= norm;
    a[2] /= norm;
    
    return CEV_SUCCESS;
}

static inline REAL8
CalculateDotProduct( const REAL8 a[], const REAL8 b[] )
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static int
CalculateRotationMatrix(
                        gsl_matrix *rotMatrix,  /**< OUTPUT, rotation matrix */
                        gsl_matrix *rotInverse, /**< OUTPUT, rotation matrix inversed */
                        REAL8       r[],        /**< position vector */
                        REAL8       v[],        /**< velocity vector */
                        REAL8       L[]         /**< orbital angular momentum */
)
{
    
    /** a, b, g are the angles alpha, beta and gamma */
    /* Use a, b and g to avoid shadowing gamma and getting a warning */
    REAL8 a, b, g;
    REAL8 cosa, sina, cosb, sinb, cosg, sing;
    
    /* Calculate the Euclidean Euler angles */
    
    /* We need to avoid a singularity if L is along z axis */
    if ( L[2] > 0.9999 )
    {
        a = b = g = 0.0;
    }
    else
    {
        a = GET_ATAN2( L[0], -L[1] );
        b = GET_ATAN2( sqrt( L[0]*L[0] + L[1]*L[1] ), L[2] );
        g = GET_ATAN2( r[2], v[2] );
    }
    
    if ( ( cosa = cos( a ) ) < 1.0e-16 ) cosa = 0.0;
    if ( ( sina = sin( a ) ) < 1.0e-16 ) sina = 0.0;
    if ( ( cosb = cos( b ) ) < 1.0e-16 ) cosb = 0.0;
    if ( ( sinb = sin( b ) ) < 1.0e-16 ) sinb = 0.0;
    if ( ( cosg = cos( g ) ) < 1.0e-16 ) cosg = 0.0;
    if ( ( sing = sin( g ) ) < 1.0e-16 ) sing = 0.0;
    
    /**
     * Implement the Rotation Matrix following the "x-convention"
     * 1. rotate about the z-axis by an angle a, rotation matrix Rz(a);
     * 2. rotate about the former x-axis (now x') by an angle b, rotation matrix Rx(b);
     * 3. rotate about the former z-axis (now z') by an algle g, rotation matrix Rz(g);
     */
    /* populate the matrix */
    // gsl_matrix_set( rotMatrix, 0, 0, cosg*cosa - cosb*sina*sing );
    // gsl_matrix_set( rotMatrix, 0, 1, cosg*sina + cosb*cosa*sing );
    // gsl_matrix_set( rotMatrix, 0, 2, sing*sinb );
    // gsl_matrix_set( rotMatrix, 1, 0, -sing*cosa - cosb*sina*cosg );
    // gsl_matrix_set( rotMatrix, 1, 1, -sing*sina + cosb*cosa*cosg );
    // gsl_matrix_set( rotMatrix, 1, 2, cosg*sinb );
    // gsl_matrix_set( rotMatrix, 2, 0, sinb*sina );
    // gsl_matrix_set( rotMatrix, 2, 1, -sinb*cosa );
    // gsl_matrix_set( rotMatrix, 2, 2, cosb );
    
    gsl_matrix_set(rotMatrix, 0, 0, r[0]);
	gsl_matrix_set(rotMatrix, 0, 1, r[1]);
	gsl_matrix_set(rotMatrix, 0, 2, r[2]);
	gsl_matrix_set(rotMatrix, 1, 0, v[0]);
	gsl_matrix_set(rotMatrix, 1, 1, v[1]);
	gsl_matrix_set(rotMatrix, 1, 2, v[2]);
	gsl_matrix_set(rotMatrix, 2, 0, L[0]);
	gsl_matrix_set(rotMatrix, 2, 1, L[1]);
	gsl_matrix_set(rotMatrix, 2, 2, L[2]);

    /* Now populate the transpose (which should be the inverse) */
    gsl_matrix_transpose_memcpy( rotInverse, rotMatrix );
    
    
    return CEV_SUCCESS;
}

static int
ApplyRotationMatrix(
                    gsl_matrix *rotMatrix, /**< rotation matrix */
                    REAL8      a[]         /**< OUTPUT, vector rotated */
)
{
    
    // gsl_vector *tmpVec1 = gsl_vector_alloc( 3 );
    // gsl_vector *tmpVec2 = gsl_vector_calloc( 3 );
    
    // gsl_vector_set( tmpVec1, 0, a[0] );
    // gsl_vector_set( tmpVec1, 1, a[1] );
    // gsl_vector_set( tmpVec1, 2, a[2] );
    
    // gsl_blas_dgemv( CblasNoTrans, 1.0,  rotMatrix, tmpVec1, 0.0, tmpVec2 );
    
    // a[0] = gsl_vector_get( tmpVec2, 0 );
    // a[1] = gsl_vector_get( tmpVec2, 1 );
    // a[2] = gsl_vector_get( tmpVec2, 2 );
    
    // gsl_vector_free( tmpVec1 );
    // gsl_vector_free( tmpVec2 );
	double b[3];
	gsl_vector_view tmpView1 = gsl_vector_view_array(a, 3);
	gsl_vector_view tmpView2 = gsl_vector_view_array(b, 3);
	gsl_vector     *tmpVec1 = &tmpView1.vector;
	gsl_vector     *tmpVec2 = &tmpView2.vector;

	gsl_blas_dgemv(CblasNoTrans, 1.0, rotMatrix, tmpVec1, 0.0, tmpVec2);

	gsl_vector_memcpy(tmpVec1, tmpVec2);

    return CEV_SUCCESS;
}

static int
XLALFindSphericalOrbit(const gsl_vector *x, /**<< Parameters requested by gsl root finder */
                       void *params,        /**<< Spin EOB parameters */
                       gsl_vector *f        /**<< Function values for the given parameters */
)
{
    SEOBRootParams *rootParams = (SEOBRootParams *) params;
    
    REAL8 py, pz, r, ptheta, pphi;
    
    /* Numerical derivative of Hamiltonian wrt given value */
    REAL8 dHdx, dHdpy, dHdpz;
    REAL8 dHdr, dHdptheta, dHdpphi;
    REAL8 e0;
    /* Populate the appropriate values */
    /* In the special theta=pi/2 phi=0 case, r is x */
    rootParams->values[0] = r  = gsl_vector_get( x, 0 );
    rootParams->values[4] = py = gsl_vector_get( x, 1 );
    rootParams->values[5] = pz = gsl_vector_get( x, 2 );
    e0 = rootParams->e0;
    // printf( "Values r = %.16e, py = %.16e, pz = %.16e\n", r, py, pz );
    
    ptheta = - r * pz;
    pphi   = r * py;
    
    /* dHdR */
    dHdx = XLALSpinHcapNumDerivWRTParam( 0, rootParams->values, rootParams->params );
    if ( IS_REAL8_FAIL_NAN( dHdx ) )
    {
        return CEV_FAILURE;
    }
    //printf( "dHdx = %.16e\n", dHdx );
    
    /* dHdPphi (I think we can use dHdPy in this coord system) */
    /* TODO: Check this is okay */
    dHdpy = XLALSpinHcapNumDerivWRTParam( 4, rootParams->values, rootParams->params );
    if ( IS_REAL8_FAIL_NAN( dHdpy ) )
    {
        return CEV_FAILURE;
    }
    
    /* dHdPtheta (I think we can use dHdPz in this coord system) */
    /* TODO: Check this is okay */
    dHdpz = XLALSpinHcapNumDerivWRTParam( 5, rootParams->values, rootParams->params );
    if ( IS_REAL8_FAIL_NAN( dHdpz ) )
    {
        return CEV_FAILURE;
    }
    
    /* Now convert to spherical polars */
    dHdr      = dHdx - dHdpy * pphi / (r*r) + dHdpz * ptheta / (r*r);
    dHdptheta = - dHdpz / r;
    dHdpphi   = dHdpy / r;
    
    /* populate the function vector */
    gsl_vector_set( f, 0, dHdr + e0/r/r);
    gsl_vector_set( f, 1, dHdptheta );
    gsl_vector_set( f, 2, dHdpphi - rootParams->omega );
    
    //printf( "Current funcvals = %.16e %.16e %.16e\n", gsl_vector_get( f, 0 ), gsl_vector_get( f, 1 ),
    //  gsl_vector_get( f, 2 )/*dHdpphi*/ );
    
    return CEV_SUCCESS;
}


/**
 * Wrapper for calculating specific derivative of the Hamiltonian in spherical co-ordinates,
 * dH/dr, dH/dptheta and dH/dpphi.
 * It only works for the specific co-ord system we use here
 */
static double GSLSpinHamiltonianDerivWrapper( double x,    /**<< Derivative at x */
                                              void  *params /**<< Function parameters */)
{
    HcapSphDeriv2Params *dParams = (HcapSphDeriv2Params *) params;

    REAL8 sphValues[12];
    REAL8 cartValues[12];

    REAL8 dHdr, dHdx, dHdpy, dHdpz;
    REAL8 r, ptheta, pphi;

    memcpy( sphValues, dParams->sphValues, sizeof( sphValues ) );
    sphValues[dParams->varyParam1] = x;

    SphericalToCartesian( cartValues, cartValues+3, sphValues, sphValues+3 );
    memcpy( cartValues+6, sphValues+6, 6*sizeof(REAL8) );

    r      = sphValues[0];
    ptheta = sphValues[4];
    pphi   = sphValues[5];

    /* Return the appropriate derivative according to varyParam2 */
    switch ( dParams->varyParam2 )
    {
        case 0:
            /* dHdr */
            dHdx  = XLALSpinHcapNumDerivWRTParam( 0, cartValues, dParams->params );
            dHdpy = XLALSpinHcapNumDerivWRTParam( 4, cartValues, dParams->params );
            dHdpz = XLALSpinHcapNumDerivWRTParam( 5, cartValues, dParams->params );

            dHdr      = dHdx - dHdpy * pphi / (r*r) + dHdpz * ptheta / (r*r);
            //printf( "dHdr = %.16e\n", dHdr );
            return dHdr;

            break;
        case 4:
            /* dHdptheta */
            dHdpz = XLALSpinHcapNumDerivWRTParam( 5, cartValues, dParams->params );
            return - dHdpz / r;
            break;
        case 5:
            /* dHdpphi */
            dHdpy = XLALSpinHcapNumDerivWRTParam( 4, cartValues, dParams->params );
            return dHdpy / r;
            break;
        default:
            PRINT_LOG_INFO(LOG_CRITICAL, "This option is not supported in the second derivative function!" );
            return REAL8_FAIL_NAN;
        break;
    }
}

/* Function to calculate the second derivative of the Hamiltonian. */
/* The derivatives are taken with respect to indices idx1, idx2    */
REAL8 XLALCalculateSphHamiltonianDeriv2(
                 const int      idx1,     /**<< Derivative w.r.t. index 1 */
                 const int      idx2,     /**<< Derivative w.r.t. index 2 */
                 const REAL8    values[], /**<< Dynamical variables in spherical coordinates */
                 SpinEOBParams *params    /**<< Spin EOB Parameters */
                 )
{
    static const REAL8 STEP_SIZE = 1.0e-4;

    REAL8 result;
    REAL8 absErr;

    HcapSphDeriv2Params dParams;

    gsl_function F;
    INT gslStatus;

    dParams.sphValues  = values;
    dParams.varyParam1 = idx1;
    dParams.varyParam2 = idx2;
    dParams.params     = params;

    // printf( " In second deriv function: values\n" );
    // for ( int i = 0; i < 12; i++ )
    // {
    //     printf( "%.16e ", values[i] );
    // }
    // printf( "\n" );
    
    F.function = GSLSpinHamiltonianDerivWrapper;
    F.params   = &dParams;

    /* GSL seemed to give weird answers - try my own code */
    // result = GSLSpinHamiltonianDerivWrapper( values[idx1] + STEP_SIZE, &dParams )
    //         - GSLSpinHamiltonianDerivWrapper( values[idx1] - STEP_SIZE, &dParams );
    // printf( "%.16e - %.16e / 2h\n", GSLSpinHamiltonianDerivWrapper( values[idx1] + STEP_SIZE, &dParams ), GSLSpinHamiltonianDerivWrapper( values[idx1] - STEP_SIZE, &dParams ) );

    //result = result / ( 2.*STEP_SIZE );
    
    
    gslStatus = gsl_deriv_central( &F, values[idx1], 
                        STEP_SIZE, &result, &absErr );
    // printf( "Second deriv abs err = %.16e\n", absErr );

    // printf( "RESULT = %.16e\n", result );

    if ( gslStatus != GSL_SUCCESS )
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Failure in GSL function" );
        return REAL8_FAIL_NAN;
    }

    return result;
}


/* Initial Condition Core function */
/* For nonprecessing case */
INT EOBInitialConditions(REAL8Vector    *initConds,
                         const REAL8    mass1,
                         const REAL8    mass2,
                         const REAL8    fMin,
                         const REAL8    ecc,
                         const REAL8    inc,
                         const REAL8    spin1[],
                         const REAL8    spin2[],
                         SpinEOBParams  *params)
{
    if (!initConds)
        return CEV_FAILURE;
    static const int lMax = 8;
    INT i;
    int tmpTortoise;
    /* non-zero eccentricity - Start frequency correction */
    // REAL8 ecc = params->eccentricity;
        
    REAL8 mTotal;
    REAL8 eta;
    REAL8 omega, v0;   /* Initial velocity and angular frequency */
    
    REAL8 ham;      /* Hamiltonian */
    
    REAL8 LnHat[3]; /* Initial orientation of angular momentum */
    REAL8 rHat[3];  /* Initial orientation of radial vector */
    REAL8 vHat[3];  /* Initial orientation of velocity vector */
    REAL8 Lhat[3];  /* Direction of relativistic ang mom */
    REAL8 qHat[3];
    REAL8 pHat[3];
    
    /* q and p vectors in Cartesian and spherical coords */
    REAL8 qCart[3], pCart[3];
    REAL8 qSph[3], pSph[3];
    
    /* We will need to manipulate the spin vectors */
    /* We will use temporary vectors to do this */
    REAL8 tmpS1[3];
    REAL8 tmpS2[3];
    REAL8 tmpS1Norm[3];
    REAL8 tmpS2Norm[3];
    
    REAL8Vector qCartVec, pCartVec;
    REAL8Vector s1Vec, s2Vec, s1VecNorm, s2VecNorm;
    REAL8Vector sKerr, sStar;
    REAL8       sKerrData[3], sStarData[3];
    REAL8       a = 0.; //, chiS, chiA;
    //REAL8       chi1, chi2;
    
    /* We will need a full values vector for calculating derivs of Hamiltonian */
    REAL8 sphValues[12];
    REAL8 cartValues[12];
    
    /* Matrices for rotating to the new basis set. */
    /* It is more convenient to calculate the ICs in a simpler basis */
    gsl_matrix *rotMatrix  = NULL;
    gsl_matrix *invMatrix  = NULL;
    gsl_matrix *rotMatrix2 = NULL;
    gsl_matrix *invMatrix2 = NULL;
    
    /* Root finding stuff for finding the spherical orbit */
    SEOBRootParams rootParams;
    const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver *rootSolver = NULL;
    
    gsl_multiroot_function rootFunction;
    gsl_vector *initValues  = NULL;
    gsl_vector *finalValues = NULL;
    int gslStatus;
    const int maxIter = 100;
    
    memset( &rootParams, 0, sizeof( rootParams ) );
    
    mTotal = mass1 + mass2;
    eta    = mass1 * mass2 / (mTotal * mTotal);
    memcpy( tmpS1, spin1, sizeof(tmpS1) );
    memcpy( tmpS2, spin2, sizeof(tmpS2) );
    memcpy( tmpS1Norm, spin1, sizeof(tmpS1Norm) );
    memcpy( tmpS2Norm, spin2, sizeof(tmpS2Norm) );
    for ( i = 0; i < 3; i++ )
    {
        tmpS1Norm[i] /= mTotal * mTotal;
        tmpS2Norm[i] /= mTotal * mTotal;
    }
    // eobVersion = params->seobCoeffs->eobVersion;
    /* We compute the ICs for the non-tortoise p, and convert at the end */
    tmpTortoise      = params->tortoise;
    params->tortoise = 0;
    
    EOBNonQCCoeffs *nqcCoeffs = NULL;
    nqcCoeffs = params->nqcCoeffs;
    
    /* STEP 1) Rotate to LNhat0 along z-axis and N0 along x-axis frame, where LNhat0 and N0 are initial normal to
     *         orbital plane and initial orbital separation;
     */
    
    /* Set the initial orbital ang mom direction. Taken from STPN code */
    LnHat[0] = GET_SIN(inc);
    LnHat[1] = 0.;
    LnHat[2] = GET_COS(inc);
    
    /* Set the radial direction - need to take care to avoid singularity if L is along z axis */
    if ( LnHat[2] > 0.9999 )
    {
        rHat[0] = 1.;
        rHat[1] = rHat[2] = 0.;
    }
    else
    {
        REAL8 theta0 = GET_ATAN( - LnHat[2] / LnHat[0] ); /* theta0 is between 0 and Pi */
        rHat[0] = GET_SIN( theta0 );
        rHat[1] = 0;
        rHat[2] = GET_COS( theta0 );
    }
    
    /* Now we can complete the triad */
    vHat[0] = CalculateCrossProduct( 0, LnHat, rHat );
    vHat[1] = CalculateCrossProduct( 1, LnHat, rHat );
    vHat[2] = CalculateCrossProduct( 2, LnHat, rHat );
    
    NormalizeVector( vHat );
    
    /* XXX Test code XXX */
    /*for ( i = 0; i < 3; i++ )
     {
     printf ( " LnHat[%d] = %.16e, rHat[%d] = %.16e, vHat[%d] = %.16e\n", i, LnHat[i], i, rHat[i], i, vHat[i] );
     }
     
     printf("\n\n" );
     for ( i = 0; i < 3; i++ )
     {
     printf ( " s1[%d] = %.16e, s2[%d] = %.16e\n", i, tmpS1[i], i, tmpS2[i] );
     }*/
    
    /* Allocate and compute the rotation matrices */
    rotMatrix = gsl_matrix_alloc( 3, 3 );
    invMatrix = gsl_matrix_alloc( 3, 3 );
    if ( !rotMatrix || !invMatrix )
    {
        if ( rotMatrix ) gsl_matrix_free( rotMatrix );
        if ( invMatrix ) gsl_matrix_free( invMatrix );
        return CEV_FAILURE;
    }
    
    if ( CalculateRotationMatrix( rotMatrix, invMatrix, rHat, vHat, LnHat ) == CEV_FAILURE )
    {
        gsl_matrix_free( rotMatrix );
        gsl_matrix_free( invMatrix );
        return CEV_FAILURE;
    }
    
    /* Rotate the orbital vectors and spins */
    ApplyRotationMatrix( rotMatrix, rHat );
    ApplyRotationMatrix( rotMatrix, vHat );
    ApplyRotationMatrix( rotMatrix, LnHat );
    ApplyRotationMatrix( rotMatrix, tmpS1 );
    ApplyRotationMatrix( rotMatrix, tmpS2 );
    ApplyRotationMatrix( rotMatrix, tmpS1Norm );
    ApplyRotationMatrix( rotMatrix, tmpS2Norm );
    
    /* XXX Test code XXX */
    /*printf( "\nAfter applying rotation matrix:\n\n" );
     for ( i = 0; i < 3; i++ )
     {
     printf ( " LnHat[%d] = %.16e, rHat[%d] = %.16e, vHat[%d] = %.16e\n", i, LnHat[i], i, rHat[i], i, vHat[i] );
     }
     
     printf("\n\n" );
     for ( i = 0; i < 3; i++ )
     {
     printf ( " s1[%d] = %.16e, s2[%d] = %.16e\n", i, tmpS1[i], i, tmpS2[i] );
     }*/
    
    /* STEP 2) After rotation in STEP 1, in spherical coordinates, phi0 and theta0 are given directly in Eq. (4.7),
     *         r0, pr0, ptheta0 and pphi0 are obtained by solving Eqs. (4.8) and (4.9) (using gsl_multiroot_fsolver).
     *         At this step, we find initial conditions for a spherical orbit without radiation reaction.
     */
    REAL8 fMinE = fMin;
    rootParams.e0 = 0.0;
    if (CODE_VERSION == 1)
        rootParams.e0 = ecc;
    else if (CODE_VERSION == 2)
    {
        fMinE /= pow(1-ecc*ecc, 1.5);
        //fMinE = fMinE * sqrt(1-ecc) / pow(1-ecc*ecc, 1.5)
    }

    /* Calculate the initial velocity from the given initial frequency */
    omega = CST_PI * mTotal * CST_MTSUN_SI * fMinE;
    v0    = GET_CBRT( omega );
    //print_log("omega0 = %e, v0 = %e\n", omega, v0);
    /* Given this, we can start to calculate the initial conditions */
    /* for spherical coords in the new basis */
    rootParams.omega  = omega;
    rootParams.params = params;

    // rootParams.values[0] = 1./(v0*v0);  /* Initial r */
    // rootParams.values[4] = v0;    /* Initial p */
    rootParams.values[0] = (1.-ecc) / (v0*v0);
    rootParams.values[4] = pow(1.-ecc, 0.5) * v0;
    memcpy( rootParams.values+6, tmpS1, sizeof( tmpS1 ) );
    memcpy( rootParams.values+9, tmpS2, sizeof( tmpS2 ) );
    
    //print_log( "ICs guess: r = %.16e, p = %.16e\n", rootParams.values[0], rootParams.values[4] );
    
    /* Initialise the gsl stuff */
    rootSolver = gsl_multiroot_fsolver_alloc( T, 3 );
    if ( !rootSolver )
    {
        gsl_matrix_free( rotMatrix );
        gsl_matrix_free( invMatrix );
        return CEV_FAILURE;
    }
    
    initValues = gsl_vector_calloc( 3 );
    if ( !initValues )
    {
        gsl_multiroot_fsolver_free( rootSolver );
        gsl_matrix_free( rotMatrix );
        gsl_matrix_free( invMatrix );
        return CEV_FAILURE;
    }
    
    gsl_vector_set( initValues, 0, rootParams.values[0] );
    gsl_vector_set( initValues, 1, rootParams.values[4] );
    
    rootFunction.f      = XLALFindSphericalOrbit;
    // rootFunction.f      = XLALFindEllipticalSphericalOrbit;
    rootFunction.n      = 3;
    rootFunction.params = &rootParams;
    
    gsl_multiroot_fsolver_set( rootSolver, &rootFunction, initValues );
    
    /* We are now ready to iterate to find the solution */
    i = 0;
    do
    {
        gslStatus = gsl_multiroot_fsolver_iterate( rootSolver );
        if ( gslStatus != GSL_SUCCESS )
        {
            print_warning( "Error in GSL iteration function!\n" );
            gsl_multiroot_fsolver_free( rootSolver );
            gsl_vector_free( initValues );
            gsl_matrix_free( rotMatrix );
            gsl_matrix_free( invMatrix );
            return CEV_FAILURE;
        }
        
        gslStatus = gsl_multiroot_test_residual( rootSolver->f, 1.0e-10 );
        i++;
    }
    while ( gslStatus == GSL_CONTINUE && i <= maxIter );
    
    if ( i > maxIter && gslStatus != GSL_SUCCESS )
    {
        gsl_multiroot_fsolver_free( rootSolver );
        gsl_vector_free( initValues );
        gsl_matrix_free( rotMatrix );
        gsl_matrix_free( invMatrix );
        return CEV_FAILURE;
    }
    
    finalValues = gsl_multiroot_fsolver_root( rootSolver );
    
    // print_log( "Spherical orbit conditions here given by the following:\n" );
    //  print_err( " x = %.16e, py = %.16e, pz = %.16e\n", gsl_vector_get( finalValues, 0 ),
    //  gsl_vector_get( finalValues, 1 ), gsl_vector_get( finalValues, 2 ) );
    
    memset( qCart, 0, sizeof(qCart) );
    memset( pCart, 0, sizeof(pCart) );
    
    qCart[0] = gsl_vector_get( finalValues, 0 );
    pCart[1] = gsl_vector_get( finalValues, 1 );
    pCart[2] = gsl_vector_get( finalValues, 2 );
    
    /* Free the GSL root finder, since we're done with it */
    gsl_multiroot_fsolver_free( rootSolver );
    gsl_vector_free( initValues );
    
    /* STEP 3) Rotate to L0 along z-axis and N0 along x-axis frame, where L0 is the initial orbital angular momentum
     *         and L0 is calculated using initial position and linear momentum obtained in STEP 2.
     */
    
    /* Now we can calculate the relativistic L and rotate to a new basis */
    memcpy( qHat, qCart, sizeof(qCart) );
    memcpy( pHat, pCart, sizeof(pCart) );
    
    NormalizeVector( qHat );
    NormalizeVector( pHat );
    
    Lhat[0] = CalculateCrossProduct( 0, qHat, pHat );
    Lhat[1] = CalculateCrossProduct( 1, qHat, pHat );
    Lhat[2] = CalculateCrossProduct( 2, qHat, pHat );
    
    NormalizeVector( Lhat );
    
    rotMatrix2 = gsl_matrix_alloc( 3, 3 );
    invMatrix2 = gsl_matrix_alloc( 3, 3 );
    
    if ( CalculateRotationMatrix( rotMatrix2, invMatrix2, qHat, pHat, Lhat ) == CEV_FAILURE )
    {
        gsl_matrix_free( rotMatrix );
        gsl_matrix_free( invMatrix );
        return CEV_FAILURE;
    }
    
    ApplyRotationMatrix( rotMatrix2, rHat );
    ApplyRotationMatrix( rotMatrix2, vHat );
    ApplyRotationMatrix( rotMatrix2, LnHat );
    ApplyRotationMatrix( rotMatrix2, tmpS1 );
    ApplyRotationMatrix( rotMatrix2, tmpS2 );
    ApplyRotationMatrix( rotMatrix2, tmpS1Norm );
    ApplyRotationMatrix( rotMatrix2, tmpS2Norm );
    ApplyRotationMatrix( rotMatrix2, qCart );
    ApplyRotationMatrix( rotMatrix2, pCart );
    
    /* STEP 4) In the L0-N0 frame, we calculate (dE/dr)|sph using Eq. (4.14), then initial dr/dt using Eq. (4.10),
     *         and finally pr0 using Eq. (4.15).
     */
    
    /* Now we can calculate the flux. Change to spherical co-ords */
    CartesianToSpherical( qSph, pSph, qCart, pCart );
    memcpy( sphValues, qSph, sizeof( qSph ) );
    memcpy( sphValues+3, pSph, sizeof( pSph ) );
    memcpy( sphValues+6, tmpS1, sizeof(tmpS1) );
    memcpy( sphValues+9, tmpS2, sizeof(tmpS2) );
    
    memcpy( cartValues, qCart, sizeof(qCart) );
    memcpy( cartValues+3, pCart, sizeof(pCart) );
    memcpy( cartValues+6, tmpS1, sizeof(tmpS1) );
    memcpy( cartValues+9, tmpS2, sizeof(tmpS2) );
    
    REAL8 dHdpphi, d2Hdr2, d2Hdrdpphi;
    REAL8 rDot, dHdpr, flux, dEdr;
    
    d2Hdr2 = XLALCalculateSphHamiltonianDeriv2( 0, 0, sphValues, params );
    d2Hdrdpphi = XLALCalculateSphHamiltonianDeriv2( 0, 5, sphValues, params );
    dHdpphi = XLALSpinHcapNumDerivWRTParam( 4, cartValues, params ) / sphValues[0];
    dEdr  = - dHdpphi * d2Hdr2 / d2Hdrdpphi;
        
    if ( d2Hdr2 != 0.0 && ecc==0.0 )
    {
        /* We will need to calculate the Hamiltonian to get the flux */
        s1Vec.length = s2Vec.length = s1VecNorm.length = s2VecNorm.length = sKerr.length = sStar.length = 3;
        s1Vec.data = tmpS1;
        s2Vec.data = tmpS2;
        s1VecNorm.data = tmpS1Norm;
        s2VecNorm.data = tmpS2Norm;
        sKerr.data = sKerrData;
        sStar.data = sStarData;
        
        qCartVec.length = pCartVec.length = 3;
        qCartVec.data   = qCart;
        pCartVec.data   = pCart;
        
        //chi1 = tmpS1[0]*LnHat[0] + tmpS1[1]*LnHat[1] + tmpS1[2]*LnHat[2];
        //chi2 = tmpS2[0]*LnHat[0] + tmpS2[1]*LnHat[1] + tmpS2[2]*LnHat[2];
        
        // print_debug( "m1 = %g, m2 = %g\n", mass1, mass2 );
        // print_debug( "eta = %g\n", eta);
        // print_debug( "qCartVec = (%g, %g, %g)\n", qCart[0], qCart[1], qCart[2]);
        // print_debug( "pCartVec = (%g, %g, %g)\n", pCart[0], pCart[1], pCart[2]);
        // print_debug( "s1Vec = (%g, %g, %g)\n", tmpS1[0], tmpS1[1], tmpS1[2]);
        // print_debug( "s2Vec = (%g, %g, %g)\n", tmpS2[0], tmpS2[1], tmpS2[2]);
        
        //chiS = 0.5 * ( chi1 / (mass1*mass1) + chi2 / (mass2*mass2) );
        //chiA = 0.5 * ( chi1 / (mass1*mass1) - chi2 / (mass2*mass2) );
        
        XLALSimIMRSpinEOBCalculateSigmaStar( &sKerr, mass1, mass2, &s1Vec, &s2Vec );
        XLALSimIMRSpinEOBCalculateSigmaKerr( &sStar, mass1, mass2, &s1Vec, &s2Vec );
        // print_debug( "sigStar = (%g, %g, %g)\n", sStarData[0], sStarData[1], sStarData[2]);
        // print_debug( "sigKerr = (%g, %g, %g)\n", sStarData[0], sStarData[1], sStarData[2]);

        /* The a in the flux has been set to zero, but not in the Hamiltonian */
        a = sqrt(sKerr.data[0]*sKerr.data[0] + sKerr.data[1]*sKerr.data[1] + sKerr.data[2]*sKerr.data[2]);
        //XLALSimIMREOBCalcSpinFacWaveformCoefficients( params->eobParams->hCoeffs, mass1, mass2, eta, /*a*/0.0, chiS, chiA );
        //XLALSimIMRCalculateSpinEOBHCoeffs( params->seobCoeffs, eta, a );
        ham = EOBHamiltonian( eta, &qCartVec, &pCartVec, &s1VecNorm, &s2VecNorm, &sKerr, &sStar, params->tortoise, params->seobCoeffs );
        // print_debug( "hamiltonian at this point is %.16e\n", ham );
        
        /* And now, finally, the flux */
        REAL8Vector polarDynamics;
        REAL8       polarData[4];
        REAL8Vector cartDynamics;
        REAL8       cartData[6];
        
        polarDynamics.length = 4;
        polarDynamics.data = polarData;
        cartDynamics.length = 6;
        cartDynamics.data = cartData;
        
        polarData[0] = qSph[0];
        polarData[1] = 0.;
        polarData[2] = pSph[0];
        polarData[3] = pSph[2];
        memcpy(cartData, qCart, sizeof(qCart));
        memcpy(cartData+3, pCart, sizeof(pCart));
        // print_debug("polarData = (%g, %g, %g, %g)\n", 
        //     polarData[0], polarData[1], polarData[2], polarData[3]);
        // print_debug("ham = %g\n", ham);
        REAL8 tmpdvalues[4] = {0.,omega,0.,0.};
        flux  = XLALInspiralSpinFactorizedFlux( &polarDynamics, &cartDynamics, nqcCoeffs, omega, 0, qSph[0]*omega, params, ham, lMax);
        // flux  = InspiralSpinFactorizedFlux_elip( &polarDynamics, polarData, tmpdvalues, nqcCoeffs, omega, params, ham, lMax, eobVersion );
        flux  = flux / eta;
        if (ecc != 0.0)
            rDot = 0.0;
        else
            rDot  = - flux / dEdr;
        // print_debug("flux = %g, dEdr = %g, rDot = %g\n", flux, dEdr, rDot);
        /* We now need dHdpr - we take it that it is safely linear up to a pr of 1.0e-3 */
        cartValues[3] = 1.0e-3;
        dHdpr = XLALSpinHcapNumDerivWRTParam( 3, cartValues, params );
        /*printf( "Ingredients going into prDot:\n" );
         printf( "flux = %.16e, dEdr = %.16e, dHdpr = %.16e\n", flux, dEdr, dHdpr );*/
        
        /* We can now calculate what pr should be taking into account the flux */
        pSph[0] = rDot / (dHdpr / cartValues[3] );
    }
    else
    {
        /* Since d2Hdr2 has evaluated to zero, we cannot do the above. Just set pr to zero */
        //printf( "d2Hdr2 is zero!\n" );
        pSph[0] = 0;
    }
    
    /* Now we are done - convert back to cartesian coordinates ) */
    SphericalToCartesian( qCart, pCart, qSph, pSph );
    PRINT_LOG_INFO(LOG_DEBUG, "Sph initial condition : r = (%e,%e,%e), p = (%e,%e,%e)", qSph[0], qSph[1], qSph[2], pSph[0], pSph[1], pSph[2]);
    /* STEP 5) Rotate back to the original inertial frame by inverting the rotation of STEP 3 and then
     *         inverting the rotation of STEP 1.
     */
    
    /* Undo rotations to get back to the original basis */
    /* Second rotation */
    ApplyRotationMatrix( invMatrix2, rHat );
    ApplyRotationMatrix( invMatrix2, vHat );
    ApplyRotationMatrix( invMatrix2, LnHat );
    ApplyRotationMatrix( invMatrix2, tmpS1 );
    ApplyRotationMatrix( invMatrix2, tmpS2 );
    ApplyRotationMatrix( invMatrix2, tmpS1Norm );
    ApplyRotationMatrix( invMatrix2, tmpS2Norm );
    ApplyRotationMatrix( invMatrix2, qCart );
    ApplyRotationMatrix( invMatrix2, pCart );
    
    /* First rotation */
    ApplyRotationMatrix( invMatrix, rHat );
    ApplyRotationMatrix( invMatrix, vHat );
    ApplyRotationMatrix( invMatrix, LnHat );
    ApplyRotationMatrix( invMatrix, tmpS1 );
    ApplyRotationMatrix( invMatrix, tmpS2 );
    ApplyRotationMatrix( invMatrix, tmpS1Norm );
    ApplyRotationMatrix( invMatrix, tmpS2Norm );
    ApplyRotationMatrix( invMatrix, qCart );
    ApplyRotationMatrix( invMatrix, pCart );
    

    /* If required, apply the tortoise transform */
    if ( tmpTortoise )
    {
        REAL8 r = sqrt(qCart[0]*qCart[0] + qCart[1]*qCart[1] + qCart[2]*qCart[2] );
        REAL8 deltaR = XLALSimIMRSpinEOBHamiltonianDeltaR( params->seobCoeffs, r, eta, a );
        REAL8 deltaT = XLALSimIMRSpinEOBHamiltonianDeltaT( params->seobCoeffs, r, eta, a );
        REAL8 csi    = sqrt( deltaT * deltaR )/(r*r + a*a);
        
        REAL8 pr = (qCart[0]*pCart[0] + qCart[1]*pCart[1] + qCart[2]*pCart[2])/r;
        params->tortoise = tmpTortoise;
        
        //printf( "Applying the tortoise to p (csi = %.26e)\n", csi );
        
        for ( i = 0; i < 3; i++ )
        {
            pCart[i] = pCart[i] + qCart[i] * pr * (csi - 1.) / r;
        }
    }

    if (CODE_VERSION == 2)
    {
        qCart[0] /= 1+ecc;
        pCart[1] *= (1+ecc);
    }

    /* Now copy the initial conditions back to the return vector */
    memcpy( initConds->data, qCart, sizeof(qCart) );
    memcpy( initConds->data+3, pCart, sizeof(pCart) );
    memcpy( initConds->data+6, tmpS1Norm, sizeof(tmpS1Norm) );
    memcpy( initConds->data+9, tmpS2Norm, sizeof(tmpS2Norm) );
    
    gsl_matrix_free(rotMatrix2);
    gsl_matrix_free(invMatrix2);
    
    gsl_matrix_free(rotMatrix);
    gsl_matrix_free(invMatrix);
    
    //printf( "THE FINAL INITIAL CONDITIONS:\n");
    /*printf( " %.16e %.16e %.16e\n%.16e %.16e %.16e\n%.16e %.16e %.16e\n%.16e %.16e %.16e\n", initConds->data[0], initConds->data[1], initConds->data[2],
     initConds->data[3], initConds->data[4], initConds->data[5], initConds->data[6], initConds->data[7], initConds->data[8],
     initConds->data[9], initConds->data[10], initConds->data[11] );*/
    
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
static REAL8 scale1 = 1, scale2 = 2, scale3 = 200;

static int XLALRobustDerivative(const gsl_function * F, double x, double h,
				  double *result, double *absErr)
{
    // Take the derivative
    gsl_deriv_central(F, x,h, result, absErr);
    // We check the estimate of the error
    REAL8 frac= 0.01; // Adjust this to change how accurtate we demand the derivative to be
    if (fabs(*absErr)>frac*fabs(*result))
    {
        REAL8 temp1 = 0.0;
        REAL8 temp2 =0.0;
        REAL8 absErr1 = 0.0;
        REAL8 absErr2 = 0.0;
        UINT deriv_ok = 0;
        UINT n = 1;
        while (!deriv_ok && n<=10)
        {
            PRINT_LOG_INFO(LOG_WARNING, "Warning: second derivative computation went wrong. Trying again");
            gsl_deriv_central(F, x,2*n*h, &temp1, &absErr1);
            gsl_deriv_central(F, x,h/(2*n), &temp2, &absErr2);
            if (fabs(absErr1)/fabs(temp1)<fabs(absErr2)/fabs(temp2) && fabs(absErr1)<frac*fabs(temp1))
            {
                *result = temp1;
                *absErr = absErr1;
                deriv_ok = 1;
                break;
            }
            else if (fabs(absErr1)/fabs(temp1)>fabs(absErr2)/fabs(temp2) && fabs(absErr2)<frac*fabs(temp2))
            {
                *result = temp2;
                *absErr = absErr2;
                deriv_ok = 1;
                break;
            }
            else
            {
                n+=1;
            }
        }
        if (!deriv_ok){
        PRINT_LOG_INFO(LOG_CRITICAL, "The computation of the second derivative of H in initial data has failed!");
        return CEV_FAILURE;
        }
    }
    return CEV_FAILURE;
}

static double
GSLSpinHamiltonianDerivWrapperPrec(double x,	/**<< Derivative at x */
                                   void *params /**<< Function parameters */
)
{
	// int		debugPK = 0;
    int i;
	HcapSphDeriv2Params *dParams = (HcapSphDeriv2Params *) params;
	REAL8		mTotal = dParams->params->m1 + dParams->params->m2;
	REAL8		sphValues[12];
	REAL8		cartValues[12];

	REAL8		dHdr    , dHdx, dHdpy, dHdpz;
	REAL8		r       , ptheta, pphi;

	memcpy(sphValues, dParams->sphValues, sizeof(sphValues));
	sphValues[dParams->varyParam1] = x;

	SphericalToCartesian(cartValues, cartValues + 3, sphValues, sphValues + 3);
	memcpy(cartValues + 6, sphValues + 6, 6 * sizeof(REAL8));

	r = sphValues[0];
	ptheta = sphValues[4];
	pphi = sphValues[5];

	/* New code to compute derivatives w.r.t. cartesian variables */
	REAL8		tmpDValues[14];
	int status;
	for ( i = 0; i < 3; i++) {
		cartValues[i + 6] /= mTotal * mTotal;
		cartValues[i + 9] /= mTotal * mTotal;
	}
    UINT oldignoreflux = dParams->params->ignoreflux;
    dParams->params->ignoreflux = 1;
    status = XLALSpinPrecHcapNumericalDerivative(0, cartValues, tmpDValues, dParams->params);
    dParams->params->ignoreflux = oldignoreflux;
	for ( i = 0; i < 3; i++) {
		cartValues[i + 6] *= mTotal * mTotal;
		cartValues[i + 9] *= mTotal * mTotal;
	}

	double		rvec    [3] = {cartValues[0], cartValues[1], cartValues[2]};
	double		pvec    [3] = {cartValues[3], cartValues[4], cartValues[5]};
	double		chi1vec [3] = {cartValues[6], cartValues[7], cartValues[8]};
	double		chi2vec [3] = {cartValues[9], cartValues[10], cartValues[11]};
	double		Lvec    [3] = {CalculateCrossProduct(0, rvec, pvec), CalculateCrossProduct(1, rvec, pvec), CalculateCrossProduct(2, rvec, pvec)};
	double		theta1 = acos(CalculateDotProduct(chi1vec, Lvec) / sqrt(CalculateDotProduct(chi1vec, chi1vec)) / sqrt(CalculateDotProduct(Lvec, Lvec)));
	double		theta2 = acos(CalculateDotProduct(chi2vec, Lvec) / sqrt(CalculateDotProduct(chi2vec, chi2vec)) / sqrt(CalculateDotProduct(Lvec, Lvec)));
	// if (debugPK) {
	// 	XLAL_PRINT_INFO("In GSLSpinHamiltonianDerivWrapperPrec");
	// 	XLAL_PRINT_INFO("rvec = %.16e %.16e %.16e\n", rvec[0], rvec[1], rvec[2]);
	// 	XLAL_PRINT_INFO("pvec = %.16e %.16e %.16e\n", pvec[0], pvec[1], pvec[2]);
	// 	XLAL_PRINT_INFO("theta1 = %.16e\n", theta1);
	// 	XLAL_PRINT_INFO("theta2 = %.16e\n", theta2);
	// }
	if (theta1 > 1.0e-6 && theta2 >= 1.0e-6) 
    {
		switch (dParams->varyParam2) {
		case 0:
			/* dHdr */
			dHdx = -tmpDValues[3];
			//XLALSpinPrecHcapNumDerivWRTParam(0, cartValues, dParams->params);
			dHdpy = tmpDValues[1];
			//XLALSpinPrecHcapNumDerivWRTParam(4, cartValues, dParams->params);
			dHdpz = tmpDValues[2];
			//XLALSpinPrecHcapNumDerivWRTParam(5, cartValues, dParams->params);

			dHdr = dHdx - dHdpy * pphi / (r * r) + dHdpz * ptheta / (r * r);
			//XLAL_PRINT_INFO("dHdr = %.16e\n", dHdr);
			return dHdr;

			break;
		case 4:
			/* dHdptheta */
			dHdpz = tmpDValues[2];
			//XLALSpinPrecHcapNumDerivWRTParam(5, cartValues, dParams->params);
			return -dHdpz / r;
			break;
		case 5:
			/* dHdpphi */
			dHdpy = tmpDValues[1];
			//XLALSpinPrecHcapNumDerivWRTParam(4, cartValues, dParams->params);
			return dHdpy / r;
			break;
		default:
			PRINT_LOG_INFO(LOG_CRITICAL, "This option is not supported in the second derivative function!");
			return REAL8_FAIL_NAN;
			break;
		}
	} else {
		switch (dParams->varyParam2) {
		case 0:
			/* dHdr */
            dHdx = XLALSpinPrecHcapNumDerivWRTParam(0, cartValues, dParams->params);
			dHdpy = XLALSpinPrecHcapNumDerivWRTParam(4, cartValues, dParams->params);
			dHdpz = XLALSpinPrecHcapNumDerivWRTParam(5, cartValues, dParams->params);
			dHdr = dHdx - dHdpy * pphi / (r * r) + dHdpz * ptheta / (r * r);
			//XLAL_PRINT_INFO("dHdr = %.16e\n", dHdr);
			return dHdr;

			break;
		case 4:
			/* dHdptheta */
            dHdpz = XLALSpinPrecHcapNumDerivWRTParam(5, cartValues, dParams->params);
			return -dHdpz / r;
			break;
		case 5:
			/* dHdpphi */
		        dHdpy = XLALSpinPrecHcapNumDerivWRTParam(4, cartValues, dParams->params);
			return dHdpy / r;
			break;
		default:
			PRINT_LOG_INFO(LOG_CRITICAL, "This option is not supported in the second derivative function!");
			return REAL8_FAIL_NAN;
			break;
		}
	}
}


/**
 * Function to calculate the second derivative of the Hamiltonian.
* The derivatives are taken with respect to indices idx1, idx2    */
static REAL8
XLALCalculateSphHamiltonianDeriv2Prec(
				  const int idx1,	/**<< Derivative w.r.t. index 1 */
				  const int idx2,	/**<< Derivative w.r.t. index 2 */
				  const REAL8 values[],	/**<< Dynamical variables in spherical coordinates */
				  SpinEOBParams * params	/**<< Spin EOB Parameters */
)
{

	static const REAL8 STEP_SIZE = 3.0e-3;
	REAL8		result;
	REAL8 absErr;

	HcapSphDeriv2Params dParams;

	gsl_function	F;
	INT4 gslStatus;

	dParams.sphValues = values;
	dParams.varyParam1 = idx1;
	dParams.varyParam2 = idx2;
	dParams.params = params;
        // dParams.use_optimized = use_optimized;

	/*
	 * XLAL_PRINT_INFO( " In second deriv function: values\n" ); for ( int i = 0;
	 * i < 12; i++ ) { XLAL_PRINT_INFO( "%.16e ", values[i] ); } XLAL_PRINT_INFO( "\n" );
	 */
	F.function = GSLSpinHamiltonianDerivWrapperPrec;
	F.params = &dParams;

	/* GSL seemed to give weird answers - try my own code */
	/*
	 * result = GSLSpinHamiltonianDerivWrapperPrec( values[idx1] + STEP_SIZE,
	 * &dParams ) - GSLSpinHamiltonianDerivWrapperPrec( values[idx1] -
	 * STEP_SIZE, &dParams ); XLAL_PRINT_INFO( "%.16e - %.16e / 2h\n",
	 * GSLSpinHamiltonianDerivWrapperPrec( values[idx1] + STEP_SIZE, &dParams
	 * ), GSLSpinHamiltonianDerivWrapperPrec( values[idx1] - STEP_SIZE,
	 * &dParams ) );
	 *
	 * result = result / ( 2.*STEP_SIZE );
	 */

	//XLAL_CALLGSL(gslStatus = gsl_deriv_central(&F, values[idx1],
	//				      STEP_SIZE, &result, &absErr));
	XLALRobustDerivative(&F, values[idx1],
			     STEP_SIZE, &result, &absErr);
	/*
	if (gslStatus != GSL_SUCCESS) {
		XLALPrintError("XLAL Error %s - Failure in GSL function\n", __func__);
		XLAL_ERROR_REAL8(XLAL_EDOM);
	}
	*/
	//XLAL_PRINT_INFO("Second deriv abs err = %.16e\n", absErr);

	//XLAL_PRINT_INFO("RESULT = %.16e\n", result);
	return result;
}


static int
XLALFindSphericalOrbitPrec(
                           const gsl_vector * x,	/**<< Parameters requested by gsl root finder */
                           void *params,	        /**<< Spin EOB parameters */
                           gsl_vector * f	      /**<< Function values for the given parameters*/
)
{
    SEOBRootParams *rootParams = (SEOBRootParams *) params;
    REAL8		mTotal = rootParams->params->m1 + rootParams->params->m2;
    REAL8		py      , pz, r, ptheta, pphi;

    /* Numerical derivative of Hamiltonian wrt given value */
    REAL8		dHdx    , dHdpy, dHdpz;
    REAL8		dHdr    , dHdptheta, dHdpphi;
    REAL8 e0;
    int i;
    /* Populate the appropriate values */
    /* In the special theta=pi/2 phi=0 case, r is x */

    REAL8 temp = gsl_vector_get(x, 0)/scale1;
    REAL8 prefactor = 1.0;
    if (temp  < 0.0)
    {
        prefactor=-1.0;
    }

    rootParams->values[0] = r =  sqrt(temp*temp+36.0);
    rootParams->values[4] = py = gsl_vector_get(x, 1)/scale2;
    rootParams->values[5] = pz = gsl_vector_get(x, 2)/scale3;
    //printf("r is %.17f\n",r);
    if(isnan(rootParams->values[0])) 
    {
        rootParams->values[0] = 100.;
    }
    if(isnan(rootParams->values[4])) 
    {
        rootParams->values[4] = 0.1;
    }
    if(isnan(rootParams->values[5])) 
    {
        rootParams->values[5] = 0.01;
    }
    // XLAL_PRINT_INFO("%3.10f %3.10f %3.10f %3.10f %3.10f %3.10f\n",
    //     rootParams->values[0], rootParams->values[1], rootParams->values[2],
    //     rootParams->values[3], rootParams->values[4], rootParams->values[5]);
    // XLAL_PRINT_INFO("Values r = %.16e, py = %.16e, pz = %.16e\n", r, py, pz);
    // fflush(NULL);

    ptheta = -r * pz;
    pphi = r * py;

    // XLAL_PRINT_INFO("Input Values r = %.16e, py = %.16e, pz = %.16e\n pthetha = %.16e pphi = %.16e\n", r, py, pz, ptheta, pphi);

    /* dH by dR and dP */
    REAL8	tmpDValues[14];
    int status;
    for ( i = 0; i < 3; i++) 
    {
        rootParams->values[i + 6] /= mTotal * mTotal;
        rootParams->values[i + 9] /= mTotal * mTotal;
    }
    UINT oldignoreflux = rootParams->params->ignoreflux;
    rootParams->params->ignoreflux = 1;
    status = XLALSpinPrecHcapNumericalDerivative(0, rootParams->values, tmpDValues, rootParams->params);
    rootParams->params->ignoreflux = oldignoreflux;
    for ( i = 0; i < 3; i++) 
    {
        rootParams->values[i + 6] *= mTotal * mTotal;
        rootParams->values[i + 9] *= mTotal * mTotal;
    }
    REAL8 rvec[3] =
        {rootParams->values[0], rootParams->values[1], rootParams->values[2]};
    REAL8 pvec[3] =
        {rootParams->values[3], rootParams->values[4], rootParams->values[5]};
    REAL8 chi1vec[3] =
        {rootParams->values[6], rootParams->values[7], rootParams->values[8]};
    REAL8 chi2vec[3] =
        {rootParams->values[9], rootParams->values[10], rootParams->values[11]};
    REAL8 Lvec[3] = {CalculateCrossProduct(0, rvec, pvec),
    CalculateCrossProduct(1, rvec, pvec), CalculateCrossProduct(2, rvec, pvec)};
    REAL8 theta1 = acos(CalculateDotProduct(chi1vec, Lvec)
                / sqrt(CalculateDotProduct(chi1vec, chi1vec))
                / sqrt(CalculateDotProduct(Lvec, Lvec)));
    REAL8 theta2 = acos(CalculateDotProduct(chi2vec, Lvec)
                / sqrt(CalculateDotProduct(chi2vec, chi2vec))
                / sqrt(CalculateDotProduct(Lvec, Lvec)));

    // XLAL_PRINT_INFO("rvec = %.16e %.16e %.16e\n", rvec[0], rvec[1], rvec[2]);
    // XLAL_PRINT_INFO("pvec = %.16e %.16e %.16e\n", pvec[0], pvec[1], pvec[2]);
    // XLAL_PRINT_INFO("theta1 = %.16e\n", theta1);
    // XLAL_PRINT_INFO("theta2 = %.16e\n", theta2);

    if (theta1 > 1.0e-6 && theta2 >= 1.0e-6) 
    {
        dHdx = -tmpDValues[3];
        dHdpy = tmpDValues[1];
        dHdpz = tmpDValues[2];
    } else 
    {
        //rootParams->values[5] = 0.;
        //rootParams->values[6] = 0.;
        //rootParams->values[7] = 0.;
        //rootParams->values[8] = sqrt(CalculateDotProductPrec(chi1vec, chi1vec));
        //rootParams->values[9] = 0.;
        //rootParams->values[10]= 0.;
        //rootParams->values[11]= sqrt(CalculateDotProductPrec(chi2vec, chi2vec));
        dHdx = XLALSpinPrecHcapNumDerivWRTParam(0,
                            rootParams->values, rootParams->params);
        dHdpy = XLALSpinPrecHcapNumDerivWRTParam(4,
                                rootParams->values, rootParams->params);
        dHdpz = XLALSpinPrecHcapNumDerivWRTParam(5,
                                rootParams->values, rootParams->params);
    }
    if (IS_REAL8_FAIL_NAN(dHdx)) { return CEV_FAILURE; }
    if (IS_REAL8_FAIL_NAN(dHdpy)) { return CEV_FAILURE; }
    if (IS_REAL8_FAIL_NAN(dHdpz)) { return CEV_FAILURE; }
    // XLAL_PRINT_INFO("dHdx = %.16e, dHdpy = %.16e, dHdpz = %.16e\n", dHdx, dHdpy, dHdpz);

    /* Now convert to spherical polars */
    dHdr      = dHdx - dHdpy * pphi / (r * r) + dHdpz * ptheta / (r * r);
    dHdptheta = -dHdpz / r;
    dHdpphi   = dHdpy / r;
    e0 = rootParams->e0;
    // XLAL_PRINT_INFO("dHdr = %.16e dHdptheta = %.16e dHdpphi = %.16e\n",
    //         dHdr, dHdptheta, dHdpphi);
    /* populate the function vector */
    gsl_vector_set(f, 0, dHdr + e0/r/r);
    gsl_vector_set(f, 1, dHdptheta);
    gsl_vector_set(f, 2, dHdpphi - rootParams->omega);

    // XLAL_PRINT_INFO("Current funcvals = %.16e %.16e %.16e\n",
    // gsl_vector_get(f, 0), gsl_vector_get(f, 1), gsl_vector_get(f, 2));

    /* Rescale back */


    rootParams->values[0] = prefactor*scale1*(sqrt(rootParams->values[0]*rootParams->values[0]-36.0));
    rootParams->values[4] *= scale2;
    rootParams->values[5] *= scale3;

	return CEV_SUCCESS;
}

INT
EOBInitialConditionsPrec(
                    REAL8Vector * initConds,	/**<< OUTPUT, Initial dynamical variables */
                    const REAL8 mass1,	/**<< mass 1 */
                    const REAL8 mass2,	/**<< mass 2 */
                    const REAL8 fMin,	/**<< Initial frequency (given) */
                    const REAL8 e0,
                    const REAL8 inc,	/**<< Inclination */
                    const REAL8 spin1[],	/**<< Initial spin vector 1 */
                    const REAL8 spin2[],	/**<< Initial spin vector 2 */
                    SpinEOBParams * params	/**<< Spin EOB parameters */
)
{
    // print_debug("Inside the XLALSimIMRSpinEOBInitialConditionsPrec function!\n");
    // print_debug(
    // "Inputs: m1 = %.16e, m2 = %.16e, fMin = %.16e, inclination = %.16e\n",
    //   mass1, mass2, (double)fMin, (double)inc);
    // print_debug("Inputs: mSpin1 = {%.16e, %.16e, %.16e}\n",
    //   spin1[0], spin1[1], spin1[2]);
    // print_debug("Inputs: mSpin2 = {%.16e, %.16e, %.16e}\n",
    //   spin2[0], spin2[1], spin2[2]);
	// 	fflush(NULL);
	static const int lMax = 8;

	int		i;

	/* Variable to keep track of whether the user requested the tortoise */
	int		tmpTortoise;

	UINT		SpinAlignedEOBversion;

	REAL8		mTotal;
	REAL8		eta;
	REAL8		omega   , v0;	/* Initial velocity and angular
					 * frequency */

	REAL8		ham;	/* Hamiltonian */

	REAL8		LnHat    [3];	/* Initial orientation of angular
					 * momentum */
	REAL8		rHat     [3];	/* Initial orientation of radial
					 * vector */
	REAL8		vHat     [3];	/* Initial orientation of velocity
					 * vector */
	REAL8		Lhat     [3];	/* Direction of relativistic ang mom */
	REAL8		qHat     [3];
	REAL8		pHat     [3];

	/* q and p vectors in Cartesian and spherical coords */
	REAL8		qCart    [3], pCart[3];
	REAL8		qSph     [3], pSph[3];

	/* We will need to manipulate the spin vectors */
	/* We will use temporary vectors to do this */
	REAL8		tmpS1    [3];
	REAL8		tmpS2    [3];
	REAL8		tmpS1Norm[3];
	REAL8		tmpS2Norm[3];

	REAL8Vector	qCartVec, pCartVec;
	REAL8Vector	s1Vec, s2Vec, s1VecNorm, s2VecNorm;
	REAL8Vector	sKerr, sStar;
	REAL8		sKerrData[3], sStarData[3];
	REAL8		a = 0.;
	//, chiS, chiA;
	//REAL8 chi1, chi2;

	/*
	 * We will need a full values vector for calculating derivs of
	 * Hamiltonian
	 */
	REAL8		sphValues[12];
	REAL8		cartValues[12];

	/* Matrices for rotating to the new basis set. */
	/* It is more convenient to calculate the ICs in a simpler basis */
	gsl_matrix     *rotMatrix = NULL;
	gsl_matrix     *invMatrix = NULL;
	gsl_matrix     *rotMatrix2 = NULL;
	gsl_matrix     *invMatrix2 = NULL;

	/* Root finding stuff for finding the spherical orbit */
	SEOBRootParams	rootParams;
	const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrids;
	gsl_multiroot_fsolver *rootSolver = NULL;

	gsl_multiroot_function rootFunction;
	gsl_vector     *initValues = NULL;
	gsl_vector     *finalValues = NULL;
	INT gslStatus;
    INT cntGslNoProgress = 0, MAXcntGslNoProgress = 5;
    //INT cntGslNoProgress = 0, MAXcntGslNoProgress = 50;
    REAL8 multFacGslNoProgress = 3./5.;
	//const int	maxIter = 2000;
	const int	maxIter = 10000;

	memset(&rootParams, 0, sizeof(rootParams));

	mTotal = mass1 + mass2;
	eta = mass1 * mass2 / (mTotal * mTotal);
	memcpy(tmpS1, spin1, sizeof(tmpS1));
	memcpy(tmpS2, spin2, sizeof(tmpS2));
	memcpy(tmpS1Norm, spin1, sizeof(tmpS1Norm));
	memcpy(tmpS2Norm, spin2, sizeof(tmpS2Norm));
	for (i = 0; i < 3; i++) {
		tmpS1Norm[i] /= mTotal * mTotal;
		tmpS2Norm[i] /= mTotal * mTotal;
	}
	SpinAlignedEOBversion = 4;
	/* We compute the ICs for the non-tortoise p, and convert at the end */
	tmpTortoise = params->tortoise;
	params->tortoise = 0;

	EOBNonQCCoeffs *nqcCoeffs = NULL;
	nqcCoeffs = params->nqcCoeffs;

	/*
	 * STEP 1) Rotate to LNhat0 along z-axis and N0 along x-axis frame,
	 * where LNhat0 and N0 are initial normal to orbital plane and
	 * initial orbital separation;
	 */

	/* Set the initial orbital ang mom direction. Taken from STPN code */
	LnHat[0] = sin(inc);
	LnHat[1] = 0.;
	LnHat[2] = cos(inc);

	/*
	 * Set the radial direction - need to take care to avoid singularity
	 * if L is along z axis
	 */
	if (LnHat[2] > 0.9999) {
		rHat[0] = 1.;
		rHat[1] = rHat[2] = 0.;
	} else {
		REAL8		theta0 = atan(-LnHat[2] / LnHat[0]);	/* theta0 is between 0
									 * and Pi */
		rHat[0] = sin(theta0);
		rHat[1] = 0;
		rHat[2] = cos(theta0);
	}

	/* Now we can complete the triad */
	vHat[0] = CalculateCrossProduct(0, LnHat, rHat);
	vHat[1] = CalculateCrossProduct(1, LnHat, rHat);
	vHat[2] = CalculateCrossProduct(2, LnHat, rHat);

	NormalizeVector(vHat);

	/* Vectors BEFORE rotation */
#if 0
		for (i = 0; i < 3; i++)
			print_debug(" LnHat[%d] = %.16e, rHat[%d] = %.16e, vHat[%d] = %.16e\n",
        i, LnHat[i], i, rHat[i], i, vHat[i]);

		print_debug("\n\n");
		for (i = 0; i < 3; i++)
			print_debug(" s1[%d] = %.16e, s2[%d] = %.16e\n", i, tmpS1[i], i, tmpS2[i]);
    //fflush(NULL);
#endif

	/* Allocate and compute the rotation matrices */
	rotMatrix = gsl_matrix_alloc(3, 3);
	invMatrix = gsl_matrix_alloc(3, 3);
	if (!rotMatrix || !invMatrix) 
    {
		if (rotMatrix)
			gsl_matrix_free(rotMatrix);
		if (invMatrix)
			gsl_matrix_free(invMatrix);
		return CEV_FAILURE;
	}
	if (CalculateRotationMatrix(rotMatrix, invMatrix, rHat, vHat, LnHat) == CEV_FAILURE) {
		gsl_matrix_free(rotMatrix);
		gsl_matrix_free(invMatrix);
		return CEV_FAILURE;
	}
	/* Rotate the orbital vectors and spins */
	ApplyRotationMatrix(rotMatrix, rHat);
	ApplyRotationMatrix(rotMatrix, vHat);
	ApplyRotationMatrix(rotMatrix, LnHat);
	ApplyRotationMatrix(rotMatrix, tmpS1);
	ApplyRotationMatrix(rotMatrix, tmpS2);
	ApplyRotationMatrix(rotMatrix, tmpS1Norm);
	ApplyRotationMatrix(rotMatrix, tmpS2Norm);

	/* See if Vectors have been rotated fine */
#if 0
		print_debug("\nAfter applying rotation matrix:\n\n");
		for (i = 0; i < 3; i++)
			print_debug(" LnHat[%d] = %.16e, rHat[%d] = %.16e, vHat[%d] = %.16e\n",
                i, LnHat[i], i, rHat[i], i, vHat[i]);

		print_debug("\n");
		for (i = 0; i < 3; i++)
			print_debug(" s1[%d] = %.16e, s2[%d] = %.16e\n", i, tmpS1[i], i, tmpS2[i]);

    //fflush(NULL);
#endif
	/*
	 * STEP 2) After rotation in STEP 1, in spherical coordinates, phi0
	 * and theta0 are given directly in Eq. (4.7), r0, pr0, ptheta0 and
	 * pphi0 are obtained by solving Eqs. (4.8) and (4.9) (using
	 * gsl_multiroot_fsolver). At this step, we find initial conditions
	 * for a spherical orbit without radiation reaction.
	 */

  /* Initialise the gsl stuff */
	rootSolver = gsl_multiroot_fsolver_alloc(T, 3);
	if (!rootSolver) {
		gsl_matrix_free(rotMatrix);
		gsl_matrix_free(invMatrix);
		return CEV_FAILURE;
	}
	initValues = gsl_vector_calloc(3);
	if (!initValues) {
		gsl_multiroot_fsolver_free(rootSolver);
		gsl_matrix_free(rotMatrix);
		gsl_matrix_free(invMatrix);
		return CEV_FAILURE;
	}

	rootFunction.f = XLALFindSphericalOrbitPrec;
	rootFunction.n = 3;
	rootFunction.params = &rootParams;

    /* Set to use optimized or unoptimized code */
    // rootParams.use_optimized = 0;

	/* Calculate the initial velocity from the given initial frequency */
    REAL8 fMinE = fMin;
    rootParams.e0 = 0.0;
    if (CODE_VERSION == 1)
        rootParams.e0 = e0;
    else if (CODE_VERSION == 2)
    {
        fMinE /= pow(1-e0*e0, 1.5);
    }
	omega = CST_PI * mTotal * CST_MTSUN_SI * fMinE;
	v0 = cbrt(omega);

	/* Given this, we can start to calculate the initial conditions */
	/* for spherical coords in the new basis */
	rootParams.omega = omega;
	rootParams.params = params;
	/* To start with, we will just assign Newtonian-ish ICs to the system */


	rootParams.values[0] = scale1 * sqrt((1.-e0)*(1.-e0) / (v0*v0*v0*v0) - 36.0);	/* Initial r */
	rootParams.values[4] = scale2 * pow(1.-e0, 0.5) * v0;	            /* Initial p */
	rootParams.values[5] = scale3 * 1e-3;
	//PK
    memcpy(rootParams.values + 6, tmpS1, sizeof(tmpS1));
	memcpy(rootParams.values + 9, tmpS2, sizeof(tmpS2));

#if 0
    print_debug("omega = %.16e\n", omega);
    print_debug("ICs guess: x = %.16e, py = %.16e, pz = %.16e\n",
      rootParams.values[0]/scale1, rootParams.values[4]/scale2,
      rootParams.values[5]/scale3);
    //fflush(NULL);
#endif
	gsl_vector_set(initValues, 0, rootParams.values[0]);
	gsl_vector_set(initValues, 1, rootParams.values[4]);
    gsl_vector_set(initValues, 2, rootParams.values[5]);

	gsl_multiroot_fsolver_set(rootSolver, &rootFunction, initValues);

	/* We are now ready to iterate to find the solution */
	i = 0;

    INT4 jittered=0;
	REAL8 r_now = 0.0;
	REAL8 temp = 0.0;
	do 
    {
            gslStatus = gsl_multiroot_fsolver_iterate(rootSolver);
            if (gslStatus == GSL_ENOPROG || gslStatus == GSL_ENOPROGJ) 
            {
                PRINT_LOG_INFO(LOG_ERROR, "NO PROGRESS being made by Spherical orbit root solver\n");

                /* Print Residual Function values whose roots we are trying to find */
                finalValues = gsl_multiroot_fsolver_f(rootSolver);
                PRINT_LOG_INFO(LOG_ERROR, "Function value here given by the following:\n");
                PRINT_LOG_INFO(LOG_ERROR, " F1 = %.16e, F2 = %.16e, F3 = %.16e\n",
                gsl_vector_get(finalValues, 0),
                gsl_vector_get(finalValues, 1), gsl_vector_get(finalValues, 2));

                /* Print Step sizes in each of function variables */
                finalValues = gsl_multiroot_fsolver_dx(rootSolver);
                // XLAL_PRINT_INFO("Stepsizes in each dimension:\n");
                // XLAL_PRINT_INFO(" x = %.16e, py = %.16e, pz = %.16e\n",
                // gsl_vector_get(finalValues, 0)/scale1,
                // gsl_vector_get(finalValues, 1)/scale2,
                // gsl_vector_get(finalValues, 2)/scale3);

                /* Only allow this flag to be caught MAXcntGslNoProgress no. of times */
                cntGslNoProgress += 1;
                if (cntGslNoProgress >= MAXcntGslNoProgress) 
                {
                    cntGslNoProgress = 0;

                    if(multFacGslNoProgress < 1.) { multFacGslNoProgress *= 1.02; }
                    else { multFacGslNoProgress /= 1.01; }

                }
                /* Now that no progress is being made, we need to reset the initial guess
                * for the (r,pPhi, pTheta) and reset the integrator */
                rootParams.values[0] = scale1 * sqrt(1. / (v0 * v0)*1./(v0*v0) -36.0);	/* Initial r */
                rootParams.values[4] = scale2 * v0;	            /* Initial p */
                if( cntGslNoProgress % 2 )
                    rootParams.values[5] = scale3 * 1e-3 / multFacGslNoProgress;
                else
                    rootParams.values[5] = scale3 * 1e-3 * multFacGslNoProgress;
                //PK
                memcpy(rootParams.values + 6, tmpS1, sizeof(tmpS1));
                memcpy(rootParams.values + 9, tmpS2, sizeof(tmpS2));
                // XLAL_PRINT_INFO("New ICs guess: x = %.16e, py = %.16e, pz = %.16e\n",
                //         rootParams.values[0]/scale1, rootParams.values[4]/scale2,
                //         rootParams.values[5]/scale3);
                // fflush(NULL);
                gsl_vector_set(initValues, 0, rootParams.values[0]);
                gsl_vector_set(initValues, 1, rootParams.values[4]);
                gsl_vector_set(initValues, 2, rootParams.values[5]);
                gsl_multiroot_fsolver_set(rootSolver, &rootFunction, initValues);
            }
            else if (gslStatus == GSL_EBADFUNC) 
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "Inf or Nan encountered in evaluluation of spherical orbit Eqn");
                gsl_multiroot_fsolver_free(rootSolver);
                gsl_vector_free(initValues);
                gsl_matrix_free(rotMatrix);
                gsl_matrix_free(invMatrix);
                return CEV_FAILURE;
            }
            else if (gslStatus != GSL_SUCCESS) 
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "Error in GSL iteration function!");
                gsl_multiroot_fsolver_free(rootSolver);
                gsl_vector_free(initValues);
                gsl_matrix_free(rotMatrix);
                gsl_matrix_free(invMatrix);
                return CEV_FAILURE;
            }

        /* different ways to test convergence of the method */
		gslStatus = gsl_multiroot_test_residual(rootSolver->f, 1.0e-9);
        /*XLAL_CALLGSL(gslStatus= gsl_multiroot_test_delta(
          gsl_multiroot_fsolver_dx(rootSolver),
          gsl_multiroot_fsolver_root(rootSolver),
          1.e-8, 1.e-5));*/

		if (jittered==0) 
        {
            finalValues = gsl_multiroot_fsolver_dx(rootSolver);
            if (isnan(gsl_vector_get(finalValues, 1))) 
            {
                rootParams.values[0] = scale1 * sqrt(1. / (v0 * v0)*1/(v0*v0)*(1.+1.e-8) - 36.0);	/* Jitter on initial r */
                rootParams.values[4] = scale2 * v0*(1.-1.e-8);	            /* Jitter on initial p */
                rootParams.values[5] = scale3 * 1e-3;
                memcpy(rootParams.values + 6, tmpS1, sizeof(tmpS1));
                memcpy(rootParams.values + 9, tmpS2, sizeof(tmpS2));
                gsl_vector_set(initValues, 0, rootParams.values[0]);
                gsl_vector_set(initValues, 1, rootParams.values[4]);
                gsl_vector_set(initValues, 2, rootParams.values[5]);
                gsl_multiroot_fsolver_set(rootSolver, &rootFunction, initValues);
                jittered=1;
            }
		}
		i++;
	}
	while (gslStatus == GSL_CONTINUE && i <= maxIter);

    // if(debugPK) { fflush(NULL); fclose(out); }

	if (i > maxIter && gslStatus != GSL_SUCCESS) 
    {
		gsl_multiroot_fsolver_free(rootSolver);
		gsl_vector_free(initValues);
		gsl_matrix_free(rotMatrix);
		gsl_matrix_free(invMatrix);
		//XLAL_ERROR(XLAL_EMAXITER);
		return CEV_FAILURE;
	}
	finalValues = gsl_multiroot_fsolver_root(rootSolver);
#if 0
    print_debug("Spherical orbit conditions here given by the following:\n");
    print_debug(" x = %.16e, py = %.16e, pz = %.16e\n",
    gsl_vector_get(finalValues, 0)/scale1,
    gsl_vector_get(finalValues, 1)/scale2,
    gsl_vector_get(finalValues, 2)/scale3);
#endif
	memset(qCart, 0, sizeof(qCart));
	memset(pCart, 0, sizeof(pCart));

	qCart[0] = sqrt(gsl_vector_get(finalValues, 0)*gsl_vector_get(finalValues, 0)+36.0);
	pCart[1] = gsl_vector_get(finalValues, 1)/scale2;
	pCart[2] = gsl_vector_get(finalValues, 2)/scale3;
    if (CODE_VERSION == 2)
    {
        qCart[0] /= 1. + e0;
        pCart[1] *= 1. + e0;
    }

	/* Free the GSL root finder, since we're done with it */
	gsl_multiroot_fsolver_free(rootSolver);
	gsl_vector_free(initValues);


	/*
	 * STEP 3) Rotate to L0 along z-axis and N0 along x-axis frame, where
	 * L0 is the initial orbital angular momentum and L0 is calculated
	 * using initial position and linear momentum obtained in STEP 2.
	 */

	/* Now we can calculate the relativistic L and rotate to a new basis */
	memcpy(qHat, qCart, sizeof(qCart));
	memcpy(pHat, pCart, sizeof(pCart));

	NormalizeVector(qHat);
	NormalizeVector(pHat);

	Lhat[0] = CalculateCrossProduct(0, qHat, pHat);
	Lhat[1] = CalculateCrossProduct(1, qHat, pHat);
	Lhat[2] = CalculateCrossProduct(2, qHat, pHat);

	NormalizeVector(Lhat);

	rotMatrix2 = gsl_matrix_alloc(3, 3);
	invMatrix2 = gsl_matrix_alloc(3, 3);

	if (CalculateRotationMatrix(rotMatrix2, invMatrix2, qHat, pHat, Lhat) == CEV_FAILURE) 
    {
		gsl_matrix_free(rotMatrix);
		gsl_matrix_free(invMatrix);
		return CEV_FAILURE;
	}
	ApplyRotationMatrix(rotMatrix2, rHat);
	ApplyRotationMatrix(rotMatrix2, vHat);
	ApplyRotationMatrix(rotMatrix2, LnHat);
	ApplyRotationMatrix(rotMatrix2, tmpS1);
	ApplyRotationMatrix(rotMatrix2, tmpS2);
	ApplyRotationMatrix(rotMatrix2, tmpS1Norm);
	ApplyRotationMatrix(rotMatrix2, tmpS2Norm);
	ApplyRotationMatrix(rotMatrix2, qCart);
	ApplyRotationMatrix(rotMatrix2, pCart);

    gsl_matrix_free(rotMatrix);
    gsl_matrix_free(rotMatrix2);

    // XLAL_PRINT_INFO("qCart after rotation2 %3.10f %3.10f %3.10f\n", qCart[0], qCart[1], qCart[2]);
    // XLAL_PRINT_INFO("pCart after rotation2 %3.10f %3.10f %3.10f\n", pCart[0], pCart[1], pCart[2]);
    // XLAL_PRINT_INFO("S1 after rotation2 %3.10f %3.10f %3.10f\n", tmpS1Norm[0], tmpS1Norm[1], tmpS1Norm[2]);
    // XLAL_PRINT_INFO("S2 after rotation2 %3.10f %3.10f %3.10f\n", tmpS2Norm[0], tmpS2Norm[1], tmpS2Norm[2]);
	/*
	 * STEP 4) In the L0-N0 frame, we calculate (dE/dr)|sph using Eq.
	 * (4.14), then initial dr/dt using Eq. (4.10), and finally pr0 using
	 * Eq. (4.15).
	 */

	/* Now we can calculate the flux. Change to spherical co-ords */
	CartesianToSpherical(qSph, pSph, qCart, pCart);
	memcpy(sphValues, qSph, sizeof(qSph));
	memcpy(sphValues + 3, pSph, sizeof(pSph));
	memcpy(sphValues + 6, tmpS1, sizeof(tmpS1));
	memcpy(sphValues + 9, tmpS2, sizeof(tmpS2));

	memcpy(cartValues, qCart, sizeof(qCart));
	memcpy(cartValues + 3, pCart, sizeof(pCart));
	memcpy(cartValues + 6, tmpS1, sizeof(tmpS1));
	memcpy(cartValues + 9, tmpS2, sizeof(tmpS2));

	REAL8		dHdpphi , d2Hdr2, d2Hdrdpphi;
	REAL8		rDot    , dHdpr, flux, dEdr;

	d2Hdr2 = XLALCalculateSphHamiltonianDeriv2Prec(0, 0, sphValues, params);
	d2Hdrdpphi = XLALCalculateSphHamiltonianDeriv2Prec(0, 5, sphValues, params);

    // XLAL_PRINT_INFO("d2Hdr2 = %.16e, d2Hdrdpphi = %.16e\n", d2Hdr2, d2Hdrdpphi);

	/* New code to compute derivatives w.r.t. cartesian variables */

	REAL8		tmpDValues[14];
	int status;
	for (i = 0; i < 3; i++) {
		cartValues[i + 6] /= mTotal * mTotal;
		cartValues[i + 9] /= mTotal * mTotal;
	}
    UINT oldignoreflux = params->ignoreflux;
    params->ignoreflux = 1;
    status = XLALSpinPrecHcapNumericalDerivative(0, cartValues, tmpDValues, params);
    params->ignoreflux = oldignoreflux;
	for (i = 0; i < 3; i++) {
		cartValues[i + 6] *= mTotal * mTotal;
		cartValues[i + 9] *= mTotal * mTotal;
	}

	dHdpphi = tmpDValues[1] / sqrt(cartValues[0] * cartValues[0] + cartValues[1] * cartValues[1] + cartValues[2] * cartValues[2]);
	//XLALSpinPrecHcapNumDerivWRTParam(4, cartValues, params) / sphValues[0];

	dEdr = -dHdpphi * d2Hdr2 / d2Hdrdpphi;

    // XLAL_PRINT_INFO("d2Hdr2 = %.16e d2Hdrdpphi = %.16e dHdpphi = %.16e\n",
    //     d2Hdr2, d2Hdrdpphi, dHdpphi);

	if (d2Hdr2 != 0.0 && e0 == 0.0) 
    {
		/* We will need to calculate the Hamiltonian to get the flux */
		s1Vec.length = s2Vec.length = s1VecNorm.length = s2VecNorm.length = sKerr.length = sStar.length = 3;
		s1Vec.data = tmpS1;
		s2Vec.data = tmpS2;
		s1VecNorm.data = tmpS1Norm;
		s2VecNorm.data = tmpS2Norm;
		sKerr.data = sKerrData;
		sStar.data = sStarData;

		qCartVec.length = pCartVec.length = 3;
		qCartVec.data = qCart;
		pCartVec.data = pCart;

		//chi1 = tmpS1[0] * LnHat[0] + tmpS1[1] * LnHat[1] + tmpS1[2] * LnHat[2];
		//chi2 = tmpS2[0] * LnHat[0] + tmpS2[1] * LnHat[1] + tmpS2[2] * LnHat[2];

		//if (debugPK)
			//XLAL_PRINT_INFO("magS1 = %.16e, magS2 = %.16e\n", chi1, chi2);

		//chiS = 0.5 * (chi1 / (mass1 * mass1) + chi2 / (mass2 * mass2));
		//chiA = 0.5 * (chi1 / (mass1 * mass1) - chi2 / (mass2 * mass2));

		EOBCalculateSigmaKerr(&sKerr, &s1VecNorm, &s2VecNorm);
		EOBCalculateSigmaStar(&sStar, mass1, mass2, &s1VecNorm, &s2VecNorm);

		/*
		 * The a in the flux has been set to zero, but not in the
		 * Hamiltonian
		 */
		a = sqrt(sKerr.data[0] * sKerr.data[0] + sKerr.data[1] * sKerr.data[1] + sKerr.data[2] * sKerr.data[2]);
		//XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients(params->eobParams->hCoeffs, mass1, mass2, eta, /* a */ 0.0, chiS, chiA);
		//XLALSimIMRCalculateSpinPrecEOBHCoeffs(params->seobCoeffs, eta, a);
		ham = XLALSimIMRSpinPrecEOBHamiltonian(eta, &qCartVec, &pCartVec, &s1VecNorm, &s2VecNorm, &sKerr, &sStar, params->tortoise, params->seobCoeffs, params->hParams);

        // XLAL_PRINT_INFO("Stas: hamiltonian in ICs at this point is %.16e\n", ham);

		/* And now, finally, the flux */
		REAL8Vector	polarDynamics, cartDynamics;
		REAL8		polarData[4], cartData[12];

		polarDynamics.length = 4;
		polarDynamics.data = polarData;

		polarData[0] = qSph[0];
		polarData[1] = 0.;
		polarData[2] = pSph[0];
		polarData[3] = pSph[2];

		cartDynamics.length = 12;
		cartDynamics.data = cartData;

		memcpy(cartData, qCart, 3 * sizeof(REAL8));
		memcpy(cartData + 3, pCart, 3 * sizeof(REAL8));
		memcpy(cartData + 6, tmpS1Norm, 3 * sizeof(REAL8));
		memcpy(cartData + 9, tmpS2Norm, 3 * sizeof(REAL8));

		//XLAL_PRINT_INFO("Stas: starting FLux calculations\n");
        flux = XLALInspiralPrecSpinFactorizedFlux(&polarDynamics, &cartDynamics, nqcCoeffs, omega, 0, qSph[0]*omega, params, ham, lMax, SpinAlignedEOBversion);
		/*
		 * flux  = XLALInspiralSpinFactorizedFlux( &polarDynamics,
		 * nqcCoeffs, omega, params, ham, lMax, SpinAlignedEOBversion
		 * );
		 */
		PRINT_LOG_INFO(LOG_DEBUG, "Initial Conditions: Stas flux = %.16e", flux);
		//exit(0);
		flux = flux / eta;
        if (e0 != 0.0)
            rDot = 0.0;
        else
		    rDot = -flux / dEdr;
		/*
		 * We now need dHdpr - we take it that it is safely linear up
		 * to a pr of 1.0e-3 PK: Ideally, the pr should be of the
		 * order of other momenta, in order for its contribution to
		 * the Hamiltonian to not get buried in the numerical noise
		 * in the numerically larger momenta components
		 */
		cartValues[3] = 1.0e-3;
		for (i = 0; i < 3; i++) 
        {
			cartValues[i + 6] /= mTotal * mTotal;
			cartValues[i + 9] /= mTotal * mTotal;
		}
        oldignoreflux = params->ignoreflux;
        params->ignoreflux = 1;
        params->seobCoeffs->updateHCoeffs = 1;
	    status = XLALSpinPrecHcapNumericalDerivative(0, cartValues, tmpDValues, params);
        params->ignoreflux = oldignoreflux;
		for (i = 0; i < 3; i++) 
        {
			cartValues[i + 6] *= mTotal * mTotal;
			cartValues[i + 9] *= mTotal * mTotal;
		}
        REAL8   csi = sqrt(XLALSimIMRSpinPrecEOBHamiltonianDeltaT(params->seobCoeffs, qSph[0], eta, a)*XLALSimIMRSpinPrecEOBHamiltonianDeltaR(params->seobCoeffs, qSph[0], eta, a)) / (qSph[0] * qSph[0] + a * a);

		dHdpr = csi*tmpDValues[0];
		//XLALSpinPrecHcapNumDerivWRTParam(3, cartValues, params);

        // XLAL_PRINT_INFO("Ingredients going into prDot:\n");
        // XLAL_PRINT_INFO("flux = %.16e, dEdr = %.16e, dHdpr = %.16e, dHdpr/pr = %.16e\n", flux, dEdr, dHdpr, dHdpr / cartValues[3]);
		/*
		 * We can now calculate what pr should be taking into account
		 * the flux
		 */
		pSph[0] = rDot / (dHdpr / cartValues[3]);
	} else 
    {
		/*
		 * Since d2Hdr2 has evaluated to zero, we cannot do the
		 * above. Just set pr to zero
		 */
		//XLAL_PRINT_INFO("d2Hdr2 is zero!\n");
		pSph[0] = 0;
	}

	/* Now we are done - convert back to cartesian coordinates ) */
	SphericalToCartesian(qCart, pCart, qSph, pSph);
    PRINT_LOG_INFO(LOG_DEBUG, "Sph initial condition : r = (%e,%e,%e), p = (%e,%e,%e)", qSph[0], qSph[1], qSph[2], pSph[0], pSph[1], pSph[2]);

	/*
	 * STEP 5) Rotate back to the original inertial frame by inverting
	 * the rotation of STEP 3 and then  inverting the rotation of STEP 1.
	 */

	/* Undo rotations to get back to the original basis */
	/* Second rotation */
	ApplyRotationMatrix(invMatrix2, rHat);
	ApplyRotationMatrix(invMatrix2, vHat);
	ApplyRotationMatrix(invMatrix2, LnHat);
	ApplyRotationMatrix(invMatrix2, tmpS1);
	ApplyRotationMatrix(invMatrix2, tmpS2);
	ApplyRotationMatrix(invMatrix2, tmpS1Norm);
	ApplyRotationMatrix(invMatrix2, tmpS2Norm);
	ApplyRotationMatrix(invMatrix2, qCart);
	ApplyRotationMatrix(invMatrix2, pCart);

	/* First rotation */
	ApplyRotationMatrix(invMatrix, rHat);
	ApplyRotationMatrix(invMatrix, vHat);
	ApplyRotationMatrix(invMatrix, LnHat);
	ApplyRotationMatrix(invMatrix, tmpS1);
	ApplyRotationMatrix(invMatrix, tmpS2);
	ApplyRotationMatrix(invMatrix, tmpS1Norm);
	ApplyRotationMatrix(invMatrix, tmpS2Norm);
	ApplyRotationMatrix(invMatrix, qCart);
	ApplyRotationMatrix(invMatrix, pCart);

    gsl_matrix_free(invMatrix);
    gsl_matrix_free(invMatrix2);

        /* If required, apply the tortoise transform */
	if (tmpTortoise) 
    {
		REAL8		r = sqrt(qCart[0] * qCart[0] + qCart[1] * qCart[1] + qCart[2] * qCart[2]);
		REAL8		deltaR = XLALSimIMRSpinPrecEOBHamiltonianDeltaR(params->seobCoeffs, r, eta, a);
		REAL8		deltaT = XLALSimIMRSpinPrecEOBHamiltonianDeltaT(params->seobCoeffs, r, eta, a);
		REAL8		csi = sqrt(deltaT * deltaR) / (r * r + a * a);

		REAL8		pr = (qCart[0] * pCart[0] + qCart[1] * pCart[1] + qCart[2] * pCart[2]) / r;

		params->tortoise = tmpTortoise;

        // XLAL_PRINT_INFO("Applying the tortoise to p (csi = %.26e)\n", csi);
        // XLAL_PRINT_INFO("pCart = %3.10f %3.10f %3.10f\n", pCart[0], pCart[1], pCart[2]);
		for (i = 0; i < 3; i++) 
        {
			pCart[i] = pCart[i] + qCart[i] * pr * (csi - 1.) / r;
		}
	}


    /* Now copy the initial conditions back to the return vector */
	memcpy(initConds->data, qCart, sizeof(qCart));
	memcpy(initConds->data + 3, pCart, sizeof(pCart));
	memcpy(initConds->data + 6, tmpS1Norm, sizeof(tmpS1Norm));
	memcpy(initConds->data + 9, tmpS2Norm, sizeof(tmpS2Norm));
    for (i=0; i<12; i++) 
    {
        if (fabs(initConds->data[i]) <=1.0e-15) 
        {
            initConds->data[i] = 0.;
        }
    }
    // REAL8 initr = sqrt(qCart[0]*qCart[0]+qCart[1]*qCart[1]+qCart[2]*qCart[2]);
    // if ( initr< 3.)
    // {
    //     PRINT_LOG_INFO(LOG_CRITICAL, "the initial polar r = %.8e is too small, abort.", initr);
    //     return CEV_FAILURE;
    // }

    // XLAL_PRINT_INFO("THE FINAL INITIAL CONDITIONS:\n");
    // XLAL_PRINT_INFO(" %.16e %.16e %.16e\n%.16e %.16e %.16e\n%.16e %.16e %.16e\n%.16e %.16e %.16e\n", initConds->data[0], initConds->data[1], initConds->data[2],
    //         initConds->data[3], initConds->data[4], initConds->data[5], initConds->data[6], initConds->data[7], initConds->data[8],
    //         initConds->data[9], initConds->data[10], initConds->data[11]);
    return CEV_SUCCESS;
}


INT
EOBInitialConditionsPrec_Conserve(
                    REAL8Vector * initConds,	/**<< OUTPUT, Initial dynamical variables */
                    const REAL8 mass1,	/**<< mass 1 */
                    const REAL8 mass2,	/**<< mass 2 */
                    const REAL8 fMin,	/**<< Initial frequency (given) */
                    const REAL8 e0,
                    const REAL8 inc,	/**<< Inclination */
                    const REAL8 spin1[],	/**<< Initial spin vector 1 */
                    const REAL8 spin2[],	/**<< Initial spin vector 2 */
                    SpinEOBParams * params	/**<< Spin EOB parameters */
)
{
	// 	XLAL_PRINT_INFO("Inside the XLALSimIMRSpinEOBInitialConditionsPrec function!\n");
	// 	XLAL_PRINT_INFO(
    // "Inputs: m1 = %.16e, m2 = %.16e, fMin = %.16e, inclination = %.16e\n",
    //   mass1, mass2, (double)fMin, (double)inc);
	// 	XLAL_PRINT_INFO("Inputs: mSpin1 = {%.16e, %.16e, %.16e}\n",
    //   spin1[0], spin1[1], spin1[2]);
	// 	XLAL_PRINT_INFO("Inputs: mSpin2 = {%.16e, %.16e, %.16e}\n",
    //   spin2[0], spin2[1], spin2[2]);
	// 	fflush(NULL);
	static const int lMax = 8;

	int		i;

	/* Variable to keep track of whether the user requested the tortoise */
	int		tmpTortoise;

	UINT		SpinAlignedEOBversion;

	REAL8		mTotal;
	REAL8		eta;
	REAL8		omega   , v0;	/* Initial velocity and angular
					 * frequency */

	REAL8		ham;	/* Hamiltonian */

	REAL8		LnHat    [3];	/* Initial orientation of angular
					 * momentum */
	REAL8		rHat     [3];	/* Initial orientation of radial
					 * vector */
	REAL8		vHat     [3];	/* Initial orientation of velocity
					 * vector */
	REAL8		Lhat     [3];	/* Direction of relativistic ang mom */
	REAL8		qHat     [3];
	REAL8		pHat     [3];

	/* q and p vectors in Cartesian and spherical coords */
	REAL8		qCart    [3], pCart[3];
	REAL8		qSph     [3], pSph[3];

	/* We will need to manipulate the spin vectors */
	/* We will use temporary vectors to do this */
	REAL8		tmpS1    [3];
	REAL8		tmpS2    [3];
	REAL8		tmpS1Norm[3];
	REAL8		tmpS2Norm[3];

	REAL8Vector	qCartVec, pCartVec;
	REAL8Vector	s1Vec, s2Vec, s1VecNorm, s2VecNorm;
	REAL8Vector	sKerr, sStar;
	REAL8		sKerrData[3], sStarData[3];
	REAL8		a = 0.;
	//, chiS, chiA;
	//REAL8 chi1, chi2;

	/*
	 * We will need a full values vector for calculating derivs of
	 * Hamiltonian
	 */
	REAL8		sphValues[12];
	REAL8		cartValues[12];

	/* Matrices for rotating to the new basis set. */
	/* It is more convenient to calculate the ICs in a simpler basis */
	gsl_matrix     *rotMatrix = NULL;
	gsl_matrix     *invMatrix = NULL;
	gsl_matrix     *rotMatrix2 = NULL;
	gsl_matrix     *invMatrix2 = NULL;

	/* Root finding stuff for finding the spherical orbit */
	SEOBRootParams	rootParams;
	const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrids;
	gsl_multiroot_fsolver *rootSolver = NULL;

	gsl_multiroot_function rootFunction;
	gsl_vector     *initValues = NULL;
	gsl_vector     *finalValues = NULL;
	INT gslStatus;
    INT cntGslNoProgress = 0, MAXcntGslNoProgress = 5;
    //INT cntGslNoProgress = 0, MAXcntGslNoProgress = 50;
    REAL8 multFacGslNoProgress = 3./5.;
	//const int	maxIter = 2000;
	const int	maxIter = 10000;

	memset(&rootParams, 0, sizeof(rootParams));

	mTotal = mass1 + mass2;
	eta = mass1 * mass2 / (mTotal * mTotal);
	memcpy(tmpS1, spin1, sizeof(tmpS1));
	memcpy(tmpS2, spin2, sizeof(tmpS2));
	memcpy(tmpS1Norm, spin1, sizeof(tmpS1Norm));
	memcpy(tmpS2Norm, spin2, sizeof(tmpS2Norm));
	for (i = 0; i < 3; i++) {
		tmpS1Norm[i] /= mTotal * mTotal;
		tmpS2Norm[i] /= mTotal * mTotal;
	}
	SpinAlignedEOBversion = 4;
	/* We compute the ICs for the non-tortoise p, and convert at the end */
	tmpTortoise = params->tortoise;
	params->tortoise = 0;

	EOBNonQCCoeffs *nqcCoeffs = NULL;
	nqcCoeffs = params->nqcCoeffs;

	/*
	 * STEP 1) Rotate to LNhat0 along z-axis and N0 along x-axis frame,
	 * where LNhat0 and N0 are initial normal to orbital plane and
	 * initial orbital separation;
	 */

	/* Set the initial orbital ang mom direction. Taken from STPN code */
	LnHat[0] = sin(inc);
	LnHat[1] = 0.;
	LnHat[2] = cos(inc);

	/*
	 * Set the radial direction - need to take care to avoid singularity
	 * if L is along z axis
	 */
	if (LnHat[2] > 0.9999) {
		rHat[0] = 1.;
		rHat[1] = rHat[2] = 0.;
	} else {
		REAL8		theta0 = atan(-LnHat[2] / LnHat[0]);	/* theta0 is between 0
									 * and Pi */
		rHat[0] = sin(theta0);
		rHat[1] = 0;
		rHat[2] = cos(theta0);
	}

	/* Now we can complete the triad */
	vHat[0] = CalculateCrossProduct(0, LnHat, rHat);
	vHat[1] = CalculateCrossProduct(1, LnHat, rHat);
	vHat[2] = CalculateCrossProduct(2, LnHat, rHat);

	NormalizeVector(vHat);

	/* Vectors BEFORE rotation */
#if 0
		for (i = 0; i < 3; i++)
			XLAL_PRINT_INFO(" LnHat[%d] = %.16e, rHat[%d] = %.16e, vHat[%d] = %.16e\n",
        i, LnHat[i], i, rHat[i], i, vHat[i]);

		XLAL_PRINT_INFO("\n\n");
		for (i = 0; i < 3; i++)
			XLAL_PRINT_INFO(" s1[%d] = %.16e, s2[%d] = %.16e\n", i, tmpS1[i], i, tmpS2[i]);
    fflush(NULL);
#endif

	/* Allocate and compute the rotation matrices */
	rotMatrix = gsl_matrix_alloc(3, 3);
	invMatrix = gsl_matrix_alloc(3, 3);
	if (!rotMatrix || !invMatrix) 
    {
		if (rotMatrix)
			gsl_matrix_free(rotMatrix);
		if (invMatrix)
			gsl_matrix_free(invMatrix);
		return CEV_FAILURE;
	}
	if (CalculateRotationMatrix(rotMatrix, invMatrix, rHat, vHat, LnHat) == CEV_FAILURE) {
		gsl_matrix_free(rotMatrix);
		gsl_matrix_free(invMatrix);
		return CEV_FAILURE;
	}
	/* Rotate the orbital vectors and spins */
	ApplyRotationMatrix(rotMatrix, rHat);
	ApplyRotationMatrix(rotMatrix, vHat);
	ApplyRotationMatrix(rotMatrix, LnHat);
	ApplyRotationMatrix(rotMatrix, tmpS1);
	ApplyRotationMatrix(rotMatrix, tmpS2);
	ApplyRotationMatrix(rotMatrix, tmpS1Norm);
	ApplyRotationMatrix(rotMatrix, tmpS2Norm);

	/* See if Vectors have been rotated fine */
#if 0
		XLAL_PRINT_INFO("\nAfter applying rotation matrix:\n\n");
		for (i = 0; i < 3; i++)
			XLAL_PRINT_INFO(" LnHat[%d] = %.16e, rHat[%d] = %.16e, vHat[%d] = %.16e\n",
                i, LnHat[i], i, rHat[i], i, vHat[i]);

		XLAL_PRINT_INFO("\n");
		for (i = 0; i < 3; i++)
			XLAL_PRINT_INFO(" s1[%d] = %.16e, s2[%d] = %.16e\n", i, tmpS1[i], i, tmpS2[i]);

    fflush(NULL);
#endif
	/*
	 * STEP 2) After rotation in STEP 1, in spherical coordinates, phi0
	 * and theta0 are given directly in Eq. (4.7), r0, pr0, ptheta0 and
	 * pphi0 are obtained by solving Eqs. (4.8) and (4.9) (using
	 * gsl_multiroot_fsolver). At this step, we find initial conditions
	 * for a spherical orbit without radiation reaction.
	 */

  /* Initialise the gsl stuff */
	rootSolver = gsl_multiroot_fsolver_alloc(T, 3);
	if (!rootSolver) {
		gsl_matrix_free(rotMatrix);
		gsl_matrix_free(invMatrix);
		return CEV_FAILURE;
	}
	initValues = gsl_vector_calloc(3);
	if (!initValues) {
		gsl_multiroot_fsolver_free(rootSolver);
		gsl_matrix_free(rotMatrix);
		gsl_matrix_free(invMatrix);
		return CEV_FAILURE;
	}

	rootFunction.f = XLALFindSphericalOrbitPrec;
	rootFunction.n = 3;
	rootFunction.params = &rootParams;

    /* Set to use optimized or unoptimized code */
    // rootParams.use_optimized = 0;

	/* Calculate the initial velocity from the given initial frequency */
	omega = CST_PI * mTotal * CST_MTSUN_SI * fMin / pow(1.+e0*e0, 1.5);
	v0 = cbrt(omega);

	/* Given this, we can start to calculate the initial conditions */
	/* for spherical coords in the new basis */
	rootParams.omega = omega;
	rootParams.params = params;
    rootParams.e0 = e0;
	/* To start with, we will just assign Newtonian-ish ICs to the system */


	rootParams.values[0] = scale1 * sqrt((1.+e0) / (v0 * v0)*(1.+e0)/(v0*v0) - 36.0);	/* Initial r */
	rootParams.values[4] = scale2 * pow(1.+e0, 0.5) * v0;	            /* Initial p */
	rootParams.values[5] = scale3 * 1e-3;
	//PK
  memcpy(rootParams.values + 6, tmpS1, sizeof(tmpS1));
	memcpy(rootParams.values + 9, tmpS2, sizeof(tmpS2));

#if 0
    XLAL_PRINT_INFO("ICs guess: x = %.16e, py = %.16e, pz = %.16e\n",
      rootParams.values[0]/scale1, rootParams.values[4]/scale2,
      rootParams.values[5]/scale3);
    fflush(NULL);
#endif
	gsl_vector_set(initValues, 0, rootParams.values[0]);
	gsl_vector_set(initValues, 1, rootParams.values[4]);
    gsl_vector_set(initValues, 2, rootParams.values[5]);

	gsl_multiroot_fsolver_set(rootSolver, &rootFunction, initValues);

	/* We are now ready to iterate to find the solution */
	i = 0;

    INT4 jittered=0;
	REAL8 r_now = 0.0;
	REAL8 temp = 0.0;
	do 
    {
            gslStatus = gsl_multiroot_fsolver_iterate(rootSolver);
            if (gslStatus == GSL_ENOPROG || gslStatus == GSL_ENOPROGJ) 
            {
                PRINT_LOG_INFO(LOG_ERROR, "NO PROGRESS being made by Spherical orbit root solver\n");

                /* Print Residual Function values whose roots we are trying to find */
                finalValues = gsl_multiroot_fsolver_f(rootSolver);
                PRINT_LOG_INFO(LOG_ERROR, "Function value here given by the following:\n");
                PRINT_LOG_INFO(LOG_ERROR, " F1 = %.16e, F2 = %.16e, F3 = %.16e\n",
                gsl_vector_get(finalValues, 0),
                gsl_vector_get(finalValues, 1), gsl_vector_get(finalValues, 2));

                /* Print Step sizes in each of function variables */
                finalValues = gsl_multiroot_fsolver_dx(rootSolver);
                // XLAL_PRINT_INFO("Stepsizes in each dimension:\n");
                // XLAL_PRINT_INFO(" x = %.16e, py = %.16e, pz = %.16e\n",
                // gsl_vector_get(finalValues, 0)/scale1,
                // gsl_vector_get(finalValues, 1)/scale2,
                // gsl_vector_get(finalValues, 2)/scale3);

                /* Only allow this flag to be caught MAXcntGslNoProgress no. of times */
                cntGslNoProgress += 1;
                if (cntGslNoProgress >= MAXcntGslNoProgress) 
                {
                    cntGslNoProgress = 0;

                    if(multFacGslNoProgress < 1.) { multFacGslNoProgress *= 1.02; }
                    else { multFacGslNoProgress /= 1.01; }

                }
                /* Now that no progress is being made, we need to reset the initial guess
                * for the (r,pPhi, pTheta) and reset the integrator */
                rootParams.values[0] = scale1 * sqrt(1. / (v0 * v0)*1./(v0*v0) -36.0);	/* Initial r */
                rootParams.values[4] = scale2 * v0;	            /* Initial p */
                if( cntGslNoProgress % 2 )
                    rootParams.values[5] = scale3 * 1e-3 / multFacGslNoProgress;
                else
                    rootParams.values[5] = scale3 * 1e-3 * multFacGslNoProgress;
                //PK
                memcpy(rootParams.values + 6, tmpS1, sizeof(tmpS1));
                memcpy(rootParams.values + 9, tmpS2, sizeof(tmpS2));
                // XLAL_PRINT_INFO("New ICs guess: x = %.16e, py = %.16e, pz = %.16e\n",
                //         rootParams.values[0]/scale1, rootParams.values[4]/scale2,
                //         rootParams.values[5]/scale3);
                // fflush(NULL);
                gsl_vector_set(initValues, 0, rootParams.values[0]);
                gsl_vector_set(initValues, 1, rootParams.values[4]);
                gsl_vector_set(initValues, 2, rootParams.values[5]);
                gsl_multiroot_fsolver_set(rootSolver, &rootFunction, initValues);
            }
            else if (gslStatus == GSL_EBADFUNC) 
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "Inf or Nan encountered in evaluluation of spherical orbit Eqn");
                gsl_multiroot_fsolver_free(rootSolver);
                gsl_vector_free(initValues);
                gsl_matrix_free(rotMatrix);
                gsl_matrix_free(invMatrix);
                return CEV_FAILURE;
            }
            else if (gslStatus != GSL_SUCCESS) 
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "Error in GSL iteration function!");
                gsl_multiroot_fsolver_free(rootSolver);
                gsl_vector_free(initValues);
                gsl_matrix_free(rotMatrix);
                gsl_matrix_free(invMatrix);
                return CEV_FAILURE;
            }

        /* different ways to test convergence of the method */
		gslStatus = gsl_multiroot_test_residual(rootSolver->f, 1.0e-9);
        /*XLAL_CALLGSL(gslStatus= gsl_multiroot_test_delta(
          gsl_multiroot_fsolver_dx(rootSolver),
          gsl_multiroot_fsolver_root(rootSolver),
          1.e-8, 1.e-5));*/

		if (jittered==0) 
        {
            finalValues = gsl_multiroot_fsolver_dx(rootSolver);
            if (isnan(gsl_vector_get(finalValues, 1))) 
            {
                rootParams.values[0] = scale1 * sqrt(1. / (v0 * v0)*1/(v0*v0)*(1.+1.e-8) - 36.0);	/* Jitter on initial r */
                rootParams.values[4] = scale2 * v0*(1.-1.e-8);	            /* Jitter on initial p */
                rootParams.values[5] = scale3 * 1e-3;
                memcpy(rootParams.values + 6, tmpS1, sizeof(tmpS1));
                memcpy(rootParams.values + 9, tmpS2, sizeof(tmpS2));
                gsl_vector_set(initValues, 0, rootParams.values[0]);
                gsl_vector_set(initValues, 1, rootParams.values[4]);
                gsl_vector_set(initValues, 2, rootParams.values[5]);
                gsl_multiroot_fsolver_set(rootSolver, &rootFunction, initValues);
                jittered=1;
            }
		}
		i++;
	}
	while (gslStatus == GSL_CONTINUE && i <= maxIter);

    // if(debugPK) { fflush(NULL); fclose(out); }

	if (i > maxIter && gslStatus != GSL_SUCCESS) 
    {
		gsl_multiroot_fsolver_free(rootSolver);
		gsl_vector_free(initValues);
		gsl_matrix_free(rotMatrix);
		gsl_matrix_free(invMatrix);
		//XLAL_ERROR(XLAL_EMAXITER);
		return CEV_FAILURE;
	}
	finalValues = gsl_multiroot_fsolver_root(rootSolver);

    // XLAL_PRINT_INFO("Spherical orbit conditions here given by the following:\n");
    // XLAL_PRINT_INFO(" x = %.16e, py = %.16e, pz = %.16e\n",
    // gsl_vector_get(finalValues, 0)/scale1,
    // gsl_vector_get(finalValues, 1)/scale2,
    // gsl_vector_get(finalValues, 2)/scale3);
	memset(qCart, 0, sizeof(qCart));
	memset(pCart, 0, sizeof(pCart));

	qCart[0] = sqrt(gsl_vector_get(finalValues, 0)*gsl_vector_get(finalValues, 0)+36.0);
	pCart[1] = gsl_vector_get(finalValues, 1)/scale2;
	pCart[2] = gsl_vector_get(finalValues, 2)/scale3;


	/* Free the GSL root finder, since we're done with it */
	gsl_multiroot_fsolver_free(rootSolver);
	gsl_vector_free(initValues);


	/*
	 * STEP 3) Rotate to L0 along z-axis and N0 along x-axis frame, where
	 * L0 is the initial orbital angular momentum and L0 is calculated
	 * using initial position and linear momentum obtained in STEP 2.
	 */

	/* Now we can calculate the relativistic L and rotate to a new basis */
	memcpy(qHat, qCart, sizeof(qCart));
	memcpy(pHat, pCart, sizeof(pCart));

	NormalizeVector(qHat);
	NormalizeVector(pHat);

	Lhat[0] = CalculateCrossProduct(0, qHat, pHat);
	Lhat[1] = CalculateCrossProduct(1, qHat, pHat);
	Lhat[2] = CalculateCrossProduct(2, qHat, pHat);

	NormalizeVector(Lhat);

	rotMatrix2 = gsl_matrix_alloc(3, 3);
	invMatrix2 = gsl_matrix_alloc(3, 3);

	if (CalculateRotationMatrix(rotMatrix2, invMatrix2, qHat, pHat, Lhat) == CEV_FAILURE) 
    {
		gsl_matrix_free(rotMatrix);
		gsl_matrix_free(invMatrix);
		return CEV_FAILURE;
	}
	ApplyRotationMatrix(rotMatrix2, rHat);
	ApplyRotationMatrix(rotMatrix2, vHat);
	ApplyRotationMatrix(rotMatrix2, LnHat);
	ApplyRotationMatrix(rotMatrix2, tmpS1);
	ApplyRotationMatrix(rotMatrix2, tmpS2);
	ApplyRotationMatrix(rotMatrix2, tmpS1Norm);
	ApplyRotationMatrix(rotMatrix2, tmpS2Norm);
	ApplyRotationMatrix(rotMatrix2, qCart);
	ApplyRotationMatrix(rotMatrix2, pCart);

    gsl_matrix_free(rotMatrix);
    gsl_matrix_free(rotMatrix2);

    // XLAL_PRINT_INFO("qCart after rotation2 %3.10f %3.10f %3.10f\n", qCart[0], qCart[1], qCart[2]);
    // XLAL_PRINT_INFO("pCart after rotation2 %3.10f %3.10f %3.10f\n", pCart[0], pCart[1], pCart[2]);
    // XLAL_PRINT_INFO("S1 after rotation2 %3.10f %3.10f %3.10f\n", tmpS1Norm[0], tmpS1Norm[1], tmpS1Norm[2]);
    // XLAL_PRINT_INFO("S2 after rotation2 %3.10f %3.10f %3.10f\n", tmpS2Norm[0], tmpS2Norm[1], tmpS2Norm[2]);
	/*
	 * STEP 4) In the L0-N0 frame, we calculate (dE/dr)|sph using Eq.
	 * (4.14), then initial dr/dt using Eq. (4.10), and finally pr0 using
	 * Eq. (4.15).
	 */

	/* Now we can calculate the flux. Change to spherical co-ords */
	CartesianToSpherical(qSph, pSph, qCart, pCart);
	memcpy(sphValues, qSph, sizeof(qSph));
	memcpy(sphValues + 3, pSph, sizeof(pSph));
	memcpy(sphValues + 6, tmpS1, sizeof(tmpS1));
	memcpy(sphValues + 9, tmpS2, sizeof(tmpS2));

	memcpy(cartValues, qCart, sizeof(qCart));
	memcpy(cartValues + 3, pCart, sizeof(pCart));
	memcpy(cartValues + 6, tmpS1, sizeof(tmpS1));
	memcpy(cartValues + 9, tmpS2, sizeof(tmpS2));

	REAL8		dHdpphi , d2Hdr2, d2Hdrdpphi;
	REAL8		rDot    , dHdpr, flux, dEdr;

	d2Hdr2 = XLALCalculateSphHamiltonianDeriv2Prec(0, 0, sphValues, params);
	d2Hdrdpphi = XLALCalculateSphHamiltonianDeriv2Prec(0, 5, sphValues, params);

    // XLAL_PRINT_INFO("d2Hdr2 = %.16e, d2Hdrdpphi = %.16e\n", d2Hdr2, d2Hdrdpphi);

	/* New code to compute derivatives w.r.t. cartesian variables */

	REAL8		tmpDValues[14];
	int status;
	for (i = 0; i < 3; i++) {
		cartValues[i + 6] /= mTotal * mTotal;
		cartValues[i + 9] /= mTotal * mTotal;
	}
    UINT oldignoreflux = params->ignoreflux;
    params->ignoreflux = 1;
    status = XLALSpinPrecHcapNumericalDerivative(0, cartValues, tmpDValues, params);
    params->ignoreflux = oldignoreflux;
	for (i = 0; i < 3; i++) {
		cartValues[i + 6] *= mTotal * mTotal;
		cartValues[i + 9] *= mTotal * mTotal;
	}

	dHdpphi = tmpDValues[1] / sqrt(cartValues[0] * cartValues[0] + cartValues[1] * cartValues[1] + cartValues[2] * cartValues[2]);
	//XLALSpinPrecHcapNumDerivWRTParam(4, cartValues, params) / sphValues[0];

	dEdr = -dHdpphi * d2Hdr2 / d2Hdrdpphi;

    // XLAL_PRINT_INFO("d2Hdr2 = %.16e d2Hdrdpphi = %.16e dHdpphi = %.16e\n",
    //     d2Hdr2, d2Hdrdpphi, dHdpphi);

	if (d2Hdr2 != 0.0 && e0 == 0.0) 
    {
		/* We will need to calculate the Hamiltonian to get the flux */
		s1Vec.length = s2Vec.length = s1VecNorm.length = s2VecNorm.length = sKerr.length = sStar.length = 3;
		s1Vec.data = tmpS1;
		s2Vec.data = tmpS2;
		s1VecNorm.data = tmpS1Norm;
		s2VecNorm.data = tmpS2Norm;
		sKerr.data = sKerrData;
		sStar.data = sStarData;

		qCartVec.length = pCartVec.length = 3;
		qCartVec.data = qCart;
		pCartVec.data = pCart;

		//chi1 = tmpS1[0] * LnHat[0] + tmpS1[1] * LnHat[1] + tmpS1[2] * LnHat[2];
		//chi2 = tmpS2[0] * LnHat[0] + tmpS2[1] * LnHat[1] + tmpS2[2] * LnHat[2];

		//if (debugPK)
			//XLAL_PRINT_INFO("magS1 = %.16e, magS2 = %.16e\n", chi1, chi2);

		//chiS = 0.5 * (chi1 / (mass1 * mass1) + chi2 / (mass2 * mass2));
		//chiA = 0.5 * (chi1 / (mass1 * mass1) - chi2 / (mass2 * mass2));

		EOBCalculateSigmaKerr(&sKerr, &s1VecNorm, &s2VecNorm);
		EOBCalculateSigmaStar(&sStar, mass1, mass2, &s1VecNorm, &s2VecNorm);

		/*
		 * The a in the flux has been set to zero, but not in the
		 * Hamiltonian
		 */
		a = sqrt(sKerr.data[0] * sKerr.data[0] + sKerr.data[1] * sKerr.data[1] + sKerr.data[2] * sKerr.data[2]);
		//XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients(params->eobParams->hCoeffs, mass1, mass2, eta, /* a */ 0.0, chiS, chiA);
		//XLALSimIMRCalculateSpinPrecEOBHCoeffs(params->seobCoeffs, eta, a);
		ham = XLALSimIMRSpinPrecEOBHamiltonian(eta, &qCartVec, &pCartVec, &s1VecNorm, &s2VecNorm, &sKerr, &sStar, params->tortoise, params->seobCoeffs, params->hParams);

        // XLAL_PRINT_INFO("Stas: hamiltonian in ICs at this point is %.16e\n", ham);

		/* And now, finally, the flux */
		REAL8Vector	polarDynamics, cartDynamics;
		REAL8		polarData[4], cartData[12];

		polarDynamics.length = 4;
		polarDynamics.data = polarData;

		polarData[0] = qSph[0];
		polarData[1] = 0.;
		polarData[2] = pSph[0];
		polarData[3] = pSph[2];

		cartDynamics.length = 12;
		cartDynamics.data = cartData;

		memcpy(cartData, qCart, 3 * sizeof(REAL8));
		memcpy(cartData + 3, pCart, 3 * sizeof(REAL8));
		memcpy(cartData + 6, tmpS1Norm, 3 * sizeof(REAL8));
		memcpy(cartData + 9, tmpS2Norm, 3 * sizeof(REAL8));

		//XLAL_PRINT_INFO("Stas: starting FLux calculations\n");
        flux = XLALInspiralPrecSpinFactorizedFlux(&polarDynamics, &cartDynamics, nqcCoeffs, omega, 0, qSph[0]*omega, params, ham, lMax, SpinAlignedEOBversion);
		/*
		 * flux  = XLALInspiralSpinFactorizedFlux( &polarDynamics,
		 * nqcCoeffs, omega, params, ham, lMax, SpinAlignedEOBversion
		 * );
		 */
		PRINT_LOG_INFO(LOG_DEBUG, "Initial Conditions: Stas flux = %.16e", flux);
		//exit(0);
		// flux = flux / eta;
        flux = 0.0;
        rDot = 0.0;
        // if (e0 != 0.0)
        //     rDot = 0.0;
        // else
		//     rDot = -flux / dEdr;
		/*
		 * We now need dHdpr - we take it that it is safely linear up
		 * to a pr of 1.0e-3 PK: Ideally, the pr should be of the
		 * order of other momenta, in order for its contribution to
		 * the Hamiltonian to not get buried in the numerical noise
		 * in the numerically larger momenta components
		 */
		cartValues[3] = 1.0e-3;
		for (i = 0; i < 3; i++) 
        {
			cartValues[i + 6] /= mTotal * mTotal;
			cartValues[i + 9] /= mTotal * mTotal;
		}
        oldignoreflux = params->ignoreflux;
        params->ignoreflux = 1;
        params->seobCoeffs->updateHCoeffs = 1;
	    status = XLALSpinPrecHcapNumericalDerivative(0, cartValues, tmpDValues, params);
        params->ignoreflux = oldignoreflux;
		for (i = 0; i < 3; i++) 
        {
			cartValues[i + 6] *= mTotal * mTotal;
			cartValues[i + 9] *= mTotal * mTotal;
		}
        REAL8   csi = sqrt(XLALSimIMRSpinPrecEOBHamiltonianDeltaT(params->seobCoeffs, qSph[0], eta, a)*XLALSimIMRSpinPrecEOBHamiltonianDeltaR(params->seobCoeffs, qSph[0], eta, a)) / (qSph[0] * qSph[0] + a * a);

		dHdpr = csi*tmpDValues[0];
		//XLALSpinPrecHcapNumDerivWRTParam(3, cartValues, params);

        // XLAL_PRINT_INFO("Ingredients going into prDot:\n");
        // XLAL_PRINT_INFO("flux = %.16e, dEdr = %.16e, dHdpr = %.16e, dHdpr/pr = %.16e\n", flux, dEdr, dHdpr, dHdpr / cartValues[3]);
		/*
		 * We can now calculate what pr should be taking into account
		 * the flux
		 */
		pSph[0] = rDot / (dHdpr / cartValues[3]);
	} else 
    {
		/*
		 * Since d2Hdr2 has evaluated to zero, we cannot do the
		 * above. Just set pr to zero
		 */
		//XLAL_PRINT_INFO("d2Hdr2 is zero!\n");
		pSph[0] = 0;
	}

	/* Now we are done - convert back to cartesian coordinates ) */
	SphericalToCartesian(qCart, pCart, qSph, pSph);
    PRINT_LOG_INFO(LOG_DEBUG, "Sph initial condition : r = (%e,%e,%e), p = (%e,%e,%e)", qSph[0], qSph[1], qSph[2], pSph[0], pSph[1], pSph[2]);

	/*
	 * STEP 5) Rotate back to the original inertial frame by inverting
	 * the rotation of STEP 3 and then  inverting the rotation of STEP 1.
	 */

	/* Undo rotations to get back to the original basis */
	/* Second rotation */
	ApplyRotationMatrix(invMatrix2, rHat);
	ApplyRotationMatrix(invMatrix2, vHat);
	ApplyRotationMatrix(invMatrix2, LnHat);
	ApplyRotationMatrix(invMatrix2, tmpS1);
	ApplyRotationMatrix(invMatrix2, tmpS2);
	ApplyRotationMatrix(invMatrix2, tmpS1Norm);
	ApplyRotationMatrix(invMatrix2, tmpS2Norm);
	ApplyRotationMatrix(invMatrix2, qCart);
	ApplyRotationMatrix(invMatrix2, pCart);

	/* First rotation */
	ApplyRotationMatrix(invMatrix, rHat);
	ApplyRotationMatrix(invMatrix, vHat);
	ApplyRotationMatrix(invMatrix, LnHat);
	ApplyRotationMatrix(invMatrix, tmpS1);
	ApplyRotationMatrix(invMatrix, tmpS2);
	ApplyRotationMatrix(invMatrix, tmpS1Norm);
	ApplyRotationMatrix(invMatrix, tmpS2Norm);
	ApplyRotationMatrix(invMatrix, qCart);
	ApplyRotationMatrix(invMatrix, pCart);

    gsl_matrix_free(invMatrix);
    gsl_matrix_free(invMatrix2);

        /* If required, apply the tortoise transform */
	if (tmpTortoise) 
    {
		REAL8		r = sqrt(qCart[0] * qCart[0] + qCart[1] * qCart[1] + qCart[2] * qCart[2]);
		REAL8		deltaR = XLALSimIMRSpinPrecEOBHamiltonianDeltaR(params->seobCoeffs, r, eta, a);
		REAL8		deltaT = XLALSimIMRSpinPrecEOBHamiltonianDeltaT(params->seobCoeffs, r, eta, a);
		REAL8		csi = sqrt(deltaT * deltaR) / (r * r + a * a);

		REAL8		pr = (qCart[0] * pCart[0] + qCart[1] * pCart[1] + qCart[2] * pCart[2]) / r;

		params->tortoise = tmpTortoise;

        // XLAL_PRINT_INFO("Applying the tortoise to p (csi = %.26e)\n", csi);
        // XLAL_PRINT_INFO("pCart = %3.10f %3.10f %3.10f\n", pCart[0], pCart[1], pCart[2]);
		for (i = 0; i < 3; i++) 
        {
			pCart[i] = pCart[i] + qCart[i] * pr * (csi - 1.) / r;
		}
	}
    qCart[0] /= 1.-e0;
    pCart[1] *= 1.-e0;
    
    /* Now copy the initial conditions back to the return vector */
	memcpy(initConds->data, qCart, sizeof(qCart));
	memcpy(initConds->data + 3, pCart, sizeof(pCart));
	memcpy(initConds->data + 6, tmpS1Norm, sizeof(tmpS1Norm));
	memcpy(initConds->data + 9, tmpS2Norm, sizeof(tmpS2Norm));

    for (i=0; i<12; i++) 
    {
        if (fabs(initConds->data[i]) <=1.0e-15) 
        {
            initConds->data[i] = 0.;
        }
    }

    // XLAL_PRINT_INFO("THE FINAL INITIAL CONDITIONS:\n");
    // XLAL_PRINT_INFO(" %.16e %.16e %.16e\n%.16e %.16e %.16e\n%.16e %.16e %.16e\n%.16e %.16e %.16e\n", initConds->data[0], initConds->data[1], initConds->data[2],
    //         initConds->data[3], initConds->data[4], initConds->data[5], initConds->data[6], initConds->data[7], initConds->data[8],
    //         initConds->data[9], initConds->data[10], initConds->data[11]);
    return CEV_SUCCESS;
}


/*--------------------------------------------------------------*/
/*                                                              */
/*                                                              */
/*                                                              */
/*                           egw init                           */
/*                                                              */
/*                                                              */
/*                                                              */
/*--------------------------------------------------------------*/

/* Spin EOB root */
typedef
struct tagSEOBRootParams_egw
{
    // REAL8          values[12]; /**<< Dynamical variables, x, y, z, px, py, pz, S1x, S1y, S1z, S2x, S2y and S2z */
    SpinEOBParams *params;     /**<< Spin EOB parameters -- physical, pre-computed, etc. */
    REAL8          omega22;      /**<< mode 22 frequency */
    REAL8          e22;
}
SEOBRootParams_egw;

REAL8 calc_e22_from_egw(REAL8 egw);
INT find_SApphi_from_rpm(REAL8 rp, REAL8 rm, SpinEOBParams *core, REAL8 *pphi);
INT EvaluateOmega22SA_form_rpphi(REAL8 r, REAL8 pphi, SpinEOBParams *seobParams, REAL8 *omega22);

static REAL8 Xscale1 = 100, Xscale2 = 100;

static REAL8 calculate_e22_from_omegas(REAL8 omg22p, REAL8 omg22m)
{
    REAL8 sqm, sqp;
    sqm = sqrt(omg22m);
    sqp = sqrt(omg22p);
    return fabs(sqm - sqp) / (sqm + sqp);
}

/**
 * @brief Input (r, pphi), return (omega22, egw)
 */
static int
XLALFindSphericalOrbitByegw(const gsl_vector *x, /**<< Parameters requested by gsl root finder */
                       void *params,        /**<< Spin EOB parameters */
                       gsl_vector *f        /**<< Function values for the given parameters */
)
{
    SEOBRootParams_egw *rootParams = (SEOBRootParams_egw *) params;
    
    REAL8 pphi;
    REAL8 rm, rp;
    /* Numerical derivative of Hamiltonian wrt given value */
    REAL8 omg22, e22, omg22p, omg22m;
    /* Populate the appropriate values */
    /* In the special theta=pi/2 phi=0 case, r is x */
    rm  = gsl_vector_get( x, 0 )/Xscale1;
    rp = gsl_vector_get( x, 1 )/Xscale2;
    if (rm<5. || rp<5.)
        return GSL_FAILURE;
    find_SApphi_from_rpm(rp, rm, rootParams->params, &pphi);
    print_debug("rm = %.16e, rp = %.16e, pphi = %.16e\n", rm, rp, pphi);
    EvaluateOmega22SA_form_rpphi(rm, pphi, rootParams->params, &omg22m);
    EvaluateOmega22SA_form_rpphi(rm, pphi, rootParams->params, &omg22p);
    e22 = calculate_e22_from_omegas(omg22p, omg22m);
    // /* populate the function vector */
    gsl_vector_set( f, 0, (omg22 - rootParams->omega22));
    gsl_vector_set( f, 1, (e22 - rootParams->e22));
    
    return CEV_SUCCESS;
}


// v2
typedef struct {
    REAL8 omega;
    REAL8 egw;
    SpinEOBParams *params;
}EgwRootParamsV2;

static int
XLALFindSphericalOrbitSA(const gsl_vector *x, /**<< Parameters requested by gsl root finder */
                       void *params,        /**<< Spin EOB parameters */
                       gsl_vector *f        /**<< Function values for the given parameters */
)
{
    SEOBRootParams *rootParams = (SEOBRootParams *) params;
    
    REAL8 py, r, pphi;
    
    /* Numerical derivative of Hamiltonian wrt given value */
    REAL8 dHdx, dHdpy;
    REAL8 dHdr, dHdpphi;
    REAL8 e0;
    /* Populate the appropriate values */
    /* In the special theta=pi/2 phi=0 case, r is x */
    rootParams->values[0] = r  = gsl_vector_get( x, 0 );
    rootParams->values[4] = py = gsl_vector_get( x, 1 );
    e0 = rootParams->e0;
    // print_debug( "Values r = %.16e, py = %.16e\n", r, py );
    
    pphi   = r * py;
    
    /* dHdR */
    dHdx = XLALSpinHcapNumDerivWRTParam( 0, rootParams->values, rootParams->params );
    if ( IS_REAL8_FAIL_NAN( dHdx ) )
    {
        return CEV_FAILURE;
    }
    //printf( "dHdx = %.16e\n", dHdx );
    
    /* dHdPphi (I think we can use dHdPy in this coord system) */
    /* TODO: Check this is okay */
    dHdpy = XLALSpinHcapNumDerivWRTParam( 4, rootParams->values, rootParams->params );
    if ( IS_REAL8_FAIL_NAN( dHdpy ) )
    {
        return CEV_FAILURE;
    }
        
    /* Now convert to spherical polars */
    dHdr      = dHdx - dHdpy * pphi / (r*r);
    dHdpphi   = dHdpy / r;
    
    /* populate the function vector */
    gsl_vector_set( f, 0, dHdr + e0/r/r);
    gsl_vector_set( f, 1, dHdpphi - rootParams->omega );
    
    //printf( "Current funcvals = %.16e %.16e %.16e\n", gsl_vector_get( f, 0 ), gsl_vector_get( f, 1 ),
    //  gsl_vector_get( f, 2 )/*dHdpphi*/ );
    
    return CEV_SUCCESS;
}


int calc_rpphi_from_e(REAL8 e, REAL8 omega, SpinEOBParams *params, REAL8 *r, REAL8 *pphi)
{
    SEOBRootParams rootParams;
    const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver *rootSolver = NULL;

    gsl_multiroot_function rootFunction;
    gsl_vector *initValues  = NULL;
    gsl_vector *finalValues = NULL;
    int gslStatus;
    const int maxIter = 100;

    REAL8 v0    = GET_CBRT( omega );
    rootParams.e0 = e;
    rootParams.omega  = omega;
    rootParams.params = params;
    memset(rootParams.values, 0, sizeof(rootParams.values));
    rootParams.values[0] = (1.-e) / (v0*v0);
    rootParams.values[4] = pow(1.-e, 0.5) * v0;
    for (int i=0; i<3; i++)
    {
        rootParams.values[i+6] = params->s1Vec->data[i]*params->m1*params->m1;
        rootParams.values[i+9] = params->s2Vec->data[i]*params->m2*params->m2;
    }
    rootSolver = gsl_multiroot_fsolver_alloc( T, 2 );
    if ( !rootSolver )
    {
        return CEV_FAILURE;
    }
    
    initValues = gsl_vector_calloc( 2 );
    if ( !initValues )
    {
        gsl_multiroot_fsolver_free( rootSolver );
        return CEV_FAILURE;
    }
    gsl_vector_set( initValues, 0, rootParams.values[0] );
    gsl_vector_set( initValues, 1, rootParams.values[4] ); // py
    rootFunction.f      = XLALFindSphericalOrbitSA;
    rootFunction.n      = 2;
    rootFunction.params = &rootParams;
    gsl_multiroot_fsolver_set( rootSolver, &rootFunction, initValues );

    int i = 0;
    do
    {
        gslStatus = gsl_multiroot_fsolver_iterate( rootSolver );
        if ( gslStatus != GSL_SUCCESS )
        {
            print_warning( "Error in GSL iteration function!\n" );
            gsl_multiroot_fsolver_free( rootSolver );
            gsl_vector_free( initValues );
            return CEV_FAILURE;
        }
        gslStatus = gsl_multiroot_test_residual( rootSolver->f, 1.0e-10 );
        i++;
    }
    while ( gslStatus == GSL_CONTINUE && i <= maxIter );

    if ( i > maxIter && gslStatus != GSL_SUCCESS )
    {
        gsl_multiroot_fsolver_free( rootSolver );
        gsl_vector_free( initValues );
        return CEV_FAILURE;
    }
    finalValues = gsl_multiroot_fsolver_root( rootSolver );
    *r = gsl_vector_get( finalValues, 0 );
    *pphi = gsl_vector_get( finalValues, 1 ) * gsl_vector_get( finalValues, 0 );
    gsl_multiroot_fsolver_free( rootSolver );
    gsl_vector_free( initValues );
    return CEV_SUCCESS;
}

INT CalculateSAh22SeriesFromrpphi(REAL8 r, REAL8 pphi, SpinEOBParams *core,
    REAL8 *omega22, REAL8 *egw);

static double GSL_findEccentricity(const double x, void *params)
{
    EgwRootParamsV2 *rootpms = (EgwRootParamsV2*) params;
    int status;
    REAL8 r, pphi;
    status = calc_rpphi_from_e(x, rootpms->omega, rootpms->params, &r, &pphi);
    REAL8 omega22, egw;
    status = CalculateSAh22SeriesFromrpphi(r, pphi, rootpms->params, &omega22, &egw);
    if (status != CEV_SUCCESS) return 0;
    // print_debug("x = %.16e, egw = %.16e\n", x, egw);
    return rootpms->egw - egw;
}

static INT Grid_Solve_Initial_egw(EgwRootParamsV2 *rootpms, REAL8 *ret, REAL8 x_lo, REAL8 x_hi)
{
    size_t n_sample = 1000;
    // REAL8Vector *eccVec = CreateREAL8Vector(n_sample);
    REAL8 this_e, find_e, this_eps, find_eps;
    // REAL8Vector *epsVec = CreateREAL8Vector(n_sample);
    REAL8 dx = (x_hi - x_lo) / (n_sample-1);
    find_e = x_lo;
    find_eps = fabs(GSL_findEccentricity(find_e, (void*)rootpms));
    for (int i = 1; i < n_sample; i++)
    {
        this_e = x_hi + dx * i;
        this_eps = fabs(GSL_findEccentricity(this_e, (void*)rootpms));
        if (this_eps < find_eps)
        {
            find_e = this_e;
            find_eps = this_eps;
        }
    }
    *ret = find_e;
    // print_debug("find_e = %.16e\n", find_e);
    if (find_eps < 1e-3)
        return CEV_SUCCESS;
    return CEV_FAILURE;
}

/* Initial Condition Core function */
/* For nonprecessing case */
INT EOBInitialConditionsSA_egw(REAL8Vector    *initConds,
                         const REAL8    mass1,
                         const REAL8    mass2,
                         const REAL8    fMin,
                         const REAL8    ecc,
                         const REAL8    inc,
                         const REAL8    spin1[],
                         const REAL8    spin2[],
                         SpinEOBParams  *params)
{
    if (!initConds)
        return CEV_FAILURE;
    static const int lMax = 8;
    INT i;
    int tmpTortoise;
    /* non-zero eccentricity - Start frequency correction */
    // REAL8 ecc = params->eccentricity;
        
    REAL8 mTotal;
    REAL8 eta;
    REAL8 omega, v0;   /* Initial velocity and angular frequency */
    
    REAL8 ham;      /* Hamiltonian */
    
    REAL8 LnHat[3]; /* Initial orientation of angular momentum */
    REAL8 rHat[3];  /* Initial orientation of radial vector */
    REAL8 vHat[3];  /* Initial orientation of velocity vector */
    REAL8 Lhat[3];  /* Direction of relativistic ang mom */
    REAL8 qHat[3];
    REAL8 pHat[3];
    
    /* q and p vectors in Cartesian and spherical coords */
    REAL8 qCart[3], pCart[3];
    REAL8 qSph[3], pSph[3];
    
    /* We will need to manipulate the spin vectors */
    /* We will use temporary vectors to do this */
    REAL8 tmpS1[3];
    REAL8 tmpS2[3];
    REAL8 tmpS1Norm[3];
    REAL8 tmpS2Norm[3];
    
    REAL8Vector qCartVec, pCartVec;
    REAL8Vector s1Vec, s2Vec, s1VecNorm, s2VecNorm;
    REAL8Vector sKerr, sStar;
    REAL8       sKerrData[3], sStarData[3];
    REAL8       a = 0.; //, chiS, chiA;
    //REAL8       chi1, chi2;
    
    /* We will need a full values vector for calculating derivs of Hamiltonian */
    REAL8 sphValues[12];
    REAL8 cartValues[12];
    
    /* Matrices for rotating to the new basis set. */
    /* It is more convenient to calculate the ICs in a simpler basis */
    gsl_matrix *rotMatrix  = NULL;
    gsl_matrix *invMatrix  = NULL;
    gsl_matrix *rotMatrix2 = NULL;
    gsl_matrix *invMatrix2 = NULL;
    
    /* Root finding stuff for finding the spherical orbit */
    EgwRootParamsV2 rootParams;
    const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
    gsl_root_fsolver *rootSolver = gsl_root_fsolver_alloc (T);
    gsl_function F;

    int gslStatus;
    const int maxIter = 100;
    
    memset( &rootParams, 0, sizeof( rootParams ) );
    
    mTotal = mass1 + mass2;
    eta    = mass1 * mass2 / (mTotal * mTotal);
    memcpy( tmpS1, spin1, sizeof(tmpS1) );
    memcpy( tmpS2, spin2, sizeof(tmpS2) );
    memcpy( tmpS1Norm, spin1, sizeof(tmpS1Norm) );
    memcpy( tmpS2Norm, spin2, sizeof(tmpS2Norm) );
    for ( i = 0; i < 3; i++ )
    {
        tmpS1Norm[i] /= mTotal * mTotal;
        tmpS2Norm[i] /= mTotal * mTotal;
    }
    // eobVersion = params->seobCoeffs->eobVersion;
    /* We compute the ICs for the non-tortoise p, and convert at the end */
    tmpTortoise      = params->tortoise;
    params->tortoise = 0;
    
    EOBNonQCCoeffs *nqcCoeffs = NULL;
    nqcCoeffs = params->nqcCoeffs;
    
    /* STEP 1) Rotate to LNhat0 along z-axis and N0 along x-axis frame, where LNhat0 and N0 are initial normal to
     *         orbital plane and initial orbital separation;
     */
    
    /* Set the initial orbital ang mom direction. Taken from STPN code */
    LnHat[0] = GET_SIN(inc);
    LnHat[1] = 0.;
    LnHat[2] = GET_COS(inc);
    
    /* Set the radial direction - need to take care to avoid singularity if L is along z axis */
    if ( LnHat[2] > 0.9999 )
    {
        rHat[0] = 1.;
        rHat[1] = rHat[2] = 0.;
    }
    else
    {
        REAL8 theta0 = GET_ATAN( - LnHat[2] / LnHat[0] ); /* theta0 is between 0 and Pi */
        rHat[0] = GET_SIN( theta0 );
        rHat[1] = 0;
        rHat[2] = GET_COS( theta0 );
    }
    
    /* Now we can complete the triad */
    vHat[0] = CalculateCrossProduct( 0, LnHat, rHat );
    vHat[1] = CalculateCrossProduct( 1, LnHat, rHat );
    vHat[2] = CalculateCrossProduct( 2, LnHat, rHat );
    
    NormalizeVector( vHat );
    
    /* XXX Test code XXX */
    /*for ( i = 0; i < 3; i++ )
     {
     printf ( " LnHat[%d] = %.16e, rHat[%d] = %.16e, vHat[%d] = %.16e\n", i, LnHat[i], i, rHat[i], i, vHat[i] );
     }
     
     printf("\n\n" );
     for ( i = 0; i < 3; i++ )
     {
     printf ( " s1[%d] = %.16e, s2[%d] = %.16e\n", i, tmpS1[i], i, tmpS2[i] );
     }*/
    
    /* Allocate and compute the rotation matrices */
    rotMatrix = gsl_matrix_alloc( 3, 3 );
    invMatrix = gsl_matrix_alloc( 3, 3 );
    if ( !rotMatrix || !invMatrix )
    {
        if ( rotMatrix ) gsl_matrix_free( rotMatrix );
        if ( invMatrix ) gsl_matrix_free( invMatrix );
        return CEV_FAILURE;
    }
    
    if ( CalculateRotationMatrix( rotMatrix, invMatrix, rHat, vHat, LnHat ) == CEV_FAILURE )
    {
        gsl_matrix_free( rotMatrix );
        gsl_matrix_free( invMatrix );
        return CEV_FAILURE;
    }
    
    /* Rotate the orbital vectors and spins */
    ApplyRotationMatrix( rotMatrix, rHat );
    ApplyRotationMatrix( rotMatrix, vHat );
    ApplyRotationMatrix( rotMatrix, LnHat );
    ApplyRotationMatrix( rotMatrix, tmpS1 );
    ApplyRotationMatrix( rotMatrix, tmpS2 );
    ApplyRotationMatrix( rotMatrix, tmpS1Norm );
    ApplyRotationMatrix( rotMatrix, tmpS2Norm );
    
    /* XXX Test code XXX */
    /*printf( "\nAfter applying rotation matrix:\n\n" );
     for ( i = 0; i < 3; i++ )
     {
     printf ( " LnHat[%d] = %.16e, rHat[%d] = %.16e, vHat[%d] = %.16e\n", i, LnHat[i], i, rHat[i], i, vHat[i] );
     }
     
     printf("\n\n" );
     for ( i = 0; i < 3; i++ )
     {
     printf ( " s1[%d] = %.16e, s2[%d] = %.16e\n", i, tmpS1[i], i, tmpS2[i] );
     }*/
    
    /* STEP 2) After rotation in STEP 1, in spherical coordinates, phi0 and theta0 are given directly in Eq. (4.7),
     *         r0, pr0, ptheta0 and pphi0 are obtained by solving Eqs. (4.8) and (4.9) (using gsl_multiroot_fsolver).
     *         At this step, we find initial conditions for a spherical orbit without radiation reaction.
     */
    // X.L.: here we solve equation (e22, omegaMinus)(rp, rm)
    /* Calculate the initial velocity from the given initial frequency */
    // initial orbital angular velocity
    omega = CST_PI * mTotal * CST_MTSUN_SI * fMin;
#if 0
    /* Given this, we can start to calculate the initial conditions */
    /* for spherical coords in the new basis */
    rootParams.omega  = omega;
    rootParams.params = params;
    rootParams.egw = ecc;
    F.function = &GSL_findEccentricity;
    F.params = &rootParams;
    REAL8 x_lo, x_hi, root;
    x_lo = GET_MAX(0.001, ecc - 0.2);
    x_hi = GET_MIN(0.9, ecc + 0.2);
    /* Initialise the gsl stuff */
    GSL_START;
    int status, check_iter, max_check_iter = 10;
    for(check_iter=0; check_iter<max_check_iter; check_iter++)
    {
        status = gsl_root_fsolver_set (rootSolver, &F, x_lo, x_hi);
        if (status != GSL_SUCCESS)
        {
            x_hi = GET_MIN(0.9, x_hi + 0.1);
        } else 
            break;
    }
    if (check_iter==max_check_iter)
    {
        print_warning( "Error in GSL solver setting!\n" );
        gsl_root_fsolver_free( rootSolver );
        gsl_matrix_free( rotMatrix );
        gsl_matrix_free( invMatrix );
        GSL_END;
        return CEV_FAILURE;
    }
    GSL_END;

    /* We are now ready to iterate to find the solution */
    i = 0;
    do
    {
        gslStatus = gsl_root_fsolver_iterate( rootSolver );
        if ( gslStatus != GSL_SUCCESS )
        {
            // print_warning( "Error in GSL iteration function!\n" );
            // gsl_root_fsolver_free( rootSolver );
            // gsl_matrix_free( rotMatrix );
            // gsl_matrix_free( invMatrix );
            // return CEV_FAILURE;
            status = Grid_Solve_Initial_egw(&rootParams, &root, x_lo, x_hi);
            if (status != CEV_SUCCESS)
            {
                gsl_root_fsolver_free( rootSolver );
                gsl_matrix_free( rotMatrix );
                gsl_matrix_free( invMatrix );
                return CEV_FAILURE;
            }
            break;
        }
        x_lo = gsl_root_fsolver_x_lower (rootSolver);
        x_hi = gsl_root_fsolver_x_upper (rootSolver);
        gslStatus = gsl_root_test_interval (x_lo, x_hi, 1e-9, 1e-9);
        i++;
    }
    while ( gslStatus == GSL_CONTINUE && i <= maxIter );

    if ( i > maxIter && gslStatus != GSL_SUCCESS )
    {
        // gsl_root_fsolver_free( rootSolver );
        // gsl_matrix_free( rotMatrix );
        // gsl_matrix_free( invMatrix );
        // return CEV_FAILURE;
        status = Grid_Solve_Initial_egw(&rootParams, &root, x_lo, x_hi);
        if (status != CEV_SUCCESS)
        {
            gsl_root_fsolver_free( rootSolver );
            gsl_matrix_free( rotMatrix );
            gsl_matrix_free( invMatrix );
            return CEV_FAILURE;
        }
    } else 
        root = gsl_root_fsolver_root (rootSolver);


    // print_log( "Spherical orbit conditions here given by the following:\n" );
    //  print_err( " x = %.16e, py = %.16e, pz = %.16e\n", gsl_vector_get( finalValues, 0 ),
    //  gsl_vector_get( finalValues, 1 ), gsl_vector_get( finalValues, 2 ) );
    // print_debug("root = %.16e\n", root);
#endif
    memset( qCart, 0, sizeof(qCart) );
    memset( pCart, 0, sizeof(pCart) );
    // calc_rpphi_from_e(root, omega, params, &(qCart[0]), &(pCart[1]));
    calc_rpphi_from_e(ecc, omega, params, &(qCart[0]), &(pCart[1]));
    // qCart[0] = gsl_vector_get( finalValues, 0 );
    // pCart[1] = gsl_vector_get( finalValues, 1 );
    pCart[1] /= qCart[0];
    /* Free the GSL root finder, since we're done with it */
    gsl_root_fsolver_free( rootSolver );
    
    /* STEP 3) Rotate to L0 along z-axis and N0 along x-axis frame, where L0 is the initial orbital angular momentum
     *         and L0 is calculated using initial position and linear momentum obtained in STEP 2.
     */
    
    /* Now we can calculate the relativistic L and rotate to a new basis */
    memcpy( qHat, qCart, sizeof(qCart) );
    memcpy( pHat, pCart, sizeof(pCart) );
    
    NormalizeVector( qHat );
    NormalizeVector( pHat );
    
    Lhat[0] = CalculateCrossProduct( 0, qHat, pHat );
    Lhat[1] = CalculateCrossProduct( 1, qHat, pHat );
    Lhat[2] = CalculateCrossProduct( 2, qHat, pHat );
    
    NormalizeVector( Lhat );
    
    rotMatrix2 = gsl_matrix_alloc( 3, 3 );
    invMatrix2 = gsl_matrix_alloc( 3, 3 );
    
    if ( CalculateRotationMatrix( rotMatrix2, invMatrix2, qHat, pHat, Lhat ) == CEV_FAILURE )
    {
        gsl_matrix_free( rotMatrix );
        gsl_matrix_free( invMatrix );
        return CEV_FAILURE;
    }
    
    ApplyRotationMatrix( rotMatrix2, rHat );
    ApplyRotationMatrix( rotMatrix2, vHat );
    ApplyRotationMatrix( rotMatrix2, LnHat );
    ApplyRotationMatrix( rotMatrix2, tmpS1 );
    ApplyRotationMatrix( rotMatrix2, tmpS2 );
    ApplyRotationMatrix( rotMatrix2, tmpS1Norm );
    ApplyRotationMatrix( rotMatrix2, tmpS2Norm );
    ApplyRotationMatrix( rotMatrix2, qCart );
    ApplyRotationMatrix( rotMatrix2, pCart );
    
    /* STEP 4) In the L0-N0 frame, we calculate (dE/dr)|sph using Eq. (4.14), then initial dr/dt using Eq. (4.10),
     *         and finally pr0 using Eq. (4.15).
     */
    
    /* Now we can calculate the flux. Change to spherical co-ords */
    CartesianToSpherical( qSph, pSph, qCart, pCart );
    memcpy( sphValues, qSph, sizeof( qSph ) );
    memcpy( sphValues+3, pSph, sizeof( pSph ) );
    memcpy( sphValues+6, tmpS1, sizeof(tmpS1) );
    memcpy( sphValues+9, tmpS2, sizeof(tmpS2) );
    
    memcpy( cartValues, qCart, sizeof(qCart) );
    memcpy( cartValues+3, pCart, sizeof(pCart) );
    memcpy( cartValues+6, tmpS1, sizeof(tmpS1) );
    memcpy( cartValues+9, tmpS2, sizeof(tmpS2) );
    
    REAL8 dHdpphi, d2Hdr2, d2Hdrdpphi;
    REAL8 rDot, dHdpr, flux, dEdr;
    
    d2Hdr2 = XLALCalculateSphHamiltonianDeriv2( 0, 0, sphValues, params );
    d2Hdrdpphi = XLALCalculateSphHamiltonianDeriv2( 0, 5, sphValues, params );
    dHdpphi = XLALSpinHcapNumDerivWRTParam( 4, cartValues, params ) / sphValues[0];
    dEdr  = - dHdpphi * d2Hdr2 / d2Hdrdpphi;
        
    if ( d2Hdr2 != 0.0 && ecc==0.0 )
    {
        /* We will need to calculate the Hamiltonian to get the flux */
        s1Vec.length = s2Vec.length = s1VecNorm.length = s2VecNorm.length = sKerr.length = sStar.length = 3;
        s1Vec.data = tmpS1;
        s2Vec.data = tmpS2;
        s1VecNorm.data = tmpS1Norm;
        s2VecNorm.data = tmpS2Norm;
        sKerr.data = sKerrData;
        sStar.data = sStarData;
        
        qCartVec.length = pCartVec.length = 3;
        qCartVec.data   = qCart;
        pCartVec.data   = pCart;
        
        //chi1 = tmpS1[0]*LnHat[0] + tmpS1[1]*LnHat[1] + tmpS1[2]*LnHat[2];
        //chi2 = tmpS2[0]*LnHat[0] + tmpS2[1]*LnHat[1] + tmpS2[2]*LnHat[2];
        
        // print_debug( "m1 = %g, m2 = %g\n", mass1, mass2 );
        // print_debug( "eta = %g\n", eta);
        // print_debug( "qCartVec = (%g, %g, %g)\n", qCart[0], qCart[1], qCart[2]);
        // print_debug( "pCartVec = (%g, %g, %g)\n", pCart[0], pCart[1], pCart[2]);
        // print_debug( "s1Vec = (%g, %g, %g)\n", tmpS1[0], tmpS1[1], tmpS1[2]);
        // print_debug( "s2Vec = (%g, %g, %g)\n", tmpS2[0], tmpS2[1], tmpS2[2]);
        
        //chiS = 0.5 * ( chi1 / (mass1*mass1) + chi2 / (mass2*mass2) );
        //chiA = 0.5 * ( chi1 / (mass1*mass1) - chi2 / (mass2*mass2) );
        
        XLALSimIMRSpinEOBCalculateSigmaStar( &sKerr, mass1, mass2, &s1Vec, &s2Vec );
        XLALSimIMRSpinEOBCalculateSigmaKerr( &sStar, mass1, mass2, &s1Vec, &s2Vec );
        // print_debug( "sigStar = (%g, %g, %g)\n", sStarData[0], sStarData[1], sStarData[2]);
        // print_debug( "sigKerr = (%g, %g, %g)\n", sStarData[0], sStarData[1], sStarData[2]);

        /* The a in the flux has been set to zero, but not in the Hamiltonian */
        a = sqrt(sKerr.data[0]*sKerr.data[0] + sKerr.data[1]*sKerr.data[1] + sKerr.data[2]*sKerr.data[2]);
        //XLALSimIMREOBCalcSpinFacWaveformCoefficients( params->eobParams->hCoeffs, mass1, mass2, eta, /*a*/0.0, chiS, chiA );
        //XLALSimIMRCalculateSpinEOBHCoeffs( params->seobCoeffs, eta, a );
        ham = EOBHamiltonian( eta, &qCartVec, &pCartVec, &s1VecNorm, &s2VecNorm, &sKerr, &sStar, params->tortoise, params->seobCoeffs );
        // print_debug( "hamiltonian at this point is %.16e\n", ham );
        
        /* And now, finally, the flux */
        REAL8Vector polarDynamics;
        REAL8       polarData[4];
        REAL8Vector cartDynamics;
        REAL8       cartData[6];
        
        polarDynamics.length = 4;
        polarDynamics.data = polarData;
        cartDynamics.length = 6;
        cartDynamics.data = cartData;
        
        polarData[0] = qSph[0];
        polarData[1] = 0.;
        polarData[2] = pSph[0];
        polarData[3] = pSph[2];
        memcpy(cartData, qCart, sizeof(qCart));
        memcpy(cartData+3, pCart, sizeof(pCart));
        // print_debug("polarData = (%g, %g, %g, %g)\n", 
        //     polarData[0], polarData[1], polarData[2], polarData[3]);
        // print_debug("ham = %g\n", ham);
        REAL8 tmpdvalues[4] = {0.,omega,0.,0.};
        flux  = XLALInspiralSpinFactorizedFlux( &polarDynamics, &cartDynamics, nqcCoeffs, omega, 0, qSph[0]*omega, params, ham, lMax);
        // flux  = InspiralSpinFactorizedFlux_elip( &polarDynamics, polarData, tmpdvalues, nqcCoeffs, omega, params, ham, lMax, eobVersion );
        flux  = flux / eta;
        if (ecc != 0.0)
            rDot = 0.0;
        else
            rDot  = - flux / dEdr;
        // print_debug("flux = %g, dEdr = %g, rDot = %g\n", flux, dEdr, rDot);
        /* We now need dHdpr - we take it that it is safely linear up to a pr of 1.0e-3 */
        cartValues[3] = 1.0e-3;
        dHdpr = XLALSpinHcapNumDerivWRTParam( 3, cartValues, params );
        /*printf( "Ingredients going into prDot:\n" );
         printf( "flux = %.16e, dEdr = %.16e, dHdpr = %.16e\n", flux, dEdr, dHdpr );*/
        
        /* We can now calculate what pr should be taking into account the flux */
        pSph[0] = rDot / (dHdpr / cartValues[3] );
    }
    else
    {
        /* Since d2Hdr2 has evaluated to zero, we cannot do the above. Just set pr to zero */
        //printf( "d2Hdr2 is zero!\n" );
        pSph[0] = 0;
    }
    
    /* Now we are done - convert back to cartesian coordinates ) */
    SphericalToCartesian( qCart, pCart, qSph, pSph );
    PRINT_LOG_INFO(LOG_DEBUG, "Sph initial condition : r = (%e,%e,%e), p = (%e,%e,%e)", qSph[0], qSph[1], qSph[2], pSph[0], pSph[1], pSph[2]);
    /* STEP 5) Rotate back to the original inertial frame by inverting the rotation of STEP 3 and then
     *         inverting the rotation of STEP 1.
     */
    
    /* Undo rotations to get back to the original basis */
    /* Second rotation */
    ApplyRotationMatrix( invMatrix2, rHat );
    ApplyRotationMatrix( invMatrix2, vHat );
    ApplyRotationMatrix( invMatrix2, LnHat );
    ApplyRotationMatrix( invMatrix2, tmpS1 );
    ApplyRotationMatrix( invMatrix2, tmpS2 );
    ApplyRotationMatrix( invMatrix2, tmpS1Norm );
    ApplyRotationMatrix( invMatrix2, tmpS2Norm );
    ApplyRotationMatrix( invMatrix2, qCart );
    ApplyRotationMatrix( invMatrix2, pCart );
    
    /* First rotation */
    ApplyRotationMatrix( invMatrix, rHat );
    ApplyRotationMatrix( invMatrix, vHat );
    ApplyRotationMatrix( invMatrix, LnHat );
    ApplyRotationMatrix( invMatrix, tmpS1 );
    ApplyRotationMatrix( invMatrix, tmpS2 );
    ApplyRotationMatrix( invMatrix, tmpS1Norm );
    ApplyRotationMatrix( invMatrix, tmpS2Norm );
    ApplyRotationMatrix( invMatrix, qCart );
    ApplyRotationMatrix( invMatrix, pCart );
    

    /* If required, apply the tortoise transform */
    if ( tmpTortoise )
    {
        REAL8 r = sqrt(qCart[0]*qCart[0] + qCart[1]*qCart[1] + qCart[2]*qCart[2] );
        REAL8 deltaR = XLALSimIMRSpinEOBHamiltonianDeltaR( params->seobCoeffs, r, eta, a );
        REAL8 deltaT = XLALSimIMRSpinEOBHamiltonianDeltaT( params->seobCoeffs, r, eta, a );
        REAL8 csi    = sqrt( deltaT * deltaR )/(r*r + a*a);
        
        REAL8 pr = (qCart[0]*pCart[0] + qCart[1]*pCart[1] + qCart[2]*pCart[2])/r;
        params->tortoise = tmpTortoise;
        
        //printf( "Applying the tortoise to p (csi = %.26e)\n", csi );
        
        for ( i = 0; i < 3; i++ )
        {
            pCart[i] = pCart[i] + qCart[i] * pr * (csi - 1.) / r;
        }
    }

    /* Now copy the initial conditions back to the return vector */
    memcpy( initConds->data, qCart, sizeof(qCart) );
    memcpy( initConds->data+3, pCart, sizeof(pCart) );
    memcpy( initConds->data+6, tmpS1Norm, sizeof(tmpS1Norm) );
    memcpy( initConds->data+9, tmpS2Norm, sizeof(tmpS2Norm) );
    
    gsl_matrix_free(rotMatrix2);
    gsl_matrix_free(invMatrix2);
    
    gsl_matrix_free(rotMatrix);
    gsl_matrix_free(invMatrix);
    
    //printf( "THE FINAL INITIAL CONDITIONS:\n");
    /*printf( " %.16e %.16e %.16e\n%.16e %.16e %.16e\n%.16e %.16e %.16e\n%.16e %.16e %.16e\n", initConds->data[0], initConds->data[1], initConds->data[2],
     initConds->data[3], initConds->data[4], initConds->data[5], initConds->data[6], initConds->data[7], initConds->data[8],
     initConds->data[9], initConds->data[10], initConds->data[11] );*/
    
    return CEV_SUCCESS;
}

static REAL8 diffHamiltonianBype(const REAL8 p,
                                 const REAL8 e,
                                 const REAL8 pphi,
                                 SpinEOBParams *params)
{
    REAL8 Hp, Hm;
    REAL8         cartValues[6];
    memset( cartValues, 0, sizeof( cartValues ) );
    HcapDerivParams hcdpms;
    hcdpms.params  = params;
    hcdpms.varyParam = 0;
    hcdpms.values = cartValues;
    REAL8 rp = p / (1. + e);
    REAL8 rm = p / (1. - e);
    cartValues[4] = pphi/rp; // pphi / r
    cartValues[4] = pphi/rm; // pphi / r
    Hp = GSLSpinAlignedHamiltonianWrapper_SA(rp, &hcdpms);
    Hm = GSLSpinAlignedHamiltonianWrapper_SA(rm, &hcdpms);
    return Hp - Hm;
}

typedef
struct tagEAnomalyParams
{
    SpinEOBParams *params;     /**<< Spin EOB parameters -- physical, pre-computed, etc. */
    REAL8          omega;      /**<< Orbital frequency */
    REAL8          e0;
}
EAnomalyParams;

INT CalculateAOmegaFromrpphi(REAL8 r, REAL8 pphi, SpinEOBParams *core,
    REAL8 *omegaOut);

static int
XLALFindSphericalOrbitSAEAnomaly(const gsl_vector *x, /**<< Parameters requested by gsl root finder */
                       void *params,        /**<< Spin EOB parameters */
                       gsl_vector *f        /**<< Function values for the given parameters */
)
{
    EAnomalyParams *rootParams = (EAnomalyParams *) params;
    
    REAL8 pp, pphi;
    
    /* Numerical derivative of Hamiltonian wrt given value */
    REAL8 dHdx, dHdpy;
    REAL8 dHdr, dHdpphi;
    REAL8 e0, omega0;
    /* Populate the appropriate values */
    /* In the special theta=pi/2 phi=0 case, r is x */
    pp  = gsl_vector_get( x, 0 );
    pphi = gsl_vector_get( x, 1 );
    e0 = rootParams->e0;
    omega0 = rootParams->omega;
    
    REAL8 diffH = diffHamiltonianBype(pp, e0, pphi, rootParams->params);
    REAL8 omegaAns, diffw = 0.0;
    if (CalculateAOmegaFromrpphi(pp / (1. + e0), pphi, rootParams->params, &omegaAns) != CEV_SUCCESS)
    {
        diffw = -100;
    } else 
        diffw = omegaAns - omega0;
    /* populate the function vector */
    gsl_vector_set( f, 0, diffH);
    gsl_vector_set( f, 1, diffw);
    
    //printf( "Current funcvals = %.16e %.16e %.16e\n", gsl_vector_get( f, 0 ), gsl_vector_get( f, 1 ),
    //  gsl_vector_get( f, 2 )/*dHdpphi*/ );
    
    return CEV_SUCCESS;
}

int calc_rprpphi_from_e_anomaly(REAL8 e, REAL8 zeta, REAL8 omega, SpinEOBParams *params, REAL8 *r, REAL8 *pr, REAL8 *pphi)
{
    EAnomalyParams rootParams;
    const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver *rootSolver = NULL;

    gsl_multiroot_function rootFunction;
    gsl_vector *initValues  = NULL;
    gsl_vector *finalValues = NULL;
    int gslStatus;
    const int maxIter = 100;

    REAL8 v0    = GET_CBRT( omega );
    rootParams.e0 = e;
    rootParams.omega  = omega;
    rootParams.params = params;
    rootSolver = gsl_multiroot_fsolver_alloc( T, 2 );
    if ( !rootSolver )
    {
        return CEV_FAILURE;
    }
    
    initValues = gsl_vector_calloc( 2 );
    if ( !initValues )
    {
        gsl_multiroot_fsolver_free( rootSolver );
        return CEV_FAILURE;
    }
    REAL8 p0 = pow(omega, -2./3.);
    gsl_vector_set( initValues, 0, p0 );
    gsl_vector_set( initValues, 1, sqrt(p0) ); // py
    rootFunction.f      = XLALFindSphericalOrbitSAEAnomaly;
    rootFunction.n      = 2;
    rootFunction.params = &rootParams;
    gsl_multiroot_fsolver_set( rootSolver, &rootFunction, initValues );

    int i = 0;
    do
    {
        gslStatus = gsl_multiroot_fsolver_iterate( rootSolver );
        if ( gslStatus != GSL_SUCCESS )
        {
            print_warning( "Error in GSL iteration function!\n" );
            gsl_multiroot_fsolver_free( rootSolver );
            gsl_vector_free( initValues );
            return CEV_FAILURE;
        }
        gslStatus = gsl_multiroot_test_residual( rootSolver->f, 1.0e-10 );
        i++;
    }
    while ( gslStatus == GSL_CONTINUE && i <= maxIter );

    if ( i > maxIter && gslStatus != GSL_SUCCESS )
    {
        gsl_multiroot_fsolver_free( rootSolver );
        gsl_vector_free( initValues );
        return CEV_FAILURE;
    }
    finalValues = gsl_multiroot_fsolver_root( rootSolver );
    REAL8 psol, Lsol;
    psol = gsl_vector_get( finalValues, 0 );
    Lsol = gsl_vector_get( finalValues, 1 );
    gsl_multiroot_fsolver_free( rootSolver );
    gsl_vector_free( initValues );

    *r = psol / (1. - e);
    *pr = 0.0;
    *pphi = Lsol;
    return CEV_SUCCESS;
}




typedef
struct tagSEOBRootParamsV2
{
    REAL8          values[12]; /**<< Dynamical variables, x, y, z, px, py, pz, S1x, S1y, S1z, S2x, S2y and S2z */
    SpinEOBParams *params;     /**<< Spin EOB parameters -- physical, pre-computed, etc. */
    REAL8          omega;      /**<< Orbital frequency */
    REAL8          e0;
    REAL8          sz;
    REAL8          cz;
}
SEOBRootParamsV2;

void calculate_prDot_from_ezetapphi(REAL8 eta, REAL8 chi1, REAL8 chi2, 
        REAL8 e, REAL8 sz, REAL8 cz, REAL8 pf,
        REAL8 *ret_prT, REAL8 *ret_prDot);
REAL8 calculate_rDot_from_pr(REAL8 eta, REAL8 chi1, REAL8 chi2, REAL8 r, REAL8 prT, REAL8 pf);

static int
XLALFindSphericalOrbitSAWithAnomaly(const gsl_vector *x, /**<< Parameters requested by gsl root finder */
                       void *params,        /**<< Spin EOB parameters */
                       gsl_vector *f        /**<< Function values for the given parameters */
)
{
    SEOBRootParamsV2 *rootParams = (SEOBRootParamsV2 *) params;
    
    REAL8 r, py;
    REAL8 pphi, prT, prDot;

    /* Numerical derivative of Hamiltonian wrt given value */
    REAL8 dHdx, dHdpy;
    REAL8 dHdr, dHdpphi;
    REAL8 e0;
    /* Populate the appropriate values */
    /* In the special theta=pi/2 phi=0 case, r is x */
    rootParams->values[0] = r  = gsl_vector_get( x, 0 );
    rootParams->values[4] = py = gsl_vector_get( x, 1 );
    e0 = rootParams->e0;
    pphi   = r * py;
    // print_debug( "Values r = %.16e, py = %.16e\n", r, py );
    calculate_prDot_from_ezetapphi(rootParams->params->eta, rootParams->params->chi1, rootParams->params->chi2,
        e0, rootParams->sz, rootParams->cz, pphi, &prT, &prDot);
    
    rootParams->values[3] = prT;
    
    /* dHdR */
    dHdx = XLALSpinHcapNumDerivWRTParam( 0, rootParams->values, rootParams->params );
    if ( IS_REAL8_FAIL_NAN( dHdx ) )
    {
        return CEV_FAILURE;
    }
    //printf( "dHdx = %.16e\n", dHdx );
    
    /* dHdPphi (I think we can use dHdPy in this coord system) */
    /* TODO: Check this is okay */
    dHdpy = XLALSpinHcapNumDerivWRTParam( 4, rootParams->values, rootParams->params );
    if ( IS_REAL8_FAIL_NAN( dHdpy ) )
    {
        return CEV_FAILURE;
    }
        
    /* Now convert to spherical polars */
    dHdr      = dHdx - dHdpy * pphi / (r*r);
    dHdpphi   = dHdpy / r;
    
    /* populate the function vector */
    gsl_vector_set( f, 0, dHdr + prDot);
    gsl_vector_set( f, 1, dHdpphi - rootParams->omega );
    // print_debug( "pr, prDot = %.16e %.16e\n\n", prT, prDot);
    // print_debug( "Current funcvals = %.16e %.16e\n", gsl_vector_get( f, 0 ), gsl_vector_get( f, 1 ));
    //  gsl_vector_get( f, 2 )/*dHdpphi*/ );
    
    return CEV_SUCCESS;
}

typedef struct {
    SpinEOBParams *eobpms;
    REAL8 r0;
    REAL8 pphi0;
    REAL8 omega0;
    REAL8 sphvalues[12];
    REAL8 cartvalues[12];
}RootParams_SolPr;

// calculate dH/dpr - rDot == 0
static REAL8 EOBSolvingInitialPr_from_ezeta(REAL8 pr, void *params)
{
    RootParams_SolPr *fparams = (RootParams_SolPr *) params;
    SpinEOBParams *eobpms = fparams->eobpms;
    REAL8 r0, pphi0, omega0;
    REAL8 sphValues[12];
    REAL8 cartValues[12];
    r0 = fparams->r0;
    pphi0 = fparams->pphi0;
    omega0 = fparams->omega0;
    memcpy(sphValues, fparams->sphvalues, sizeof(sphValues));
    memcpy(cartValues, fparams->cartvalues, sizeof(cartValues));

    REAL8 eq;
    REAL8 rDot0, rDot1, dHdpr;
    REAL8 d2Hdrdpphi, d2Hdr2;
    sphValues[3] = pr;
    cartValues[3] = pr;
    // CODING PLACE
    d2Hdr2 = XLALCalculateSphHamiltonianDeriv2( 0, 0, sphValues, fparams->eobpms );
    d2Hdrdpphi = XLALCalculateSphHamiltonianDeriv2( 0, 5, sphValues, fparams->eobpms );
    REAL8Vector qCartVec, pCartVec;
    REAL8 qCart[3], pCart[3];
    qCartVec.length = 3;
    pCartVec.length = 3;
    qCartVec.data = qCart;
    pCartVec.data = pCart;
    memcpy(qCart, cartValues, sizeof(qCart));
    memcpy(pCart, 3+cartValues, sizeof(qCart));
    REAL8 ham;
    ham = EOBHamiltonian( fparams->eobpms->eta, &qCartVec, &pCartVec, 
            fparams->eobpms->s1Vec, fparams->eobpms->s2Vec, 
            fparams->eobpms->sigmaKerr, fparams->eobpms->sigmaStar, 
            1, fparams->eobpms->seobCoeffs);
    // print_debug( "hamiltonian at this point is %.16e\n", ham );
    
    /* And now, finally, the flux */
    REAL8Vector polarDynamics;
    REAL8       polarData[4];
    REAL8Vector cartDynamics;
    REAL8       cartData[6];
    
    polarDynamics.length = 4;
    polarDynamics.data = polarData;
    cartDynamics.length = 6;
    cartDynamics.data = cartData;
    
    polarData[0] = r0;
    polarData[1] = 0.;
    polarData[2] = pr;
    polarData[3] = pphi0;
    // memcpy(cartData, qCart, sizeof(qCart));
    // memcpy(cartData+3, pCart, sizeof(pCart));
    memcpy(cartData, cartValues, sizeof(cartData));
    // print_debug("polarData = (%g, %g, %g, %g)\n", 
    //     polarData[0], polarData[1], polarData[2], polarData[3]);
    // print_debug("ham = %g\n", ham);
    REAL8 flux;
    flux  = XLALInspiralSpinFactorizedFlux( &polarDynamics, &cartDynamics, eobpms->nqcCoeffs, omega0, 0, r0*omega0, eobpms, ham, 8);
    // flux  = InspiralSpinFactorizedFlux_elip( &polarDynamics, polarData, tmpdvalues, nqcCoeffs, omega, params, ham, lMax, eobVersion );
    flux  = flux / eobpms->eta;

    rDot0 = calculate_rDot_from_pr(eobpms->eta, eobpms->chi1, eobpms->chi2, r0, pr, pphi0);
    rDot1 = -flux * d2Hdrdpphi / (d2Hdr2 * omega0);
    dHdpr = XLALSpinHcapNumDerivWRTParam( 3, cartValues, eobpms );
    eq = dHdpr - (rDot0 + rDot1);
    // print_debug("d2Hdrdpphi = %.16e\n", d2Hdrdpphi);
    // print_debug("d2Hdr2 = %.16e\n", d2Hdr2);
    // print_debug("pr, ham, flux, d2Hdrdpphi, d2Hdr2, rDot0, rDot1, dHdpr = \n\t%.16e, %.16e, %.16e\n\t%.16e, %.16e\n\t%.16e, %.16e, %.16e\n",
    //     pr, ham, flux, 
    //     d2Hdrdpphi, d2Hdr2, 
    //     rDot0, rDot1, dHdpr);
    return 100*eq;
}

int calc_rpphi_from_eanomaly(REAL8 e, REAL8 anomaly, REAL8 omega, SpinEOBParams *params, REAL8 *r, REAL8 *pr, REAL8 *pphi)
{
    SEOBRootParamsV2 rootParams;
    const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver *rootSolver = NULL;

    gsl_multiroot_function rootFunction;
    gsl_vector *initValues  = NULL;
    gsl_vector *finalValues = NULL;
    int gslStatus;
    const int maxIter = 100;

    REAL8 v0    = GET_CBRT( omega );
    rootParams.e0 = e;
    rootParams.omega  = omega;
    rootParams.sz = sin(anomaly);
    rootParams.cz = cos(anomaly);
    // print_debug("zeta sz, cz = %.16e, %.16e, %.16e\n", anomaly, rootParams.sz, rootParams.cz);
    rootParams.params = params;
    memset(rootParams.values, 0, sizeof(rootParams.values));
    rootParams.values[0] = (1.-e) / (v0*v0);
    rootParams.values[4] = pow(1.-e, 0.5) * v0;
    REAL8 mT2 = (params->m1 + params->m2) * (params->m1 + params->m2);
    for (int i=0; i<3; i++)
    {
        rootParams.values[i+6] = params->s1Vec->data[i]*mT2;
        rootParams.values[i+9] = params->s2Vec->data[i]*mT2;
    }
    rootSolver = gsl_multiroot_fsolver_alloc( T, 2 );
    if ( !rootSolver )
    {
        return CEV_FAILURE;
    }
    
    initValues = gsl_vector_calloc( 2 );
    if ( !initValues )
    {
        gsl_multiroot_fsolver_free( rootSolver );
        return CEV_FAILURE;
    }
    gsl_vector_set( initValues, 0, rootParams.values[0] ); // r
    gsl_vector_set( initValues, 1, rootParams.values[4] ); // py
    rootFunction.f      = XLALFindSphericalOrbitSAWithAnomaly;
    rootFunction.n      = 2;
    rootFunction.params = &rootParams;
    gsl_multiroot_fsolver_set( rootSolver, &rootFunction, initValues );

    int i = 0;
    do
    {
        gslStatus = gsl_multiroot_fsolver_iterate( rootSolver );
        if ( gslStatus != GSL_SUCCESS )
        {
            print_warning( "Error in GSL iteration function!\n" );
            gsl_multiroot_fsolver_free( rootSolver );
            gsl_vector_free( initValues );
            return CEV_FAILURE;
        }
        gslStatus = gsl_multiroot_test_residual( rootSolver->f, 1.0e-10 );
        i++;
    }
    while ( gslStatus == GSL_CONTINUE && i <= maxIter );

    if ( i > maxIter && gslStatus != GSL_SUCCESS )
    {
        gsl_multiroot_fsolver_free( rootSolver );
        gsl_vector_free( initValues );
        return CEV_FAILURE;
    }
    finalValues = gsl_multiroot_fsolver_root( rootSolver );
    REAL8 r0, pphi0;
    *r = r0 = gsl_vector_get( finalValues, 0 );
    *pphi = pphi0 = gsl_vector_get( finalValues, 1 ) * gsl_vector_get( finalValues, 0 );
    gsl_multiroot_fsolver_free( rootSolver );
    gsl_vector_free( initValues );

    // Step 2, solve initial pr
    REAL8 pr0, prDot0;
    calculate_prDot_from_ezetapphi(params->eta, params->chi1, params->chi2,
        e, rootParams.sz, rootParams.cz, pphi0, &pr0, &prDot0);
    const gsl_root_fsolver_type *TT = gsl_root_fsolver_brent;
    gsl_root_fsolver *rootSolver1D = gsl_root_fsolver_alloc (TT);
    gsl_function F;
    RootParams_SolPr fparams2;
    fparams2.eobpms = params;
    fparams2.r0 = r0;
    fparams2.pphi0 = pphi0;
    fparams2.omega0 = omega;

    REAL8 qCart[3], pCart[3];
    REAL8 qSph[3], pSph[3];
    REAL8 tmpS1[3], tmpS2[3];
    REAL8 cartValues[12];
    REAL8 sphValues[12];
    memset(qSph, 0, sizeof(qSph));
    memset(pSph, 0, sizeof(pSph));
    memset(qCart, 0, sizeof(qCart));
    memset(pCart, 0, sizeof(pCart));
    qCart[0] = r0;
    pCart[0] = pr0;
    pCart[1] = pphi0/r0;

    CartesianToSpherical( qSph, pSph, qCart, pCart );
    for (int i=0; i<3; i++)
    {
        tmpS1[i] = params->s1Vec->data[i]*mT2;
        tmpS2[i] = params->s2Vec->data[i]*mT2;
    }
    memcpy( sphValues, qSph, sizeof( qSph ) );
    memcpy( sphValues+3, pSph, sizeof( pSph ) );
    memcpy( sphValues+6, tmpS1, sizeof(tmpS1) );
    memcpy( sphValues+9, tmpS2, sizeof(tmpS2) );
    
    memcpy( cartValues, qCart, sizeof(qCart) );
    memcpy( cartValues+3, pCart, sizeof(pCart) );
    memcpy( cartValues+6, tmpS1, sizeof(tmpS1) );
    memcpy( cartValues+9, tmpS2, sizeof(tmpS2) );
    memcpy(fparams2.sphvalues, sphValues, sizeof(sphValues));
    memcpy(fparams2.cartvalues, cartValues, sizeof(cartValues));
    F.function = &EOBSolvingInitialPr_from_ezeta;
    F.params = &fparams2;
    // print_debug("tmpeq(0) = %.16e\n", EOBSolvingInitialPr_from_ezeta(0, &fparams2));
    // print_debug("tmpeq(pr0) = %.16e\n", EOBSolvingInitialPr_from_ezeta(pr0, &fparams2));

    /* Given this, we can start to calculate the initial conditions */
    /* for spherical coords in the new basis */
    REAL8 x_lo, x_hi, root;
    x_lo = GET_MAX(-0.9, pr0 - 0.2);
    x_hi = GET_MIN(pr0 + 0.2, 0.9);
    gsl_root_fsolver_set (rootSolver1D, &F, x_lo, x_hi);
    // print_debug("x_lo = %.16e, x_hi = %.16e\n", x_lo, x_hi);
    // print_debug("tmpeq(x_lo) = %.16e, tmpeq(x_hi) = %.16e\n", 
    //     EOBSolvingInitialPr_from_ezeta(x_lo, &fparams2), 
    //     EOBSolvingInitialPr_from_ezeta(x_hi, &fparams2));
#if 1
    /* Initialise the gsl stuff */
    GSL_START;
    INT status;
    /* We are now ready to iterate to find the solution */
    i = 0;
    do
    {
        gslStatus = gsl_root_fsolver_iterate( rootSolver1D );
        if ( gslStatus != GSL_SUCCESS )
        {
            // PRINT_LOG_INFO(LOG_WARNING, "cannot find initial condition for pr");
            print_debug("cannot find initial condition for pr0\n");
            gsl_root_fsolver_free( rootSolver1D );
            *pr = pr0;
            return CEV_SUCCESS;
        }
        x_lo = gsl_root_fsolver_x_lower (rootSolver1D);
        x_hi = gsl_root_fsolver_x_upper (rootSolver1D);
        gslStatus = gsl_root_test_interval (x_lo, x_hi, 1e-9, 1e-9);
        // print_debug ("%5d [%.7f, %.7f] %.7f %+.7f\n",
        //     i, x_lo, x_hi,
        //     gsl_root_fsolver_root (rootSolver1D),
        //     x_hi - x_lo);
        i++;
    }
    while ( gslStatus == GSL_CONTINUE && i <= maxIter );

    if ( i > maxIter && gslStatus != GSL_SUCCESS )
    {
        // PRINT_LOG_INFO(LOG_WARNING, "cannot find initial condition for pr");
        print_debug("cannot find initial condition for pr0\n");
        root = pr0;
    }
    else 
        root = gsl_root_fsolver_root (rootSolver1D);
    *pr = root;
    GSL_END;
    // print_log( "Spherical orbit conditions here given by the following:\n" );
    //  print_err( " x = %.16e, py = %.16e, pz = %.16e\n", gsl_vector_get( finalValues, 0 ),
    //  gsl_vector_get( finalValues, 1 ), gsl_vector_get( finalValues, 2 ) );
    // print_debug("pr0 = %.16e, eq(pr0) = %.16e\n", pr0, EOBSolvingInitialPr_from_ezeta(pr0, &fparams2));
    // print_debug("root = %.16e, eq(root) = %.16e\n", root, EOBSolvingInitialPr_from_ezeta(root, &fparams2));
#else
    *pr = pr0;
#endif

    gsl_root_fsolver_free(rootSolver1D);
    return CEV_SUCCESS;
}

/* For nonprecessing case */
INT EOBInitialConditionsSA_e_anomaly(REAL8Vector    *initConds,
                         const REAL8    mass1,
                         const REAL8    mass2,
                         const REAL8    fMin,
                         const REAL8    ecc,
                         const REAL8    anomaly,
                         const REAL8    inc,
                         const REAL8    spin1[],
                         const REAL8    spin2[],
                         SpinEOBParams  *params)
{
    if (!initConds)
        return CEV_FAILURE;
    static const int lMax = 8;
    INT i;
    int tmpTortoise;
    /* non-zero eccentricity - Start frequency correction */
    // REAL8 ecc = params->eccentricity;
    
    REAL8 mTotal;
    REAL8 eta;
    REAL8 omega, v0;   /* Initial velocity and angular frequency */
    
    REAL8 ham;      /* Hamiltonian */
    
    REAL8 LnHat[3]; /* Initial orientation of angular momentum */
    REAL8 rHat[3];  /* Initial orientation of radial vector */
    REAL8 vHat[3];  /* Initial orientation of velocity vector */
    REAL8 Lhat[3];  /* Direction of relativistic ang mom */
    REAL8 qHat[3];
    REAL8 pHat[3];
    
    /* q and p vectors in Cartesian and spherical coords */
    REAL8 qCart[3], pCart[3];
    REAL8 qSph[3], pSph[3];
    
    /* We will need to manipulate the spin vectors */
    /* We will use temporary vectors to do this */
    REAL8 tmpS1[3];
    REAL8 tmpS2[3];
    REAL8 tmpS1Norm[3];
    REAL8 tmpS2Norm[3];
    
    REAL8Vector qCartVec, pCartVec;
    REAL8Vector s1Vec, s2Vec, s1VecNorm, s2VecNorm;
    REAL8Vector sKerr, sStar;
    REAL8       sKerrData[3], sStarData[3];
    REAL8       a = 0.; //, chiS, chiA;
    //REAL8       chi1, chi2;
    
    /* We will need a full values vector for calculating derivs of Hamiltonian */
    REAL8 sphValues[12];
    REAL8 cartValues[12];
    
    /* Matrices for rotating to the new basis set. */
    /* It is more convenient to calculate the ICs in a simpler basis */
    gsl_matrix *rotMatrix  = NULL;
    gsl_matrix *invMatrix  = NULL;
    gsl_matrix *rotMatrix2 = NULL;
    gsl_matrix *invMatrix2 = NULL;
    
    /* Root finding stuff for finding the spherical orbit */
    EgwRootParamsV2 rootParams;
    const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
    gsl_root_fsolver *rootSolver = gsl_root_fsolver_alloc (T);
    gsl_function F;

    int gslStatus;
    const int maxIter = 100;
    
    memset( &rootParams, 0, sizeof( rootParams ) );
    
    mTotal = mass1 + mass2;
    eta    = mass1 * mass2 / (mTotal * mTotal);
    memcpy( tmpS1, spin1, sizeof(tmpS1) );
    memcpy( tmpS2, spin2, sizeof(tmpS2) );
    memcpy( tmpS1Norm, spin1, sizeof(tmpS1Norm) );
    memcpy( tmpS2Norm, spin2, sizeof(tmpS2Norm) );
    for ( i = 0; i < 3; i++ )
    {
        tmpS1Norm[i] /= mTotal * mTotal;
        tmpS2Norm[i] /= mTotal * mTotal;
    }
    // eobVersion = params->seobCoeffs->eobVersion;
    /* We compute the ICs for the non-tortoise p, and convert at the end */
    tmpTortoise      = params->tortoise;
    params->tortoise = 0;
    
    EOBNonQCCoeffs *nqcCoeffs = NULL;
    nqcCoeffs = params->nqcCoeffs;
    
    /* STEP 1) Rotate to LNhat0 along z-axis and N0 along x-axis frame, where LNhat0 and N0 are initial normal to
     *         orbital plane and initial orbital separation;
     */
    
    /* Set the initial orbital ang mom direction. Taken from STPN code */
    LnHat[0] = GET_SIN(inc);
    LnHat[1] = 0.;
    LnHat[2] = GET_COS(inc);
    
    /* Set the radial direction - need to take care to avoid singularity if L is along z axis */
    if ( LnHat[2] > 0.9999 )
    {
        rHat[0] = 1.;
        rHat[1] = rHat[2] = 0.;
    }
    else
    {
        REAL8 theta0 = GET_ATAN( - LnHat[2] / LnHat[0] ); /* theta0 is between 0 and Pi */
        rHat[0] = GET_SIN( theta0 );
        rHat[1] = 0;
        rHat[2] = GET_COS( theta0 );
    }
    
    /* Now we can complete the triad */
    vHat[0] = CalculateCrossProduct( 0, LnHat, rHat );
    vHat[1] = CalculateCrossProduct( 1, LnHat, rHat );
    vHat[2] = CalculateCrossProduct( 2, LnHat, rHat );
    
    NormalizeVector( vHat );
    
    /* XXX Test code XXX */
    /*for ( i = 0; i < 3; i++ )
     {
     printf ( " LnHat[%d] = %.16e, rHat[%d] = %.16e, vHat[%d] = %.16e\n", i, LnHat[i], i, rHat[i], i, vHat[i] );
     }
     
     printf("\n\n" );
     for ( i = 0; i < 3; i++ )
     {
     printf ( " s1[%d] = %.16e, s2[%d] = %.16e\n", i, tmpS1[i], i, tmpS2[i] );
     }*/
    
    /* Allocate and compute the rotation matrices */
    rotMatrix = gsl_matrix_alloc( 3, 3 );
    invMatrix = gsl_matrix_alloc( 3, 3 );
    if ( !rotMatrix || !invMatrix )
    {
        if ( rotMatrix ) gsl_matrix_free( rotMatrix );
        if ( invMatrix ) gsl_matrix_free( invMatrix );
        return CEV_FAILURE;
    }
    
    if ( CalculateRotationMatrix( rotMatrix, invMatrix, rHat, vHat, LnHat ) == CEV_FAILURE )
    {
        gsl_matrix_free( rotMatrix );
        gsl_matrix_free( invMatrix );
        return CEV_FAILURE;
    }
    
    /* Rotate the orbital vectors and spins */
    ApplyRotationMatrix( rotMatrix, rHat );
    ApplyRotationMatrix( rotMatrix, vHat );
    ApplyRotationMatrix( rotMatrix, LnHat );
    ApplyRotationMatrix( rotMatrix, tmpS1 );
    ApplyRotationMatrix( rotMatrix, tmpS2 );
    ApplyRotationMatrix( rotMatrix, tmpS1Norm );
    ApplyRotationMatrix( rotMatrix, tmpS2Norm );
    
    /* XXX Test code XXX */
    /*printf( "\nAfter applying rotation matrix:\n\n" );
     for ( i = 0; i < 3; i++ )
     {
     printf ( " LnHat[%d] = %.16e, rHat[%d] = %.16e, vHat[%d] = %.16e\n", i, LnHat[i], i, rHat[i], i, vHat[i] );
     }
     
     printf("\n\n" );
     for ( i = 0; i < 3; i++ )
     {
     printf ( " s1[%d] = %.16e, s2[%d] = %.16e\n", i, tmpS1[i], i, tmpS2[i] );
     }*/
    
    /* STEP 2) After rotation in STEP 1, in spherical coordinates, phi0 and theta0 are given directly in Eq. (4.7),
     *         r0, pr0, ptheta0 and pphi0 are obtained by solving Eqs. (4.8) and (4.9) (using gsl_multiroot_fsolver).
     *         At this step, we find initial conditions for a spherical orbit without radiation reaction.
     */
    // X.L.: here we solve equation (e22, omegaMinus)(rp, rm)
    /* Calculate the initial velocity from the given initial frequency */
    // initial orbital angular velocity
    omega = CST_PI * mTotal * CST_MTSUN_SI * fMin;
    memset( qCart, 0, sizeof(qCart) );
    memset( pCart, 0, sizeof(pCart) );
    // calc_rpphi_from_e(root, omega, params, &(qCart[0]), &(pCart[1]));
    calc_rpphi_from_eanomaly(ecc, anomaly, omega, params, &(qCart[0]), &(pCart[0]), &(pCart[1]));
    // qCart[0] = gsl_vector_get( finalValues, 0 );
    // pCart[1] = gsl_vector_get( finalValues, 1 );
    pCart[1] /= qCart[0];
    /* Free the GSL root finder, since we're done with it */
    gsl_root_fsolver_free( rootSolver );
    
    /* STEP 3) Rotate to L0 along z-axis and N0 along x-axis frame, where L0 is the initial orbital angular momentum
     *         and L0 is calculated using initial position and linear momentum obtained in STEP 2.
     */
    
    /* Now we can calculate the relativistic L and rotate to a new basis */
    memcpy( qHat, qCart, sizeof(qCart) );
    memcpy( pHat, pCart, sizeof(pCart) );
    
    NormalizeVector( qHat );
    NormalizeVector( pHat );
    
    Lhat[0] = CalculateCrossProduct( 0, qHat, pHat );
    Lhat[1] = CalculateCrossProduct( 1, qHat, pHat );
    Lhat[2] = CalculateCrossProduct( 2, qHat, pHat );
    
    NormalizeVector( Lhat );
    
    rotMatrix2 = gsl_matrix_alloc( 3, 3 );
    invMatrix2 = gsl_matrix_alloc( 3, 3 );
    
    if ( CalculateRotationMatrix( rotMatrix2, invMatrix2, qHat, pHat, Lhat ) == CEV_FAILURE )
    {
        gsl_matrix_free( rotMatrix );
        gsl_matrix_free( invMatrix );
        return CEV_FAILURE;
    }
    
    ApplyRotationMatrix( rotMatrix2, rHat );
    ApplyRotationMatrix( rotMatrix2, vHat );
    ApplyRotationMatrix( rotMatrix2, LnHat );
    ApplyRotationMatrix( rotMatrix2, tmpS1 );
    ApplyRotationMatrix( rotMatrix2, tmpS2 );
    ApplyRotationMatrix( rotMatrix2, tmpS1Norm );
    ApplyRotationMatrix( rotMatrix2, tmpS2Norm );
    ApplyRotationMatrix( rotMatrix2, qCart );
    ApplyRotationMatrix( rotMatrix2, pCart );
    
    /* STEP 4) In the L0-N0 frame, we calculate (dE/dr)|sph using Eq. (4.14), then initial dr/dt using Eq. (4.10),
     *         and finally pr0 using Eq. (4.15).
     */
    
    /* Now we can calculate the flux. Change to spherical co-ords */
    CartesianToSpherical( qSph, pSph, qCart, pCart );
    memcpy( sphValues, qSph, sizeof( qSph ) );
    memcpy( sphValues+3, pSph, sizeof( pSph ) );
    memcpy( sphValues+6, tmpS1, sizeof(tmpS1) );
    memcpy( sphValues+9, tmpS2, sizeof(tmpS2) );
    
    memcpy( cartValues, qCart, sizeof(qCart) );
    memcpy( cartValues+3, pCart, sizeof(pCart) );
    memcpy( cartValues+6, tmpS1, sizeof(tmpS1) );
    memcpy( cartValues+9, tmpS2, sizeof(tmpS2) );

#if 0
    REAL8 dHdpphi, d2Hdr2, d2Hdrdpphi;
    REAL8 rDot, dHdpr, flux, dEdr;
    
    d2Hdr2 = XLALCalculateSphHamiltonianDeriv2( 0, 0, sphValues, params );
    d2Hdrdpphi = XLALCalculateSphHamiltonianDeriv2( 0, 5, sphValues, params );
    dHdpphi = XLALSpinHcapNumDerivWRTParam( 4, cartValues, params ) / sphValues[0];
    dEdr  = - dHdpphi * d2Hdr2 / d2Hdrdpphi;
        
    if ( d2Hdr2 != 0.0 && ecc==0.0 )
    {
        /* We will need to calculate the Hamiltonian to get the flux */
        s1Vec.length = s2Vec.length = s1VecNorm.length = s2VecNorm.length = sKerr.length = sStar.length = 3;
        s1Vec.data = tmpS1;
        s2Vec.data = tmpS2;
        s1VecNorm.data = tmpS1Norm;
        s2VecNorm.data = tmpS2Norm;
        sKerr.data = sKerrData;
        sStar.data = sStarData;
        
        qCartVec.length = pCartVec.length = 3;
        qCartVec.data   = qCart;
        pCartVec.data   = pCart;
        
        //chi1 = tmpS1[0]*LnHat[0] + tmpS1[1]*LnHat[1] + tmpS1[2]*LnHat[2];
        //chi2 = tmpS2[0]*LnHat[0] + tmpS2[1]*LnHat[1] + tmpS2[2]*LnHat[2];
        
        // print_debug( "m1 = %g, m2 = %g\n", mass1, mass2 );
        // print_debug( "eta = %g\n", eta);
        // print_debug( "qCartVec = (%g, %g, %g)\n", qCart[0], qCart[1], qCart[2]);
        // print_debug( "pCartVec = (%g, %g, %g)\n", pCart[0], pCart[1], pCart[2]);
        // print_debug( "s1Vec = (%g, %g, %g)\n", tmpS1[0], tmpS1[1], tmpS1[2]);
        // print_debug( "s2Vec = (%g, %g, %g)\n", tmpS2[0], tmpS2[1], tmpS2[2]);
        
        //chiS = 0.5 * ( chi1 / (mass1*mass1) + chi2 / (mass2*mass2) );
        //chiA = 0.5 * ( chi1 / (mass1*mass1) - chi2 / (mass2*mass2) );
        
        XLALSimIMRSpinEOBCalculateSigmaStar( &sKerr, mass1, mass2, &s1Vec, &s2Vec );
        XLALSimIMRSpinEOBCalculateSigmaKerr( &sStar, mass1, mass2, &s1Vec, &s2Vec );
        // print_debug( "sigStar = (%g, %g, %g)\n", sStarData[0], sStarData[1], sStarData[2]);
        // print_debug( "sigKerr = (%g, %g, %g)\n", sStarData[0], sStarData[1], sStarData[2]);

        /* The a in the flux has been set to zero, but not in the Hamiltonian */
        a = sqrt(sKerr.data[0]*sKerr.data[0] + sKerr.data[1]*sKerr.data[1] + sKerr.data[2]*sKerr.data[2]);
        //XLALSimIMREOBCalcSpinFacWaveformCoefficients( params->eobParams->hCoeffs, mass1, mass2, eta, /*a*/0.0, chiS, chiA );
        //XLALSimIMRCalculateSpinEOBHCoeffs( params->seobCoeffs, eta, a );
        ham = EOBHamiltonian( eta, &qCartVec, &pCartVec, &s1VecNorm, &s2VecNorm, &sKerr, &sStar, params->tortoise, params->seobCoeffs );
        // print_debug( "hamiltonian at this point is %.16e\n", ham );
        
        /* And now, finally, the flux */
        REAL8Vector polarDynamics;
        REAL8       polarData[4];
        REAL8Vector cartDynamics;
        REAL8       cartData[6];
        
        polarDynamics.length = 4;
        polarDynamics.data = polarData;
        cartDynamics.length = 6;
        cartDynamics.data = cartData;
        
        polarData[0] = qSph[0];
        polarData[1] = 0.;
        polarData[2] = pSph[0];
        polarData[3] = pSph[2];
        memcpy(cartData, qCart, sizeof(qCart));
        memcpy(cartData+3, pCart, sizeof(pCart));
        // print_debug("polarData = (%g, %g, %g, %g)\n", 
        //     polarData[0], polarData[1], polarData[2], polarData[3]);
        // print_debug("ham = %g\n", ham);
        REAL8 tmpdvalues[4] = {0.,omega,0.,0.};
        flux  = XLALInspiralSpinFactorizedFlux( &polarDynamics, &cartDynamics, nqcCoeffs, omega, 0, qSph[0]*omega, params, ham, lMax);
        // flux  = InspiralSpinFactorizedFlux_elip( &polarDynamics, polarData, tmpdvalues, nqcCoeffs, omega, params, ham, lMax, eobVersion );
        flux  = flux / eta;
        if (ecc != 0.0)
            rDot = 0.0;
        else
            rDot  = - flux / dEdr;
        // print_debug("flux = %g, dEdr = %g, rDot = %g\n", flux, dEdr, rDot);
        /* We now need dHdpr - we take it that it is safely linear up to a pr of 1.0e-3 */
        cartValues[3] = 1.0e-3;
        dHdpr = XLALSpinHcapNumDerivWRTParam( 3, cartValues, params );
        /*printf( "Ingredients going into prDot:\n" );
         printf( "flux = %.16e, dEdr = %.16e, dHdpr = %.16e\n", flux, dEdr, dHdpr );*/
        
        /* We can now calculate what pr should be taking into account the flux */
        pSph[0] = rDot / (dHdpr / cartValues[3] );
    }
    else
    {
        /* Since d2Hdr2 has evaluated to zero, we cannot do the above. Just set pr to zero */
        //printf( "d2Hdr2 is zero!\n" );
        pSph[0] = 0;
    }
#endif
    /* Now we are done - convert back to cartesian coordinates ) */
    SphericalToCartesian( qCart, pCart, qSph, pSph );
    PRINT_LOG_INFO(LOG_DEBUG, "Sph initial condition : r = (%e,%e,%e), p = (%e,%e,%e)", qSph[0], qSph[1], qSph[2], pSph[0], pSph[1], pSph[2]);
    /* STEP 5) Rotate back to the original inertial frame by inverting the rotation of STEP 3 and then
     *         inverting the rotation of STEP 1.
     */
    
    /* Undo rotations to get back to the original basis */
    /* Second rotation */
    ApplyRotationMatrix( invMatrix2, rHat );
    ApplyRotationMatrix( invMatrix2, vHat );
    ApplyRotationMatrix( invMatrix2, LnHat );
    ApplyRotationMatrix( invMatrix2, tmpS1 );
    ApplyRotationMatrix( invMatrix2, tmpS2 );
    ApplyRotationMatrix( invMatrix2, tmpS1Norm );
    ApplyRotationMatrix( invMatrix2, tmpS2Norm );
    ApplyRotationMatrix( invMatrix2, qCart );
    ApplyRotationMatrix( invMatrix2, pCart );
    
    /* First rotation */
    ApplyRotationMatrix( invMatrix, rHat );
    ApplyRotationMatrix( invMatrix, vHat );
    ApplyRotationMatrix( invMatrix, LnHat );
    ApplyRotationMatrix( invMatrix, tmpS1 );
    ApplyRotationMatrix( invMatrix, tmpS2 );
    ApplyRotationMatrix( invMatrix, tmpS1Norm );
    ApplyRotationMatrix( invMatrix, tmpS2Norm );
    ApplyRotationMatrix( invMatrix, qCart );
    ApplyRotationMatrix( invMatrix, pCart );
    

    /* If required, apply the tortoise transform */
    // if ( tmpTortoise )
    // {
    //     REAL8 r = sqrt(qCart[0]*qCart[0] + qCart[1]*qCart[1] + qCart[2]*qCart[2] );
    //     REAL8 deltaR = XLALSimIMRSpinEOBHamiltonianDeltaR( params->seobCoeffs, r, eta, a );
    //     REAL8 deltaT = XLALSimIMRSpinEOBHamiltonianDeltaT( params->seobCoeffs, r, eta, a );
    //     REAL8 csi    = sqrt( deltaT * deltaR )/(r*r + a*a);
        
    //     REAL8 pr = (qCart[0]*pCart[0] + qCart[1]*pCart[1] + qCart[2]*pCart[2])/r;
    //     params->tortoise = tmpTortoise;
        
    //     //printf( "Applying the tortoise to p (csi = %.26e)\n", csi );
        
    //     for ( i = 0; i < 3; i++ )
    //     {
    //         pCart[i] = pCart[i] + qCart[i] * pr * (csi - 1.) / r;
    //     }
    // }

    /* Now copy the initial conditions back to the return vector */
    memcpy( initConds->data, qCart, sizeof(qCart) );
    memcpy( initConds->data+3, pCart, sizeof(pCart) );
    memcpy( initConds->data+6, tmpS1Norm, sizeof(tmpS1Norm) );
    memcpy( initConds->data+9, tmpS2Norm, sizeof(tmpS2Norm) );
    
    gsl_matrix_free(rotMatrix2);
    gsl_matrix_free(invMatrix2);
    
    gsl_matrix_free(rotMatrix);
    gsl_matrix_free(invMatrix);
    
    //printf( "THE FINAL INITIAL CONDITIONS:\n");
    /*printf( " %.16e %.16e %.16e\n%.16e %.16e %.16e\n%.16e %.16e %.16e\n%.16e %.16e %.16e\n", initConds->data[0], initConds->data[1], initConds->data[2],
     initConds->data[3], initConds->data[4], initConds->data[5], initConds->data[6], initConds->data[7], initConds->data[8],
     initConds->data[9], initConds->data[10], initConds->data[11] );*/
    
    return CEV_SUCCESS;
}


/* For precessing case */
INT EOBInitialConditionsPrec_e_anomaly(REAL8Vector    *initConds,
                         const REAL8    mass1,
                         const REAL8    mass2,
                         const REAL8    fMin,
                         const REAL8    e0,
                         const REAL8    zeta,
                         const REAL8    xi,
                         const REAL8    inc,
                         const REAL8    spin1[],
                         const REAL8    spin2[],
                         SpinEOBParams  *params)
{
	static const int lMax = 8;

	int		i;

	/* Variable to keep track of whether the user requested the tortoise */
	int		tmpTortoise;

	UINT		SpinAlignedEOBversion;

	REAL8		mTotal;
	REAL8		eta;
	REAL8		omega   , v0;	/* Initial velocity and angular
					 * frequency */

	REAL8		ham;	/* Hamiltonian */

	REAL8		LnHat    [3];	/* Initial orientation of angular
					 * momentum */
	REAL8		rHat     [3];	/* Initial orientation of radial
					 * vector */
	REAL8		vHat     [3];	/* Initial orientation of velocity
					 * vector */
	REAL8		Lhat     [3];	/* Direction of relativistic ang mom */
	REAL8		qHat     [3];
	REAL8		pHat     [3];

	/* q and p vectors in Cartesian and spherical coords */
	REAL8		qCart    [3], pCart[3];
	REAL8		qSph     [3], pSph[3];

	/* We will need to manipulate the spin vectors */
	/* We will use temporary vectors to do this */
	REAL8		tmpS1    [3];
	REAL8		tmpS2    [3];
	REAL8		tmpS1Norm[3];
	REAL8		tmpS2Norm[3];

	REAL8Vector	qCartVec, pCartVec;
	REAL8Vector	s1Vec, s2Vec, s1VecNorm, s2VecNorm;
	REAL8Vector	sKerr, sStar;
	REAL8		sKerrData[3], sStarData[3];
	REAL8		a = 0.;
	//, chiS, chiA;
	//REAL8 chi1, chi2;

	/*
	 * We will need a full values vector for calculating derivs of
	 * Hamiltonian
	 */
	REAL8		sphValues[12];
	REAL8		cartValues[12];

	/* Matrices for rotating to the new basis set. */
	/* It is more convenient to calculate the ICs in a simpler basis */
	gsl_matrix     *rotMatrix = NULL;
	gsl_matrix     *invMatrix = NULL;
	gsl_matrix     *rotMatrix2 = NULL;
	gsl_matrix     *invMatrix2 = NULL;

	/* Root finding stuff for finding the spherical orbit */
	SEOBRootParams	rootParams;
	const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrids;
	gsl_multiroot_fsolver *rootSolver = NULL;

	gsl_multiroot_function rootFunction;
	gsl_vector     *initValues = NULL;
	gsl_vector     *finalValues = NULL;
	INT gslStatus;
    INT cntGslNoProgress = 0, MAXcntGslNoProgress = 5;
    //INT cntGslNoProgress = 0, MAXcntGslNoProgress = 50;
    REAL8 multFacGslNoProgress = 3./5.;
	//const int	maxIter = 2000;
	const int	maxIter = 10000;

	memset(&rootParams, 0, sizeof(rootParams));

	mTotal = mass1 + mass2;
	eta = mass1 * mass2 / (mTotal * mTotal);
	memcpy(tmpS1, spin1, sizeof(tmpS1));
	memcpy(tmpS2, spin2, sizeof(tmpS2));
	memcpy(tmpS1Norm, spin1, sizeof(tmpS1Norm));
	memcpy(tmpS2Norm, spin2, sizeof(tmpS2Norm));
	for (i = 0; i < 3; i++) {
		tmpS1Norm[i] /= mTotal * mTotal;
		tmpS2Norm[i] /= mTotal * mTotal;
	}
	SpinAlignedEOBversion = 4;
	/* We compute the ICs for the non-tortoise p, and convert at the end */
	tmpTortoise = params->tortoise;
	params->tortoise = 0;

	EOBNonQCCoeffs *nqcCoeffs = NULL;
	nqcCoeffs = params->nqcCoeffs;

	/*
	 * STEP 1) Rotate to LNhat0 along z-axis and N0 along x-axis frame,
	 * where LNhat0 and N0 are initial normal to orbital plane and
	 * initial orbital separation;
	 */

	/* Set the initial orbital ang mom direction. Taken from STPN code */
	LnHat[0] = sin(inc);
	LnHat[1] = 0.;
	LnHat[2] = cos(inc);

	/*
	 * Set the radial direction - need to take care to avoid singularity
	 * if L is along z axis
	 */
	if (LnHat[2] > 0.9999) {
		rHat[0] = 1.;
		rHat[1] = rHat[2] = 0.;
	} else {
		REAL8		theta0 = atan(-LnHat[2] / LnHat[0]);	/* theta0 is between 0
									 * and Pi */
		rHat[0] = sin(theta0);
		rHat[1] = 0;
		rHat[2] = cos(theta0);
	}

	/* Now we can complete the triad */
	vHat[0] = CalculateCrossProduct(0, LnHat, rHat);
	vHat[1] = CalculateCrossProduct(1, LnHat, rHat);
	vHat[2] = CalculateCrossProduct(2, LnHat, rHat);

	NormalizeVector(vHat);

	/* Vectors BEFORE rotation */
#if 0
		for (i = 0; i < 3; i++)
			print_debug(" LnHat[%d] = %.16e, rHat[%d] = %.16e, vHat[%d] = %.16e\n",
        i, LnHat[i], i, rHat[i], i, vHat[i]);

		print_debug("\n\n");
		for (i = 0; i < 3; i++)
			print_debug(" s1[%d] = %.16e, s2[%d] = %.16e\n", i, tmpS1[i], i, tmpS2[i]);
    //fflush(NULL);
#endif

	/* Allocate and compute the rotation matrices */
	rotMatrix = gsl_matrix_alloc(3, 3);
	invMatrix = gsl_matrix_alloc(3, 3);
	if (!rotMatrix || !invMatrix) 
    {
		if (rotMatrix)
			gsl_matrix_free(rotMatrix);
		if (invMatrix)
			gsl_matrix_free(invMatrix);
		return CEV_FAILURE;
	}
	if (CalculateRotationMatrix(rotMatrix, invMatrix, rHat, vHat, LnHat) == CEV_FAILURE) {
		gsl_matrix_free(rotMatrix);
		gsl_matrix_free(invMatrix);
		return CEV_FAILURE;
	}
	/* Rotate the orbital vectors and spins */
	ApplyRotationMatrix(rotMatrix, rHat);
	ApplyRotationMatrix(rotMatrix, vHat);
	ApplyRotationMatrix(rotMatrix, LnHat);
	ApplyRotationMatrix(rotMatrix, tmpS1);
	ApplyRotationMatrix(rotMatrix, tmpS2);
	ApplyRotationMatrix(rotMatrix, tmpS1Norm);
	ApplyRotationMatrix(rotMatrix, tmpS2Norm);

	/* See if Vectors have been rotated fine */
#if 0
		print_debug("\nAfter applying rotation matrix:\n\n");
		for (i = 0; i < 3; i++)
			print_debug(" LnHat[%d] = %.16e, rHat[%d] = %.16e, vHat[%d] = %.16e\n",
                i, LnHat[i], i, rHat[i], i, vHat[i]);

		print_debug("\n");
		for (i = 0; i < 3; i++)
			print_debug(" s1[%d] = %.16e, s2[%d] = %.16e\n", i, tmpS1[i], i, tmpS2[i]);

    //fflush(NULL);
#endif
	/*
	 * STEP 2) After rotation in STEP 1, in spherical coordinates, phi0
	 * and theta0 are given directly in Eq. (4.7), r0, pr0, ptheta0 and
	 * pphi0 are obtained by solving Eqs. (4.8) and (4.9) (using
	 * gsl_multiroot_fsolver). At this step, we find initial conditions
	 * for a spherical orbit without radiation reaction.
	 */

  /* Initialise the gsl stuff */
	rootSolver = gsl_multiroot_fsolver_alloc(T, 3);
	if (!rootSolver) {
		gsl_matrix_free(rotMatrix);
		gsl_matrix_free(invMatrix);
		return CEV_FAILURE;
	}
	initValues = gsl_vector_calloc(3);
	if (!initValues) {
		gsl_multiroot_fsolver_free(rootSolver);
		gsl_matrix_free(rotMatrix);
		gsl_matrix_free(invMatrix);
		return CEV_FAILURE;
	}

	rootFunction.f = XLALFindSphericalOrbitPrec;
	rootFunction.n = 3;
	rootFunction.params = &rootParams;

    /* Set to use optimized or unoptimized code */
    // rootParams.use_optimized = 0;

	/* Calculate the initial velocity from the given initial frequency */
    REAL8 fMinE = fMin;
    rootParams.e0 = 0.0;
    if (CODE_VERSION == 1)
        rootParams.e0 = e0;
    else if (CODE_VERSION == 2)
    {
        fMinE /= pow(1-e0*e0, 1.5);
    }
	omega = CST_PI * mTotal * CST_MTSUN_SI * fMinE;
	v0 = cbrt(omega);

	/* Given this, we can start to calculate the initial conditions */
	/* for spherical coords in the new basis */
	rootParams.omega = omega;
	rootParams.params = params;
	/* To start with, we will just assign Newtonian-ish ICs to the system */


	rootParams.values[0] = scale1 * sqrt((1.-e0)*(1.-e0) / (v0*v0*v0*v0) - 36.0);	/* Initial r */
	rootParams.values[4] = scale2 * pow(1.-e0, 0.5) * v0;	            /* Initial p */
	rootParams.values[5] = scale3 * 1e-3;
	//PK
    memcpy(rootParams.values + 6, tmpS1, sizeof(tmpS1));
	memcpy(rootParams.values + 9, tmpS2, sizeof(tmpS2));

#if 0
    print_debug("omega = %.16e\n", omega);
    print_debug("ICs guess: x = %.16e, py = %.16e, pz = %.16e\n",
      rootParams.values[0]/scale1, rootParams.values[4]/scale2,
      rootParams.values[5]/scale3);
    //fflush(NULL);
#endif
	gsl_vector_set(initValues, 0, rootParams.values[0]);
	gsl_vector_set(initValues, 1, rootParams.values[4]);
    gsl_vector_set(initValues, 2, rootParams.values[5]);

	gsl_multiroot_fsolver_set(rootSolver, &rootFunction, initValues);

	/* We are now ready to iterate to find the solution */
	i = 0;

    INT4 jittered=0;
	REAL8 r_now = 0.0;
	REAL8 temp = 0.0;
	do 
    {
            gslStatus = gsl_multiroot_fsolver_iterate(rootSolver);
            if (gslStatus == GSL_ENOPROG || gslStatus == GSL_ENOPROGJ) 
            {
                PRINT_LOG_INFO(LOG_ERROR, "NO PROGRESS being made by Spherical orbit root solver\n");

                /* Print Residual Function values whose roots we are trying to find */
                finalValues = gsl_multiroot_fsolver_f(rootSolver);
                PRINT_LOG_INFO(LOG_ERROR, "Function value here given by the following:\n");
                PRINT_LOG_INFO(LOG_ERROR, " F1 = %.16e, F2 = %.16e, F3 = %.16e\n",
                gsl_vector_get(finalValues, 0),
                gsl_vector_get(finalValues, 1), gsl_vector_get(finalValues, 2));

                /* Print Step sizes in each of function variables */
                finalValues = gsl_multiroot_fsolver_dx(rootSolver);
                // XLAL_PRINT_INFO("Stepsizes in each dimension:\n");
                // XLAL_PRINT_INFO(" x = %.16e, py = %.16e, pz = %.16e\n",
                // gsl_vector_get(finalValues, 0)/scale1,
                // gsl_vector_get(finalValues, 1)/scale2,
                // gsl_vector_get(finalValues, 2)/scale3);

                /* Only allow this flag to be caught MAXcntGslNoProgress no. of times */
                cntGslNoProgress += 1;
                if (cntGslNoProgress >= MAXcntGslNoProgress) 
                {
                    cntGslNoProgress = 0;

                    if(multFacGslNoProgress < 1.) { multFacGslNoProgress *= 1.02; }
                    else { multFacGslNoProgress /= 1.01; }

                }
                /* Now that no progress is being made, we need to reset the initial guess
                * for the (r,pPhi, pTheta) and reset the integrator */
                rootParams.values[0] = scale1 * sqrt(1. / (v0 * v0)*1./(v0*v0) -36.0);	/* Initial r */
                rootParams.values[4] = scale2 * v0;	            /* Initial p */
                if( cntGslNoProgress % 2 )
                    rootParams.values[5] = scale3 * 1e-3 / multFacGslNoProgress;
                else
                    rootParams.values[5] = scale3 * 1e-3 * multFacGslNoProgress;
                //PK
                memcpy(rootParams.values + 6, tmpS1, sizeof(tmpS1));
                memcpy(rootParams.values + 9, tmpS2, sizeof(tmpS2));
                // XLAL_PRINT_INFO("New ICs guess: x = %.16e, py = %.16e, pz = %.16e\n",
                //         rootParams.values[0]/scale1, rootParams.values[4]/scale2,
                //         rootParams.values[5]/scale3);
                // fflush(NULL);
                gsl_vector_set(initValues, 0, rootParams.values[0]);
                gsl_vector_set(initValues, 1, rootParams.values[4]);
                gsl_vector_set(initValues, 2, rootParams.values[5]);
                gsl_multiroot_fsolver_set(rootSolver, &rootFunction, initValues);
            }
            else if (gslStatus == GSL_EBADFUNC) 
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "Inf or Nan encountered in evaluluation of spherical orbit Eqn");
                gsl_multiroot_fsolver_free(rootSolver);
                gsl_vector_free(initValues);
                gsl_matrix_free(rotMatrix);
                gsl_matrix_free(invMatrix);
                return CEV_FAILURE;
            }
            else if (gslStatus != GSL_SUCCESS) 
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "Error in GSL iteration function!");
                gsl_multiroot_fsolver_free(rootSolver);
                gsl_vector_free(initValues);
                gsl_matrix_free(rotMatrix);
                gsl_matrix_free(invMatrix);
                return CEV_FAILURE;
            }

        /* different ways to test convergence of the method */
		gslStatus = gsl_multiroot_test_residual(rootSolver->f, 1.0e-9);
        /*XLAL_CALLGSL(gslStatus= gsl_multiroot_test_delta(
          gsl_multiroot_fsolver_dx(rootSolver),
          gsl_multiroot_fsolver_root(rootSolver),
          1.e-8, 1.e-5));*/

		if (jittered==0) 
        {
            finalValues = gsl_multiroot_fsolver_dx(rootSolver);
            if (isnan(gsl_vector_get(finalValues, 1))) 
            {
                rootParams.values[0] = scale1 * sqrt(1. / (v0 * v0)*1/(v0*v0)*(1.+1.e-8) - 36.0);	/* Jitter on initial r */
                rootParams.values[4] = scale2 * v0*(1.-1.e-8);	            /* Jitter on initial p */
                rootParams.values[5] = scale3 * 1e-3;
                memcpy(rootParams.values + 6, tmpS1, sizeof(tmpS1));
                memcpy(rootParams.values + 9, tmpS2, sizeof(tmpS2));
                gsl_vector_set(initValues, 0, rootParams.values[0]);
                gsl_vector_set(initValues, 1, rootParams.values[4]);
                gsl_vector_set(initValues, 2, rootParams.values[5]);
                gsl_multiroot_fsolver_set(rootSolver, &rootFunction, initValues);
                jittered=1;
            }
		}
		i++;
	}
	while (gslStatus == GSL_CONTINUE && i <= maxIter);

    // if(debugPK) { fflush(NULL); fclose(out); }

	if (i > maxIter && gslStatus != GSL_SUCCESS) 
    {
		gsl_multiroot_fsolver_free(rootSolver);
		gsl_vector_free(initValues);
		gsl_matrix_free(rotMatrix);
		gsl_matrix_free(invMatrix);
		//XLAL_ERROR(XLAL_EMAXITER);
		return CEV_FAILURE;
	}
	finalValues = gsl_multiroot_fsolver_root(rootSolver);
#if 0
    print_debug("Spherical orbit conditions here given by the following:\n");
    print_debug(" x = %.16e, py = %.16e, pz = %.16e\n",
    gsl_vector_get(finalValues, 0)/scale1,
    gsl_vector_get(finalValues, 1)/scale2,
    gsl_vector_get(finalValues, 2)/scale3);
#endif
	memset(qCart, 0, sizeof(qCart));
	memset(pCart, 0, sizeof(pCart));

	qCart[0] = sqrt(gsl_vector_get(finalValues, 0)*gsl_vector_get(finalValues, 0)+36.0);
	pCart[1] = gsl_vector_get(finalValues, 1)/scale2;
	pCart[2] = gsl_vector_get(finalValues, 2)/scale3;
    if (CODE_VERSION == 2)
    {
        qCart[0] /= 1. + e0;
        pCart[1] *= 1. + e0;
    }

	/* Free the GSL root finder, since we're done with it */
	gsl_multiroot_fsolver_free(rootSolver);
	gsl_vector_free(initValues);


	/*
	 * STEP 3) Rotate to L0 along z-axis and N0 along x-axis frame, where
	 * L0 is the initial orbital angular momentum and L0 is calculated
	 * using initial position and linear momentum obtained in STEP 2.
	 */

	/* Now we can calculate the relativistic L and rotate to a new basis */
	memcpy(qHat, qCart, sizeof(qCart));
	memcpy(pHat, pCart, sizeof(pCart));

	NormalizeVector(qHat);
	NormalizeVector(pHat);

	Lhat[0] = CalculateCrossProduct(0, qHat, pHat);
	Lhat[1] = CalculateCrossProduct(1, qHat, pHat);
	Lhat[2] = CalculateCrossProduct(2, qHat, pHat);

	NormalizeVector(Lhat);

	rotMatrix2 = gsl_matrix_alloc(3, 3);
	invMatrix2 = gsl_matrix_alloc(3, 3);

	if (CalculateRotationMatrix(rotMatrix2, invMatrix2, qHat, pHat, Lhat) == CEV_FAILURE) 
    {
		gsl_matrix_free(rotMatrix);
		gsl_matrix_free(invMatrix);
		return CEV_FAILURE;
	}
	ApplyRotationMatrix(rotMatrix2, rHat);
	ApplyRotationMatrix(rotMatrix2, vHat);
	ApplyRotationMatrix(rotMatrix2, LnHat);
	ApplyRotationMatrix(rotMatrix2, tmpS1);
	ApplyRotationMatrix(rotMatrix2, tmpS2);
	ApplyRotationMatrix(rotMatrix2, tmpS1Norm);
	ApplyRotationMatrix(rotMatrix2, tmpS2Norm);
	ApplyRotationMatrix(rotMatrix2, qCart);
	ApplyRotationMatrix(rotMatrix2, pCart);

    gsl_matrix_free(rotMatrix);
    gsl_matrix_free(rotMatrix2);

    // XLAL_PRINT_INFO("qCart after rotation2 %3.10f %3.10f %3.10f\n", qCart[0], qCart[1], qCart[2]);
    // XLAL_PRINT_INFO("pCart after rotation2 %3.10f %3.10f %3.10f\n", pCart[0], pCart[1], pCart[2]);
    // XLAL_PRINT_INFO("S1 after rotation2 %3.10f %3.10f %3.10f\n", tmpS1Norm[0], tmpS1Norm[1], tmpS1Norm[2]);
    // XLAL_PRINT_INFO("S2 after rotation2 %3.10f %3.10f %3.10f\n", tmpS2Norm[0], tmpS2Norm[1], tmpS2Norm[2]);
	/*
	 * STEP 4) In the L0-N0 frame, we calculate (dE/dr)|sph using Eq.
	 * (4.14), then initial dr/dt using Eq. (4.10), and finally pr0 using
	 * Eq. (4.15).
	 */

	/* Now we can calculate the flux. Change to spherical co-ords */
	CartesianToSpherical(qSph, pSph, qCart, pCart);
	memcpy(sphValues, qSph, sizeof(qSph));
	memcpy(sphValues + 3, pSph, sizeof(pSph));
	memcpy(sphValues + 6, tmpS1, sizeof(tmpS1));
	memcpy(sphValues + 9, tmpS2, sizeof(tmpS2));

	memcpy(cartValues, qCart, sizeof(qCart));
	memcpy(cartValues + 3, pCart, sizeof(pCart));
	memcpy(cartValues + 6, tmpS1, sizeof(tmpS1));
	memcpy(cartValues + 9, tmpS2, sizeof(tmpS2));

	REAL8		dHdpphi , d2Hdr2, d2Hdrdpphi;
	REAL8		rDot    , dHdpr, flux, dEdr;

	d2Hdr2 = XLALCalculateSphHamiltonianDeriv2Prec(0, 0, sphValues, params);
	d2Hdrdpphi = XLALCalculateSphHamiltonianDeriv2Prec(0, 5, sphValues, params);

    // XLAL_PRINT_INFO("d2Hdr2 = %.16e, d2Hdrdpphi = %.16e\n", d2Hdr2, d2Hdrdpphi);

	/* New code to compute derivatives w.r.t. cartesian variables */

	REAL8		tmpDValues[14];
	int status;
	for (i = 0; i < 3; i++) {
		cartValues[i + 6] /= mTotal * mTotal;
		cartValues[i + 9] /= mTotal * mTotal;
	}
    UINT oldignoreflux = params->ignoreflux;
    params->ignoreflux = 1;
    status = XLALSpinPrecHcapNumericalDerivative(0, cartValues, tmpDValues, params);
    params->ignoreflux = oldignoreflux;
	for (i = 0; i < 3; i++) {
		cartValues[i + 6] *= mTotal * mTotal;
		cartValues[i + 9] *= mTotal * mTotal;
	}

	dHdpphi = tmpDValues[1] / sqrt(cartValues[0] * cartValues[0] + cartValues[1] * cartValues[1] + cartValues[2] * cartValues[2]);
	//XLALSpinPrecHcapNumDerivWRTParam(4, cartValues, params) / sphValues[0];

	dEdr = -dHdpphi * d2Hdr2 / d2Hdrdpphi;

    // XLAL_PRINT_INFO("d2Hdr2 = %.16e d2Hdrdpphi = %.16e dHdpphi = %.16e\n",
    //     d2Hdr2, d2Hdrdpphi, dHdpphi);

	if (d2Hdr2 != 0.0 && e0 == 0.0) 
    {
		/* We will need to calculate the Hamiltonian to get the flux */
		s1Vec.length = s2Vec.length = s1VecNorm.length = s2VecNorm.length = sKerr.length = sStar.length = 3;
		s1Vec.data = tmpS1;
		s2Vec.data = tmpS2;
		s1VecNorm.data = tmpS1Norm;
		s2VecNorm.data = tmpS2Norm;
		sKerr.data = sKerrData;
		sStar.data = sStarData;

		qCartVec.length = pCartVec.length = 3;
		qCartVec.data = qCart;
		pCartVec.data = pCart;

		//chi1 = tmpS1[0] * LnHat[0] + tmpS1[1] * LnHat[1] + tmpS1[2] * LnHat[2];
		//chi2 = tmpS2[0] * LnHat[0] + tmpS2[1] * LnHat[1] + tmpS2[2] * LnHat[2];

		//if (debugPK)
			//XLAL_PRINT_INFO("magS1 = %.16e, magS2 = %.16e\n", chi1, chi2);

		//chiS = 0.5 * (chi1 / (mass1 * mass1) + chi2 / (mass2 * mass2));
		//chiA = 0.5 * (chi1 / (mass1 * mass1) - chi2 / (mass2 * mass2));

		EOBCalculateSigmaKerr(&sKerr, &s1VecNorm, &s2VecNorm);
		EOBCalculateSigmaStar(&sStar, mass1, mass2, &s1VecNorm, &s2VecNorm);

		/*
		 * The a in the flux has been set to zero, but not in the
		 * Hamiltonian
		 */
		a = sqrt(sKerr.data[0] * sKerr.data[0] + sKerr.data[1] * sKerr.data[1] + sKerr.data[2] * sKerr.data[2]);
		//XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients(params->eobParams->hCoeffs, mass1, mass2, eta, /* a */ 0.0, chiS, chiA);
		//XLALSimIMRCalculateSpinPrecEOBHCoeffs(params->seobCoeffs, eta, a);
		ham = XLALSimIMRSpinPrecEOBHamiltonian(eta, &qCartVec, &pCartVec, &s1VecNorm, &s2VecNorm, &sKerr, &sStar, params->tortoise, params->seobCoeffs, params->hParams);

        // XLAL_PRINT_INFO("Stas: hamiltonian in ICs at this point is %.16e\n", ham);

		/* And now, finally, the flux */
		REAL8Vector	polarDynamics, cartDynamics;
		REAL8		polarData[4], cartData[12];

		polarDynamics.length = 4;
		polarDynamics.data = polarData;

		polarData[0] = qSph[0];
		polarData[1] = 0.;
		polarData[2] = pSph[0];
		polarData[3] = pSph[2];

		cartDynamics.length = 12;
		cartDynamics.data = cartData;

		memcpy(cartData, qCart, 3 * sizeof(REAL8));
		memcpy(cartData + 3, pCart, 3 * sizeof(REAL8));
		memcpy(cartData + 6, tmpS1Norm, 3 * sizeof(REAL8));
		memcpy(cartData + 9, tmpS2Norm, 3 * sizeof(REAL8));

		//XLAL_PRINT_INFO("Stas: starting FLux calculations\n");
        flux = XLALInspiralPrecSpinFactorizedFlux(&polarDynamics, &cartDynamics, nqcCoeffs, omega, 0, qSph[0]*omega, params, ham, lMax, SpinAlignedEOBversion);
		/*
		 * flux  = XLALInspiralSpinFactorizedFlux( &polarDynamics,
		 * nqcCoeffs, omega, params, ham, lMax, SpinAlignedEOBversion
		 * );
		 */
		PRINT_LOG_INFO(LOG_DEBUG, "Initial Conditions: Stas flux = %.16e", flux);
		//exit(0);
		flux = flux / eta;
        if (e0 != 0.0)
            rDot = 0.0;
        else
		    rDot = -flux / dEdr;
		/*
		 * We now need dHdpr - we take it that it is safely linear up
		 * to a pr of 1.0e-3 PK: Ideally, the pr should be of the
		 * order of other momenta, in order for its contribution to
		 * the Hamiltonian to not get buried in the numerical noise
		 * in the numerically larger momenta components
		 */
		cartValues[3] = 1.0e-3;
		for (i = 0; i < 3; i++) 
        {
			cartValues[i + 6] /= mTotal * mTotal;
			cartValues[i + 9] /= mTotal * mTotal;
		}
        oldignoreflux = params->ignoreflux;
        params->ignoreflux = 1;
        params->seobCoeffs->updateHCoeffs = 1;
	    status = XLALSpinPrecHcapNumericalDerivative(0, cartValues, tmpDValues, params);
        params->ignoreflux = oldignoreflux;
		for (i = 0; i < 3; i++) 
        {
			cartValues[i + 6] *= mTotal * mTotal;
			cartValues[i + 9] *= mTotal * mTotal;
		}
        REAL8   csi = sqrt(XLALSimIMRSpinPrecEOBHamiltonianDeltaT(params->seobCoeffs, qSph[0], eta, a)*XLALSimIMRSpinPrecEOBHamiltonianDeltaR(params->seobCoeffs, qSph[0], eta, a)) / (qSph[0] * qSph[0] + a * a);

		dHdpr = csi*tmpDValues[0];
		//XLALSpinPrecHcapNumDerivWRTParam(3, cartValues, params);

        // XLAL_PRINT_INFO("Ingredients going into prDot:\n");
        // XLAL_PRINT_INFO("flux = %.16e, dEdr = %.16e, dHdpr = %.16e, dHdpr/pr = %.16e\n", flux, dEdr, dHdpr, dHdpr / cartValues[3]);
		/*
		 * We can now calculate what pr should be taking into account
		 * the flux
		 */
		pSph[0] = rDot / (dHdpr / cartValues[3]);
	} else 
    {
		/*
		 * Since d2Hdr2 has evaluated to zero, we cannot do the
		 * above. Just set pr to zero
		 */
		//XLAL_PRINT_INFO("d2Hdr2 is zero!\n");
		pSph[0] = 0;
	}

	/* Now we are done - convert back to cartesian coordinates ) */
	SphericalToCartesian(qCart, pCart, qSph, pSph);
    PRINT_LOG_INFO(LOG_DEBUG, "Sph initial condition : r = (%e,%e,%e), p = (%e,%e,%e)", qSph[0], qSph[1], qSph[2], pSph[0], pSph[1], pSph[2]);

	/*
	 * STEP 5) Rotate back to the original inertial frame by inverting
	 * the rotation of STEP 3 and then  inverting the rotation of STEP 1.
	 */

	/* Undo rotations to get back to the original basis */
	/* Second rotation */
	ApplyRotationMatrix(invMatrix2, rHat);
	ApplyRotationMatrix(invMatrix2, vHat);
	ApplyRotationMatrix(invMatrix2, LnHat);
	ApplyRotationMatrix(invMatrix2, tmpS1);
	ApplyRotationMatrix(invMatrix2, tmpS2);
	ApplyRotationMatrix(invMatrix2, tmpS1Norm);
	ApplyRotationMatrix(invMatrix2, tmpS2Norm);
	ApplyRotationMatrix(invMatrix2, qCart);
	ApplyRotationMatrix(invMatrix2, pCart);

	/* First rotation */
	ApplyRotationMatrix(invMatrix, rHat);
	ApplyRotationMatrix(invMatrix, vHat);
	ApplyRotationMatrix(invMatrix, LnHat);
	ApplyRotationMatrix(invMatrix, tmpS1);
	ApplyRotationMatrix(invMatrix, tmpS2);
	ApplyRotationMatrix(invMatrix, tmpS1Norm);
	ApplyRotationMatrix(invMatrix, tmpS2Norm);
	ApplyRotationMatrix(invMatrix, qCart);
	ApplyRotationMatrix(invMatrix, pCart);

    gsl_matrix_free(invMatrix);
    gsl_matrix_free(invMatrix2);

        /* If required, apply the tortoise transform */
	if (tmpTortoise) 
    {
		REAL8		r = sqrt(qCart[0] * qCart[0] + qCart[1] * qCart[1] + qCart[2] * qCart[2]);
		REAL8		deltaR = XLALSimIMRSpinPrecEOBHamiltonianDeltaR(params->seobCoeffs, r, eta, a);
		REAL8		deltaT = XLALSimIMRSpinPrecEOBHamiltonianDeltaT(params->seobCoeffs, r, eta, a);
		REAL8		csi = sqrt(deltaT * deltaR) / (r * r + a * a);

		REAL8		pr = (qCart[0] * pCart[0] + qCart[1] * pCart[1] + qCart[2] * pCart[2]) / r;

		params->tortoise = tmpTortoise;

        // XLAL_PRINT_INFO("Applying the tortoise to p (csi = %.26e)\n", csi);
        // XLAL_PRINT_INFO("pCart = %3.10f %3.10f %3.10f\n", pCart[0], pCart[1], pCart[2]);
		for (i = 0; i < 3; i++) 
        {
			pCart[i] = pCart[i] + qCart[i] * pr * (csi - 1.) / r;
		}
	}


    /* Now copy the initial conditions back to the return vector */
	memcpy(initConds->data, qCart, sizeof(qCart));
	memcpy(initConds->data + 3, pCart, sizeof(pCart));
	memcpy(initConds->data + 6, tmpS1Norm, sizeof(tmpS1Norm));
	memcpy(initConds->data + 9, tmpS2Norm, sizeof(tmpS2Norm));
    for (i=0; i<12; i++) 
    {
        if (fabs(initConds->data[i]) <=1.0e-15) 
        {
            initConds->data[i] = 0.;
        }
    }
    // REAL8 initr = sqrt(qCart[0]*qCart[0]+qCart[1]*qCart[1]+qCart[2]*qCart[2]);
    // if ( initr< 3.)
    // {
    //     PRINT_LOG_INFO(LOG_CRITICAL, "the initial polar r = %.8e is too small, abort.", initr);
    //     return CEV_FAILURE;
    // }

    // XLAL_PRINT_INFO("THE FINAL INITIAL CONDITIONS:\n");
    // XLAL_PRINT_INFO(" %.16e %.16e %.16e\n%.16e %.16e %.16e\n%.16e %.16e %.16e\n%.16e %.16e %.16e\n", initConds->data[0], initConds->data[1], initConds->data[2],
    //         initConds->data[3], initConds->data[4], initConds->data[5], initConds->data[6], initConds->data[7], initConds->data[8],
    //         initConds->data[9], initConds->data[10], initConds->data[11]);
    return CEV_SUCCESS;
}