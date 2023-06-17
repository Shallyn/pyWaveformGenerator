/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "mySpherical.h"
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_gamma.h>

INT SpinWeightedSphericalHarmonic(REAL8 theta,
                                  REAL8 phi,
                                  INT s,
                                  INT l,
                                  INT m,
                                  COMPLEX16 *ret)
{
    REAL8 fac;
    COMPLEX16 ans;
    
    /* sanity checks ... */
    if ( l < abs(s) )
    {
        print_err("Error -%s : Invalid mode s=%d, l=%d, m=%d - require |s| <= l\n", __func__, s, l, m);
        return CEV_FAILURE;
    }
    if ( l < abs(m) )
    {
        print_err("Error -%s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n" ,__func__, s, l, m );
        return CEV_FAILURE;
    }
    
    if ( s == -2 )
    {
        if ( l == 2 )
        {
            switch ( m )
            {
                case -2:
                    fac = GET_SQRT( 5.0 / ( 64.0 * CST_PI ) ) * ( 1.0 - GET_COS( theta ))*( 1.0 - GET_COS( theta ));
                    break;
                case -1:
                    fac = GET_SQRT( 5.0 / ( 16.0 * CST_PI ) ) * GET_SIN( theta )*( 1.0 - GET_COS( theta ));
                    break;
                    
                case 0:
                    fac = GET_SQRT( 15.0 / ( 32.0 * CST_PI ) ) * GET_SIN( theta )*GET_SIN( theta );
                    break;
                    
                case 1:
                    fac = GET_SQRT( 5.0 / ( 16.0 * CST_PI ) ) * GET_SIN( theta )*( 1.0 + GET_COS( theta ));
                    break;
                    
                case 2:
                    fac = GET_SQRT( 5.0 / ( 64.0 * CST_PI ) ) * ( 1.0 + GET_COS( theta ))*( 1.0 + GET_COS( theta ));
                    break;
                default:
                    print_err("Error -%s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
                    return CEV_FAILURE;
                    break;
            } /*  switch (m) */
        }  /* l==2*/
        else if ( l == 3 )
        {
            switch ( m )
            {
                case -3:
                    fac = GET_SQRT(21.0/(2.0*CST_PI))*GET_COS(theta/2.0)*GET_POW(GET_SIN(theta/2.0),5.0);
                    break;
                case -2:
                    fac = GET_SQRT(7.0/(4.0*CST_PI))*(2.0 + 3.0*GET_COS(theta))*GET_POW(GET_SIN(theta/2.0),4.0);
                    break;
                case -1:
                    fac = GET_SQRT(35.0/(2.0*CST_PI))*(GET_SIN(theta) + 4.0*GET_SIN(2.0*theta) - 3.0*GET_SIN(3.0*theta))/32.0;
                    break;
                case 0:
                    fac = (GET_SQRT(105.0/(2.0*CST_PI))*GET_COS(theta)*GET_POW(GET_SIN(theta),2.0))/4.0;
                    break;
                case 1:
                    fac = -GET_SQRT(35.0/(2.0*CST_PI))*(GET_SIN(theta) - 4.0*GET_SIN(2.0*theta) - 3.0*GET_SIN(3.0*theta))/32.0;
                    break;
                    
                case 2:
                    fac = GET_SQRT(7.0/CST_PI)*GET_POW(GET_COS(theta/2.0),4.0)*(-2.0 + 3.0*GET_COS(theta))/2.0;
                    break;
                    
                case 3:
                    fac = -GET_SQRT(21.0/(2.0*CST_PI))*GET_POW(GET_COS(theta/2.0),5.0)*GET_SIN(theta/2.0);
                    break;
                    
                default:
                    print_err("Error -%s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
                    return CEV_FAILURE;
                    break;
            }
        }   /* l==3 */
        else if ( l == 4 )
        {
            switch ( m )
            {
                case -4:
                    fac = 3.0*GET_SQRT(7.0/CST_PI)*GET_POW(GET_COS(theta/2.0),2.0)*GET_POW(GET_SIN(theta/2.0),6.0);
                    break;
                case -3:
                    fac = 3.0*GET_SQRT(7.0/(2.0*CST_PI))*GET_COS(theta/2.0)*(1.0 + 2.0*GET_COS(theta))*GET_POW(GET_SIN(theta/2.0),5.0);
                    break;
                    
                case -2:
                    fac = (3.0*(9.0 + 14.0*GET_COS(theta) + 7.0*GET_COS(2.0*theta))*GET_POW(GET_SIN(theta/2.0),4.0))/(4.0*GET_SQRT(CST_PI));
                    break;
                case -1:
                    fac = (3.0*(3.0*GET_SIN(theta) + 2.0*GET_SIN(2.0*theta) + 7.0*GET_SIN(3.0*theta) - 7.0*GET_SIN(4.0*theta)))/(32.0*GET_SQRT(2.0*CST_PI));
                    break;
                case 0:
                    fac = (3.0*GET_SQRT(5.0/(2.0*CST_PI))*(5.0 + 7.0*GET_COS(2.0*theta))*GET_POW(GET_SIN(theta),2.0))/16.0;
                    break;
                case 1:
                    fac = (3.0*(3.0*GET_SIN(theta) - 2.0*GET_SIN(2.0*theta) + 7.0*GET_SIN(3.0*theta) + 7.0*GET_SIN(4.0*theta)))/(32.0*GET_SQRT(2.0*CST_PI));
                    break;
                case 2:
                    fac = (3.0*GET_POW(GET_COS(theta/2.0),4.0)*(9.0 - 14.0*GET_COS(theta) + 7.0*GET_COS(2.0*theta)))/(4.0*GET_SQRT(CST_PI));
                    break;
                case 3:
                    fac = -3.0*GET_SQRT(7.0/(2.0*CST_PI))*GET_POW(GET_COS(theta/2.0),5.0)*(-1.0 + 2.0*GET_COS(theta))*GET_SIN(theta/2.0);
                    break;
                case 4:
                    fac = 3.0*GET_SQRT(7.0/CST_PI)*GET_POW(GET_COS(theta/2.0),6.0)*GET_POW(GET_SIN(theta/2.0),2.0);
                    break;
                default:
                    print_err("Error -%s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
                    return CEV_FAILURE;
                    break;
            }
        }    /* l==4 */
        else if ( l == 5 )
        {
            switch ( m )
            {
                case -5:
                    fac = GET_SQRT(330.0/CST_PI)*GET_POW(GET_COS(theta/2.0),3.0)*GET_POW(GET_SIN(theta/2.0),7.0);
                    break;
                case -4:
                    fac = GET_SQRT(33.0/CST_PI)*GET_POW(GET_COS(theta/2.0),2.0)*(2.0 + 5.0*GET_COS(theta))*GET_POW(GET_SIN(theta/2.0),6.0);
                    break;
                case -3:
                    fac = (GET_SQRT(33.0/(2.0*CST_PI))*GET_COS(theta/2.0)*(17.0 + 24.0*GET_COS(theta) + 15.0*GET_COS(2.0*theta))*GET_POW(GET_SIN(theta/2.0),5.0))/4.0;
                    break;
                case -2:
                    fac = (GET_SQRT(11.0/CST_PI)*(32.0 + 57.0*GET_COS(theta) + 36.0*GET_COS(2.0*theta) + 15.0*GET_COS(3.0*theta))*GET_POW(GET_SIN(theta/2.0),4.0))/8.0;
                    break;
                case -1:
                    fac = (GET_SQRT(77.0/CST_PI)*(2.0*GET_SIN(theta) + 8.0*GET_SIN(2.0*theta) + 3.0*GET_SIN(3.0*theta) + 12.0*GET_SIN(4.0*theta) - 15.0*GET_SIN(5.0*theta)))/256.0;
                    break;
                case 0:
                    fac = (GET_SQRT(1155.0/(2.0*CST_PI))*(5.0*GET_COS(theta) + 3.0*GET_COS(3.0*theta))*GET_POW(GET_SIN(theta),2.0))/32.0;
                    break;
                case 1:
                    fac = GET_SQRT(77.0/CST_PI)*(-2.0*GET_SIN(theta) + 8.0*GET_SIN(2.0*theta) - 3.0*GET_SIN(3.0*theta) + 12.0*GET_SIN(4.0*theta) + 15.0*GET_SIN(5.0*theta))/256.0;
                    break;
                case 2:
                    fac = GET_SQRT(11.0/CST_PI)*GET_POW(GET_COS(theta/2.0),4.0)*(-32.0 + 57.0*GET_COS(theta) - 36.0*GET_COS(2.0*theta) + 15.0*GET_COS(3.0*theta))/8.0;
                    break;
                case 3:
                    fac = -GET_SQRT(33.0/(2.0*CST_PI))*GET_POW(GET_COS(theta/2.0),5.0)*(17.0 - 24.0*GET_COS(theta) + 15.0*GET_COS(2.0*theta))*GET_SIN(theta/2.0)/4.0;
                    break;
                case 4:
                    fac = GET_SQRT(33.0/CST_PI)*GET_POW(GET_COS(theta/2.0),6.0)*(-2.0 + 5.0*GET_COS(theta))*GET_POW(GET_SIN(theta/2.0),2.0);
                    break;
                case 5:
                    fac = -GET_SQRT(330.0/CST_PI)*GET_POW(GET_COS(theta/2.0),7.0)*GET_POW(GET_SIN(theta/2.0),3.0);
                    break;
                default:
                    print_err("Error -%s: nvalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
                    return CEV_FAILURE;
                    break;
            }
        }  /* l==5 */
        else if ( l == 6 )
        {
            switch ( m )
            {
                case -6:
                    fac = (3.*GET_SQRT(715./CST_PI)*GET_POW(GET_COS(theta/2.0),4)*GET_POW(GET_SIN(theta/2.0),8))/2.0;
                    break;
                case -5:
                    fac = (GET_SQRT(2145./CST_PI)*GET_POW(GET_COS(theta/2.0),3)*(1. + 3.*GET_COS(theta))*GET_POW(GET_SIN(theta/2.0),7))/2.0;
                    break;
                case -4:
                    fac = (GET_SQRT(195./(2.0*CST_PI))*GET_POW(GET_COS(theta/2.0),2)*(35. + 44.*GET_COS(theta)
                                                                                      + 33.*GET_COS(2.*theta))*GET_POW(GET_SIN(theta/2.0),6))/8.0;
                    break;
                case -3:
                    fac = (3.*GET_SQRT(13./CST_PI)*GET_COS(theta/2.0)*(98. + 185.*GET_COS(theta) + 110.*GET_COS(2*theta)
                                                                       + 55.*GET_COS(3.*theta))*GET_POW(GET_SIN(theta/2.0),5))/32.0;
                    break;
                case -2:
                    fac = (GET_SQRT(13./CST_PI)*(1709. + 3096.*GET_COS(theta) + 2340.*GET_COS(2.*theta) + 1320.*GET_COS(3.*theta)
                                                 + 495.*GET_COS(4.*theta))*GET_POW(GET_SIN(theta/2.0),4))/256.0;
                    break;
                case -1:
                    fac = (GET_SQRT(65./(2.0*CST_PI))*GET_COS(theta/2.0)*(161. + 252.*GET_COS(theta) + 252.*GET_COS(2.*theta)
                                                                          + 132.*GET_COS(3.*theta) + 99.*GET_COS(4.*theta))*GET_POW(GET_SIN(theta/2.0),3))/64.0;
                    break;
                case 0:
                    fac = (GET_SQRT(1365./CST_PI)*(35. + 60.*GET_COS(2.*theta) + 33.*GET_COS(4.*theta))*GET_POW(GET_SIN(theta),2))/512.0;
                    break;
                case 1:
                    fac = (GET_SQRT(65./(2.0*CST_PI))*GET_POW(GET_COS(theta/2.0),3)*(161. - 252.*GET_COS(theta) + 252.*GET_COS(2.*theta)
                                                                                     - 132.*GET_COS(3.*theta) + 99.*GET_COS(4.*theta))*GET_SIN(theta/2.0))/64.0;
                    break;
                case 2:
                    fac = (GET_SQRT(13./CST_PI)*GET_POW(GET_COS(theta/2.0),4)*(1709. - 3096.*GET_COS(theta) + 2340.*GET_COS(2.*theta)
                                                                               - 1320*GET_COS(3*theta) + 495*GET_COS(4*theta)))/256.0;
                    break;
                case 3:
                    fac = (-3.*GET_SQRT(13./CST_PI)*GET_POW(GET_COS(theta/2.0),5)*(-98. + 185.*GET_COS(theta) - 110.*GET_COS(2*theta)
                                                                                   + 55.*GET_COS(3.*theta))*GET_SIN(theta/2.0))/32.0;
                    break;
                case 4:
                    fac = (GET_SQRT(195./(2.0*CST_PI))*GET_POW(GET_COS(theta/2.0),6)*(35. - 44.*GET_COS(theta)
                                                                                      + 33.*GET_COS(2*theta))*GET_POW(GET_SIN(theta/2.0),2))/8.0;
                    break;
                case 5:
                    fac = -(GET_SQRT(2145./CST_PI)*GET_POW(GET_COS(theta/2.0),7)*(-1. + 3.*GET_COS(theta))*GET_POW(GET_SIN(theta/2.0),3))/2.0;
                    break;
                case 6:
                    fac = (3.*GET_SQRT(715./CST_PI)*GET_POW(GET_COS(theta/2.0),8)*GET_POW(GET_SIN(theta/2.0),4))/2.0;
                    break;
                default:
                    print_err("Error -%s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
                    return CEV_FAILURE;
                    break;
            }
        } /* l==6 */
        else if ( l == 7 )
        {
            switch ( m )
            {
                case -7:
                    fac = GET_SQRT(15015./(2.0*CST_PI))*GET_POW(GET_COS(theta/2.0),5)*GET_POW(GET_SIN(theta/2.0),9);
                    break;
                case -6:
                    fac = (GET_SQRT(2145./CST_PI)*GET_POW(GET_COS(theta/2.0),4)*(2. + 7.*GET_COS(theta))*GET_POW(GET_SIN(theta/2.0),8))/2.0;
                    break;
                case -5:
                    fac = (GET_SQRT(165./(2.0*CST_PI))*GET_POW(GET_COS(theta/2.0),3)*(93. + 104.*GET_COS(theta)
                                                                                      + 91.*GET_COS(2.*theta))*GET_POW(GET_SIN(theta/2.0),7))/8.0;
                    break;
                case -4:
                    fac = (GET_SQRT(165./(2.0*CST_PI))*GET_POW(GET_COS(theta/2.0),2)*(140. + 285.*GET_COS(theta)
                                                                                      + 156.*GET_COS(2.*theta) + 91.*GET_COS(3.*theta))*GET_POW(GET_SIN(theta/2.0),6))/16.0;
                    break;
                case -3:
                    fac = (GET_SQRT(15./(2.0*CST_PI))*GET_COS(theta/2.0)*(3115. + 5456.*GET_COS(theta) + 4268.*GET_COS(2.*theta)
                                                                          + 2288.*GET_COS(3.*theta) + 1001.*GET_COS(4.*theta))*GET_POW(GET_SIN(theta/2.0),5))/128.0;
                    break;
                case -2:
                    fac = (GET_SQRT(15./CST_PI)*(5220. + 9810.*GET_COS(theta) + 7920.*GET_COS(2.*theta) + 5445.*GET_COS(3.*theta)
                                                 + 2860.*GET_COS(4.*theta) + 1001.*GET_COS(5.*theta))*GET_POW(GET_SIN(theta/2.0),4))/512.0;
                    break;
                case -1:
                    fac = (3.*GET_SQRT(5./(2.0*CST_PI))*GET_COS(theta/2.0)*(1890. + 4130.*GET_COS(theta) + 3080.*GET_COS(2.*theta)
                                                                            + 2805.*GET_COS(3.*theta) + 1430.*GET_COS(4.*theta) + 1001.*GET_COS(5*theta))*GET_POW(GET_SIN(theta/2.0),3))/512.0;
                    break;
                case 0:
                    fac = (3.*GET_SQRT(35./CST_PI)*GET_COS(theta)*(109. + 132.*GET_COS(2.*theta)
                                                                   + 143.*GET_COS(4.*theta))*GET_POW(GET_SIN(theta),2))/512.0;
                    break;
                case 1:
                    fac = (3.*GET_SQRT(5./(2.0*CST_PI))*GET_POW(GET_COS(theta/2.0),3)*(-1890. + 4130.*GET_COS(theta) - 3080.*GET_COS(2.*theta)
                                                                                       + 2805.*GET_COS(3.*theta) - 1430.*GET_COS(4.*theta) + 1001.*GET_COS(5.*theta))*GET_SIN(theta/2.0))/512.0;
                    break;
                case 2:
                    fac = (GET_SQRT(15./CST_PI)*GET_POW(GET_COS(theta/2.0),4)*(-5220. + 9810.*GET_COS(theta) - 7920.*GET_COS(2.*theta)
                                                                               + 5445.*GET_COS(3.*theta) - 2860.*GET_COS(4.*theta) + 1001.*GET_COS(5.*theta)))/512.0;
                    break;
                case 3:
                    fac = -(GET_SQRT(15./(2.0*CST_PI))*GET_POW(GET_COS(theta/2.0),5)*(3115. - 5456.*GET_COS(theta) + 4268.*GET_COS(2.*theta)
                                                                                      - 2288.*GET_COS(3.*theta) + 1001.*GET_COS(4.*theta))*GET_SIN(theta/2.0))/128.0;
                    break;
                case 4:
                    fac = (GET_SQRT(165./(2.0*CST_PI))*GET_POW(GET_COS(theta/2.0),6)*(-140. + 285.*GET_COS(theta) - 156.*GET_COS(2*theta)
                                                                                      + 91.*GET_COS(3.*theta))*GET_POW(GET_SIN(theta/2.0),2))/16.0;
                    break;
                case 5:
                    fac = -(GET_SQRT(165./(2.0*CST_PI))*GET_POW(GET_COS(theta/2.0),7)*(93. - 104.*GET_COS(theta)
                                                                                       + 91.*GET_COS(2.*theta))*GET_POW(GET_SIN(theta/2.0),3))/8.0;
                    break;
                case 6:
                    fac = (GET_SQRT(2145./CST_PI)*GET_POW(GET_COS(theta/2.0),8)*(-2. + 7.*GET_COS(theta))*GET_POW(GET_SIN(theta/2.0),4))/2.0;
                    break;
                case 7:
                    fac = -(GET_SQRT(15015./(2.0*CST_PI))*GET_POW(GET_COS(theta/2.0),9)*GET_POW(GET_SIN(theta/2.0),5));
                    break;
                default:
                    print_err("Error : Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
                    return CEV_FAILURE;
                    break;
            }
        } /* l==7 */
        else if ( l == 8 )
        {
            switch ( m )
            {
                case -8:
                    fac = GET_SQRT(34034./CST_PI)*GET_POW(GET_COS(theta/2.0),6)*GET_POW(GET_SIN(theta/2.0),10);
                    break;
                case -7:
                    fac = GET_SQRT(17017./(2.0*CST_PI))*GET_POW(GET_COS(theta/2.0),5)*(1. + 4.*GET_COS(theta))*GET_POW(GET_SIN(theta/2.0),9);
                    break;
                case -6:
                    fac = GET_SQRT(255255./CST_PI)*GET_POW(GET_COS(theta/2.0),4)*(1. + 2.*GET_COS(theta))
                    *GET_SIN(CST_PI/4.0 - theta/2.0)*GET_SIN(CST_PI/4.0 + theta/2.0)*GET_POW(GET_SIN(theta/2.0),8);
                    break;
                case -5:
                    fac = (GET_SQRT(12155./(2.0*CST_PI))*GET_POW(GET_COS(theta/2.0),3)*(19. + 42.*GET_COS(theta)
                                                                                        + 21.*GET_COS(2.*theta) + 14.*GET_COS(3.*theta))*GET_POW(GET_SIN(theta/2.0),7))/8.0;
                    break;
                case -4:
                    fac = (GET_SQRT(935./(2.0*CST_PI))*GET_POW(GET_COS(theta/2.0),2)*(265. + 442.*GET_COS(theta) + 364.*GET_COS(2.*theta)
                                                                                      + 182.*GET_COS(3.*theta) + 91.*GET_COS(4.*theta))*GET_POW(GET_SIN(theta/2.0),6))/32.0;
                    break;
                case -3:
                    fac = (GET_SQRT(561./(2.0*CST_PI))*GET_COS(theta/2.0)*(869. + 1660.*GET_COS(theta) + 1300.*GET_COS(2.*theta)
                                                                           + 910.*GET_COS(3.*theta) + 455.*GET_COS(4.*theta) + 182.*GET_COS(5.*theta))*GET_POW(GET_SIN(theta/2.0),5))/128.0;
                    break;
                case -2:
                    fac = (GET_SQRT(17./CST_PI)*(7626. + 14454.*GET_COS(theta) + 12375.*GET_COS(2.*theta) + 9295.*GET_COS(3.*theta)
                                                 + 6006.*GET_COS(4.*theta) + 3003.*GET_COS(5.*theta) + 1001.*GET_COS(6.*theta))*GET_POW(GET_SIN(theta/2.0),4))/512.0;
                    break;
                case -1:
                    fac = (GET_SQRT(595./(2.0*CST_PI))*GET_COS(theta/2.0)*(798. + 1386.*GET_COS(theta) + 1386.*GET_COS(2.*theta)
                                                                           + 1001.*GET_COS(3.*theta) + 858.*GET_COS(4.*theta) + 429.*GET_COS(5.*theta) + 286.*GET_COS(6.*theta))*GET_POW(GET_SIN(theta/2.0),3))/512.0;
                    break;
                case 0:
                    fac = (3.*GET_SQRT(595./CST_PI)*(210. + 385.*GET_COS(2.*theta) + 286.*GET_COS(4.*theta)
                                                     + 143.*GET_COS(6.*theta))*GET_POW(GET_SIN(theta),2))/4096.0;
                    break;
                case 1:
                    fac = (GET_SQRT(595./(2.0*CST_PI))*GET_POW(GET_COS(theta/2.0),3)*(798. - 1386.*GET_COS(theta) + 1386.*GET_COS(2.*theta)
                                                                                      - 1001.*GET_COS(3.*theta) + 858.*GET_COS(4.*theta) - 429.*GET_COS(5.*theta) + 286.*GET_COS(6.*theta))*GET_SIN(theta/2.0))/512.0;
                    break;
                case 2:
                    fac = (GET_SQRT(17./CST_PI)*GET_POW(GET_COS(theta/2.0),4)*(7626. - 14454.*GET_COS(theta) + 12375.*GET_COS(2.*theta)
                                                                               - 9295.*GET_COS(3.*theta) + 6006.*GET_COS(4.*theta) - 3003.*GET_COS(5.*theta) + 1001.*GET_COS(6.*theta)))/512.0;
                    break;
                case 3:
                    fac = -(GET_SQRT(561./(2.0*CST_PI))*GET_POW(GET_COS(theta/2.0),5)*(-869. + 1660.*GET_COS(theta) - 1300.*GET_COS(2.*theta)
                                                                                       + 910.*GET_COS(3.*theta) - 455.*GET_COS(4.*theta) + 182.*GET_COS(5.*theta))*GET_SIN(theta/2.0))/128.0;
                    break;
                case 4:
                    fac = (GET_SQRT(935./(2.0*CST_PI))*GET_POW(GET_COS(theta/2.0),6)*(265. - 442.*GET_COS(theta) + 364.*GET_COS(2.*theta)
                                                                                      - 182.*GET_COS(3.*theta) + 91.*GET_COS(4.*theta))*GET_POW(GET_SIN(theta/2.0),2))/32.0;
                    break;
                case 5:
                    fac = -(GET_SQRT(12155./(2.0*CST_PI))*GET_POW(GET_COS(theta/2.0),7)*(-19. + 42.*GET_COS(theta) - 21.*GET_COS(2.*theta)
                                                                                         + 14.*GET_COS(3.*theta))*GET_POW(GET_SIN(theta/2.0),3))/8.0;
                    break;
                case 6:
                    fac = GET_SQRT(255255./CST_PI)*GET_POW(GET_COS(theta/2.0),8)*(-1. + 2.*GET_COS(theta))*GET_SIN(CST_PI/4.0 - theta/2.0)
                    *GET_SIN(CST_PI/4.0 + theta/2.0)*GET_POW(GET_SIN(theta/2.0),4);
                    break;
                case 7:
                    fac = -(GET_SQRT(17017./(2.0*CST_PI))*GET_POW(GET_COS(theta/2.0),9)*(-1. + 4.*GET_COS(theta))*GET_POW(GET_SIN(theta/2.0),5));
                    break;
                case 8:
                    fac = GET_SQRT(34034./CST_PI)*GET_POW(GET_COS(theta/2.0),10)*GET_POW(GET_SIN(theta/2.0),6);
                    break;
                default:
                    print_err("Error -%s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
                    return CEV_FAILURE;
                    break;
            }
        } /* l==8 */
        else
        {
            print_err("Error -%s: Unsupported mode l=%d (only l in [2,8] implemented)\n", __func__, l);
            return CEV_FAILURE;
        }
    }
    else
    {
        print_err("Error -%s: Unsupported mode s=%d (only s=-2 implemented)\n", __func__, __func__, s);
        return CEV_FAILURE;
    }
    if (m)
        ans = CX16polar(1.0, m*phi) * fac;
    else
        ans = fac;
    *ret = ans;
    return CEV_SUCCESS;
}

/**
 * Computes the n-th Jacobi polynomial for polynomial weights alpha and beta.
 * The implementation here is only valid for real x -- enforced by the argument
 * type. An extension to complex values would require evaluation of several
 * gamma functions.
 *
 * See http://en.wikipedia.org/wiki/Jacobi_polynomials
 */
double XLALJacobiPolynomial(int n, int alpha, int beta, double x)
{
    double f1 = (x-1)/2.0, f2 = (x+1)/2.0;
    int s=0;
    double sum=0, val=0;
    if( n == 0 ) return 1.0;
    for( s=0; n-s >= 0; s++ )
    {
        val=1.0;
        val *= gsl_sf_choose( n+alpha, s );
        val *= gsl_sf_choose( n+beta, n-s );
        if( n-s != 0 ) val *= pow( f1, n-s );
        if( s != 0 ) val*= pow( f2, s );

        sum += val;
    }
    return sum;
}

/**
 * These functions compute the amplitude and phase of a Wigner coefficient
 * Dlmmp, given Euler angles of an active rotation.
 */
double WignerdMatrix(
                    int l,        /**< mode number l */
                    int mp,        /**< mode number m' */
                    int m,        /**< mode number m */
                    double beta  /**< euler angle (rad) */
)
{

	int k = MIN( l+m, MIN( l-m, MIN( l+mp, l-mp )));
	double a=0, lam=0;
	if(k == l+m){
		a = mp-m;
		lam = mp-m;
	} else if(k == l-m) {
		a = m-mp;
		lam = 0;
	} else if(k == l+mp) {
		a = m-mp;
		lam = 0;
	} else if(k == l-mp) {
		a = mp-m;
		lam = mp-m;
	}

	int b = 2*l-2*k-a;
	double pref = pow(-1, lam) * sqrt(gsl_sf_choose( 2*l-k, k+a )) / sqrt(gsl_sf_choose( k+b, b ));

	return pref * pow(sin(beta/2.0), a) * pow( cos(beta/2.0), b) * XLALJacobiPolynomial(k, a, b, cos(beta));
}

/**
 * Computes the full Wigner D matrix for the Euler angle alpha, beta, and gamma
 * with major index 'l' and minor index transition from m to mp.
 *
 * Uses a slightly unconventional method since the intuitive version by Wigner
 * is less suitable to algorthmic development.
 *
 * See http://en.wikipedia.org/wiki/Wigner_D-matrix
 *
 * Currently only supports the modes which are implemented for the spin
 * weighted spherical harmonics.
 */
COMPLEX16 WignerDMatrix(
                        int l,        /**< mode number l */
                        int mp,        /**< mode number m' */
                        int m,        /**< mode number m */
                        double alpha,  /**< euler angle (rad) */
                        double beta, /**< euler angle (rad) */
                        double gam  /**< euler angle (rad) */
)
{
	 return cexp( -I*mp*alpha ) *
			WignerdMatrix( l, mp, m, beta ) * 
			cexp( -I*m*gam );
}
