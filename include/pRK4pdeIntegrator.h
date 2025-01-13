/**
 * Writer: Xiaolin.liu
 * xiaolin.liu@mail.bnu.edu.cn
 *
 * This module contains basic functions for  calculation.
 * Functions list:
 * Kernel:
 * 20xx.xx.xx, LOC
 **/

#ifndef __INCLUDE_PRK4PDEINTEGRATOR__
#define __INCLUDE_PRK4PDEINTEGRATOR__

#include "pUtils.h"
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_spline.h>

typedef struct tagARKIntegrator
{
    gsl_odeiv_step *step;
    gsl_odeiv_control *control;
    gsl_odeiv_evolve *evolve;

    gsl_odeiv_system *sys;

    int (*dydt)(double t, const double y[], double dydt[], void *params);
    int (*stop)(double t, const double y[], double dydt[], void *params);

    int retries;        /* retries with smaller step when derivatives encounter
                           singularity */
    int stopontestonly; /* stop only on test, use tend to size buffers only */

    int returncode;
} ARKIntegrator;

ARKIntegrator *XLALAdaptiveRungeKutta4Init(
    int dim, int (*dydt)(double t, const double y[], double dydt[], void *params), /* These are XLAL functions! */
    int (*stop)(double t, const double y[], double dydt[], void *params), double eps_abs, double eps_rel);

void DestroyARKIntegrator(ARKIntegrator *integrator);

int XLALAdaptiveRungeKutta4NoInterpolate(ARKIntegrator *integrator, void *params, REAL8 *yinit, REAL8 tinit, REAL8 tend,
                                         REAL8 deltat_or_h0, REAL8 min_deltat_or_h0, REAL8Array **t_and_y_out);
int XLALAdaptiveRungeKutta4NoInterpolateWithDeriv(ARKIntegrator *integrator, void *params, REAL8 *yinit, REAL8 tinit,
                                                  REAL8 tend, REAL8 deltat_or_h0, REAL8 min_deltat_or_h0,
                                                  REAL8Array **t_and_y_and_dydt_out);
int XLALAdaptiveRungeKutta4NoInterpolateWithDerivPrec(ARKIntegrator *integrator, void *params, REAL8 *yinit,
                                                      REAL8 tinit, REAL8 tend, REAL8 deltat_or_h0,
                                                      REAL8 min_deltat_or_h0, REAL8Array **t_and_y_and_dydt_out);
int XLALAdaptiveRungeKutta4WithDerivPrec(ARKIntegrator *integrator, void *params, REAL8 *yinit, REAL8 tinit, REAL8 tend,
                                         REAL8 deltat_or_h0, REAL8 min_deltat_or_h0, REAL8Array **t_and_y_and_dydt_out);

int XLALAdaptiveRungeKutta4(ARKIntegrator *integrator, void *params, REAL8 *yinit, REAL8 tinit, REAL8 tend,
                            REAL8 deltat, REAL8Array **yout);
int XLALAdaptiveRungeKutta4WithDeriv(ARKIntegrator *integrator, void *params, REAL8 *yinit, REAL8 tinit, REAL8 tend,
                                     REAL8 deltat, REAL8Array **yout);

#endif
