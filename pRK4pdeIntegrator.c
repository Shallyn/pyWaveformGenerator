/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "pRK4pdeIntegrator.h"
#include "myLog.h"
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_spline.h>
#define LOCAL_USE_DEBUG 0
#define XLAL_BEGINGSL \
        { \
          gsl_error_handler_t *saveGSLErrorHandler_; \
          saveGSLErrorHandler_ = gsl_set_error_handler_off();

#define XLAL_ENDGSL \
          gsl_set_error_handler( saveGSLErrorHandler_ ); \
        }

ARKIntegrator *XLALAdaptiveRungeKutta4Init(int dim, 
    int (*dydt) (double t, const double y[], double dydt[], void *params),   /* These are XLAL functions! */
    int (*stop) (double t, const double y[], double dydt[], void *params), 
    double eps_abs, double eps_rel)
{
    ARKIntegrator *integrator;

    /* allocate our custom integrator structure */
    if (!(integrator = (ARKIntegrator *) MYCalloc(1, sizeof(ARKIntegrator)))) 
    {
        return NULL;
    }

    /* allocate the GSL ODE components */
    integrator->step = gsl_odeiv_step_alloc(gsl_odeiv_step_rkf45, dim);
    //XLAL_CALLGSL(integrator->step = gsl_odeiv_step_alloc(gsl_odeiv_step_rk8pd, dim));
    integrator->control = gsl_odeiv_control_y_new(eps_abs, eps_rel);
    integrator->evolve = gsl_odeiv_evolve_alloc(dim);

    /* allocate the GSL system (functions, etc.) */
    integrator->sys = (gsl_odeiv_system *) MYCalloc(1, sizeof(gsl_odeiv_system));

    /* if something failed to be allocated, bail out */
    if (!(integrator->step) || !(integrator->control) || !(integrator->evolve) || !(integrator->sys)) 
    {
        STRUCTFREE(integrator, ARKIntegrator);
        return NULL;
    }

    integrator->dydt = dydt;
    integrator->stop = stop;

    integrator->sys->function = dydt;
    integrator->sys->jacobian = NULL;
    integrator->sys->dimension = dim;
    integrator->sys->params = NULL;

    integrator->retries = 6;
    integrator->stopontestonly = 0;

    return integrator;
}

void DestroyARKIntegrator(ARKIntegrator *integrator)
{
    if (!integrator)
        return;
    if (integrator->evolve)
        gsl_odeiv_evolve_free(integrator->evolve);
    if (integrator->control)
        gsl_odeiv_control_free(integrator->control);
    if (integrator->step)
        gsl_odeiv_step_free(integrator->step);
    MYFree(integrator->sys);
    MYFree(integrator);
    return;
}


int XLALAdaptiveRungeKutta4NoInterpolate(ARKIntegrator * integrator, 
                void * params, REAL8 * yinit, REAL8 tinit, REAL8 tend, 
                REAL8 deltat_or_h0, REAL8 min_deltat_or_h0,
                REAL8Array **t_and_y_out)
{
    int errnum = CEV_SUCCESS;
    int status; /* used throughout */
    unsigned int i, j;
    INT EOBversion = 2;
    /* needed for the integration */
    size_t dim, outputlength=0, bufferlength, retries;
    REAL8 t, tnew, h0, h0old;
    REAL8Array *buffers = NULL;
    REAL8 *temp = NULL, *y, *y0, *dydt_in, *dydt_in0, *dydt_out, *yerr; /* aliases */

    /* note: for speed, this replaces the single CALLGSL wrapper applied before each GSL call */
    XLAL_BEGINGSL;

    /* allocate the buffers!
        * note: REAL8Array has a field dimLength (UINT4Vector) with dimensions, and a field data that points to a single memory block;
        * dimLength itself has fields length and data */
    dim = integrator->sys->dimension;
    bufferlength = (int)((tend - tinit) / deltat_or_h0) + 2;   /* allow for the initial value and possibly a final semi-step */

    UINT dimn;/* Variable for different loop indices below */
    if(EOBversion==2) dimn = dim + 1;
    else dimn = dim + 4;//v3opt: Include three derivatives

    buffers = CreateREAL8Array(2, dimn/*dim + 1*/, bufferlength); /* 2-dimensional array, ((dim+1)) x bufferlength */

    temp = MYCalloc(6 * dim, sizeof(REAL8));

    if (!buffers || !temp) 
    {
        errnum = CEV_ENOMEM;
        goto bail_out;
    }

    y = temp;
    y0 = temp + dim;
    dydt_in = temp + 2 * dim;
    dydt_in0 = temp + 3 * dim;
    dydt_out = temp + 4 * dim;
    yerr = temp + 5 * dim;      /* aliases */

    /* set up to get started */
    integrator->sys->params = params;

    integrator->returncode = 0;

    retries = integrator->retries;

    t = tinit;
    h0 = deltat_or_h0;
    h0old = h0; /* initialized so that it will not trigger the check h0<h0old at the first step */
    memcpy(y, yinit, dim * sizeof(REAL8));

    /* store the first data point */
    buffers->data[0] = t;
    for ( i = 1; i <= dim; i++)
        buffers->data[i * bufferlength] = y[i - 1];

    /* compute derivatives at the initial time (dydt_in), bail out if impossible */
    if ((status = integrator->dydt(t, y, dydt_in, params)) != GSL_SUCCESS) 
    {
        integrator->returncode = status;
        errnum = CEV_FAILURE;
        goto bail_out;
    }

    if(EOBversion==3)
    {
        for ( i = 1; i <= 3; i++) //OPTV3: include the initial derivatives
            buffers->data[(dim+i)*bufferlength] = dydt_in[i-1];
    }

    UINT loop;/*variable for different loop indices below. */
    if(EOBversion==2) loop = dim;
    else loop = dim + 3;

    while (1) 
    {
        if (!integrator->stopontestonly && t >= tend)
            break;

        if (integrator->stop) 
        {
            if ((status = integrator->stop(t, y, dydt_in, params)) != GSL_SUCCESS) 
            {
                integrator->returncode = status;
                break;
            }
        }

        /* ready to try stepping! */
        try_step:

        /* if we would be stepping beyond the final time, stop there instead... */

#if LOCAL_USE_DEBUG
        if (!integrator->stopontestonly && t + h0 > tend)
            h0 = tend - t;
        PRINT_LOG_INFO(LOG_DEBUG, "INTEGRATE CURRENT STATE (t, y):");
        print_err("\t%g", t);
        for (i=0; i< dim; i++)
        {
            print_err("\t%g", y[i]);
        }
        print_err("\n");
#endif
        memcpy(y0, y, dim * sizeof(REAL8));     /* save y to y0, dydt_in to dydt_in0 */
        memcpy(dydt_in0, dydt_in, dim * sizeof(REAL8));

        /* call the GSL stepper function */
        status = gsl_odeiv_step_apply(integrator->step, t, h0, y, yerr, dydt_in, dydt_out, integrator->sys);
        /* note: If the user-supplied functions defined in the system dydt return a status other than GSL_SUCCESS,
            * the step will be aborted. In this case, the elements of y will be restored to their pre-step values,
            * and the error code from the user-supplied function will be returned. */

        /* did the stepper report a derivative-evaluation error? */
        if (status != GSL_SUCCESS) 
        {
            if (retries--) 
            {
                h0 = h0 / 10.0; /* if we have singularity retries left, reduce the timestep and try again */
                goto try_step;
            } 
            else 
            {
                integrator->returncode = status;
                break;  /* otherwise exit the loop */
            }
        } 
        else 
        {
            retries = integrator->retries;      /* we stepped successfully, reset the singularity retries */
        }

        tnew = t + h0;

        /* call the GSL error-checking function */
        status = gsl_odeiv_control_hadjust(integrator->control, integrator->step, y, yerr, dydt_out, &h0);

        /* Enforce a minimal allowed time step */
        /* To ignore this type of constraint, set min_deltat_or_h0 = 0 */
        if (h0 < min_deltat_or_h0) h0 = min_deltat_or_h0;
        // if (status == GSL_ODEIV_HADJ_INC)
        //     print_debug("h_old = %.16e, h_new = %.16e\n", h0old, h0);

        /* did the error-checker reduce the stepsize?
            * note: previously, was using status == GSL_ODEIV_HADJ_DEC
            * other possible return codes are GSL_ODEIV_HADJ_INC if it was increased,
            * GSL_ODEIV_HADJ_NIL if it was unchanged
            * since we introduced a minimal step size we simply compare to the saved value of h0 */
        /* if (status == GSL_ODEIV_HADJ_DEC) { */
        if (h0 < h0old) 
        {
            h0old = h0;
            memcpy(y, y0, dim * sizeof(REAL8)); /* if so, undo the step, and try again */
            memcpy(dydt_in, dydt_in0, dim * sizeof(REAL8));
            goto try_step;
        }

        /* update the current time and input derivatives, save the time step */
        t = tnew;
        h0old = h0;
        memcpy(dydt_in, dydt_out, dim * sizeof(REAL8));
        outputlength++;

        /* check if interpolation buffers need to be extended */
        if (outputlength >= bufferlength) 
        {
            REAL8Array *rebuffers;

            /* sadly, we cannot use ResizeREAL8Array, because it would only work if we extended the first array dimension,
                * so we need to copy everything over and switch the buffers. Oh well. */
            if (!(rebuffers = CreateREAL8Array(2, dimn, 2 * bufferlength))) 
            {
                errnum = CEV_ENOMEM;   /* ouch, that hurt */
                goto bail_out;
            } 
            else 
            {
                for ( i = 0; i <= loop /* dim */; i++)
                        memcpy(&rebuffers->data[i * 2 * bufferlength], &buffers->data[i * bufferlength], outputlength * sizeof(REAL8));
                STRUCTFREE(buffers, REAL8Array);
                buffers = rebuffers;
                bufferlength *= 2;
            }
        }

        /* copy time and state into output buffers */
        buffers->data[outputlength] = t;
        for ( i = 1; i <= loop; i++)
            buffers->data[i * bufferlength + outputlength] = y[i - 1];   /* y does not have time */
        if(EOBversion==3)
        {
            for ( i = 1; i <= 3; i++)
                buffers->data[(dim+i) * bufferlength + outputlength] = dydt_out[i - 1];  //OPTV3: Include 3 derivatives
        }
    } // End loop
    /* copy the final state into yinit */
    PRINT_LOG_INFO(LOG_DEBUG, "loop = %d, dim = %zu, OutputLen = %zu", loop, dim, outputlength);
    memcpy(yinit, y, dim * sizeof(REAL8));

    /* if we have completed at least one step, allocate the GSL interpolation object and the output array */
    if (outputlength == 0)
        goto bail_out;

    if(EOBversion==3)
    {
        outputlength++;
        (*t_and_y_out) = CreateREAL8Array(2,(dim + 4),outputlength);//OPTV3: Include derivatives
    }
    else
    {
        (*t_and_y_out) = CreateREAL8Array(2,(dim + 3),outputlength); //OPTV3: Include derivatives
    }

    if (!(*t_and_y_out)) 
    {
        errnum = CEV_ENOMEM;   /* ouch again, ran out of memory */
        STRUCTFREE(*t_and_y_out, REAL8Array);
        outputlength = 0;
        goto bail_out;
    }

    for( j=0;j<outputlength;j++) 
    {
        (*t_and_y_out)->data[j] = buffers->data[j];
        for( i=1;i<=loop;i++) 
        {
            (*t_and_y_out)->data[i*outputlength + j] = buffers->data[i*bufferlength + j];
        }
    }
    /* deallocate stuff and return */
    bail_out:
    XLAL_ENDGSL;

    // if (buffers)
    //     XLALDestroyREAL8Array(buffers); /* let's be careful, although all these checks may not be needed */
    STRUCTFREE(buffers, REAL8Array);
    if (temp)
        MYFree(temp);

    if (errnum)
        return -1;

    return outputlength;
}

int XLALAdaptiveRungeKutta4NoInterpolateWithDeriv(ARKIntegrator * integrator, 
                void * params, REAL8 * yinit, REAL8 tinit, REAL8 tend, 
                REAL8 deltat_or_h0, REAL8 min_deltat_or_h0,
                REAL8Array **t_and_y_and_dydt_out)
{
    int errnum = CEV_SUCCESS;
    int status; /* used throughout */
    unsigned int i, j;
    /* needed for the integration */
    size_t dim, outputlength=0, bufferlength, retries;
    REAL8 t, tnew, h0, h0old;
    REAL8Array *buffers = NULL;
    REAL8 *temp = NULL, *y, *y0, *dydt_in, *dydt_in0, *dydt_out, *yerr; /* aliases */

    /* note: for speed, this replaces the single CALLGSL wrapper applied before each GSL call */
    XLAL_BEGINGSL;

    /* allocate the buffers!
        * note: REAL8Array has a field dimLength (UINT4Vector) with dimensions, and a field data that points to a single memory block;
        * dimLength itself has fields length and data */
    dim = integrator->sys->dimension;
    bufferlength = (int)((tend - tinit) / deltat_or_h0) + 2;   /* allow for the initial value and possibly a final semi-step */

    UINT dimn;/* Variable for different loop indices below */
    UINT nextra = 3;
    dimn = 1 + dim + 2 + nextra; // SEOBSA: t, (r, ph, prT, pphi), (dr, dph), (H, dpr, dpphi)
    // dimn = dim + 1;
    buffers = CreateREAL8Array(2, dimn, bufferlength); 

    temp = MYCalloc(6 * dim, sizeof(REAL8));

    if (!buffers || !temp) 
    {
        errnum = CEV_ENOMEM;
        goto bail_out;
    }

    y = temp;
    y0 = temp + dim;
    dydt_in = temp + 2 * dim;
    dydt_in0 = temp + 3 * dim;
    dydt_out = temp + 4 * dim;
    yerr = temp + 5 * dim;      /* aliases */

    /* set up to get started */
    integrator->sys->params = params;
    SpinEOBParams *seobpms = (SpinEOBParams *) params;

    integrator->returncode = 0;

    retries = integrator->retries;

    t = tinit;
    h0 = deltat_or_h0;
    h0old = h0; /* initialized so that it will not trigger the check h0<h0old at the first step */
    memcpy(y, yinit, dim * sizeof(REAL8));

    /* store the first data point */
    buffers->data[0] = t;
    for ( i = 1; i <= dim; i++)
        buffers->data[i * bufferlength] = y[i - 1];

    /* compute derivatives at the initial time (dydt_in), bail out if impossible */
    if ((status = integrator->dydt(t, y, dydt_in, params)) != GSL_SUCCESS) 
    {
        integrator->returncode = status;
        errnum = CEV_FAILURE;
        goto bail_out;
    }

    for (i = 1; i <= 2; i++)
        buffers->data[(dim+i)*bufferlength] = dydt_in[i-1];

    for (i=1; i<=nextra; i++)
        buffers->data[(i+dim+2)*bufferlength] = seobpms->cache[i-1];

    // UINT loop;/*variable for different loop indices below. */
    // if(EOBversion==2) loop = dim;
    // else loop = dim + 3;
    PRINT_LOG_INFO(LOG_INFO, "Run RK4 integrator...");
    while (1) 
    {
        if (!integrator->stopontestonly && t >= tend)
            break;

        if (integrator->stop) 
        {
            if ((status = integrator->stop(t, y, dydt_in, params)) != GSL_SUCCESS) 
            {
                integrator->returncode = status;
                break;
            }
        }

        /* ready to try stepping! */
try_step:

        /* if we would be stepping beyond the final time, stop there instead... */

#if LOCAL_USE_DEBUG
        if (!integrator->stopontestonly && t + h0 > tend)
            h0 = tend - t;
        PRINT_LOG_INFO(LOG_DEBUG, "INTEGRATE CURRENT STATE (t, y):");
        print_err("\t%g", t);
        for (i=0; i< dim; i++)
        {
            print_err("\t%g", y[i]);
        }
        print_err("\n");
#endif
        memcpy(y0, y, dim * sizeof(REAL8));     /* save y to y0, dydt_in to dydt_in0 */
        memcpy(dydt_in0, dydt_in, dim * sizeof(REAL8));

        /* call the GSL stepper function */
        status = gsl_odeiv_step_apply(integrator->step, t, h0, y, yerr, dydt_in, dydt_out, integrator->sys);
        /* note: If the user-supplied functions defined in the system dydt return a status other than GSL_SUCCESS,
            * the step will be aborted. In this case, the elements of y will be restored to their pre-step values,
            * and the error code from the user-supplied function will be returned. */

        /* did the stepper report a derivative-evaluation error? */
        if (status != GSL_SUCCESS) 
        {
            if (retries--) 
            {
                h0 = h0 / 10.0; /* if we have singularity retries left, reduce the timestep and try again */
                goto try_step;
            } 
            else 
            {
                integrator->returncode = status;
                break;  /* otherwise exit the loop */
            }
        } 
        else 
        {
            retries = integrator->retries;      /* we stepped successfully, reset the singularity retries */
        }

        tnew = t + h0;

        /* call the GSL error-checking function */
        status = gsl_odeiv_control_hadjust(integrator->control, integrator->step, y, yerr, dydt_out, &h0);

        /* Enforce a minimal allowed time step */
        /* To ignore this type of constraint, set min_deltat_or_h0 = 0 */
        if (h0 < min_deltat_or_h0) h0 = min_deltat_or_h0;
        // if (status == GSL_ODEIV_HADJ_INC)
        //     print_debug("h_old = %.16e, h_new = %.16e\n", h0old, h0);

        /* did the error-checker reduce the stepsize?
            * note: previously, was using status == GSL_ODEIV_HADJ_DEC
            * other possible return codes are GSL_ODEIV_HADJ_INC if it was increased,
            * GSL_ODEIV_HADJ_NIL if it was unchanged
            * since we introduced a minimal step size we simply compare to the saved value of h0 */
        /* if (status == GSL_ODEIV_HADJ_DEC) { */
        if (h0 < h0old) 
        {
            h0old = h0;
            // print_debug("[%zu/%zu]h_old = %.16e, h_new = %.16e\n", outputlength, bufferlength, h0old, h0);
            memcpy(y, y0, dim * sizeof(REAL8)); /* if so, undo the step, and try again */
            memcpy(dydt_in, dydt_in0, dim * sizeof(REAL8));
            goto try_step;
        }

        /* update the current time and input derivatives, save the time step */
        t = tnew;
        h0old = h0;
        memcpy(dydt_in, dydt_out, dim * sizeof(REAL8));
        outputlength++;

        /* check if interpolation buffers need to be extended */
#if 1
        if (outputlength >= bufferlength) 
        {
            REAL8Array *rebuffers;

            /* sadly, we cannot use ResizeREAL8Array, because it would only work if we extended the first array dimension,
                * so we need to copy everything over and switch the buffers. Oh well. */
            if (!(rebuffers = CreateREAL8Array(2, dimn, 2 * bufferlength))) 
            {
                errnum = CEV_ENOMEM;   /* ouch, that hurt */
                goto bail_out;
            } 
            else 
            {
                for ( i = 0; i < dimn /* dim */; i++)
                    memcpy(&rebuffers->data[i * 2 * bufferlength], &buffers->data[i * bufferlength], outputlength * sizeof(REAL8));
                STRUCTFREE(buffers, REAL8Array);
                buffers = rebuffers;
                bufferlength *= 2;
            }
        }

        /* copy time and state into output buffers */
        buffers->data[outputlength] = t;
        for ( i = 1; i <= dim; i++)
        {
            buffers->data[i * bufferlength + outputlength] = y[i - 1];   /* y does not have time */
        }
        for (i=1; i<= 2; i++)
            buffers->data[(i+dim)*bufferlength + outputlength] = dydt_out[i-1];
        for (i=1; i<=nextra; i++)
            buffers->data[(i+dim+2)*bufferlength + outputlength] = seobpms->cache[i-1];
#endif

    } // End loop
    /* copy the final state into yinit */
    PRINT_LOG_INFO(LOG_DEBUG, "dim = %zu, OutputLen = %zu", dim, outputlength);
    memcpy(yinit, y, dim * sizeof(REAL8));

    /* if we have completed at least one step, allocate the GSL interpolation object and the output array */
    if (outputlength == 0)
        goto bail_out;

    // outputlength++;
#if 1
    (*t_and_y_and_dydt_out) = CreateREAL8Array(2, dimn, outputlength);//OPTV3: Include derivatives

    if (!(*t_and_y_and_dydt_out)) 
    {
        errnum = CEV_ENOMEM;   /* ouch again, ran out of memory */
        STRUCTFREE(*t_and_y_and_dydt_out, REAL8Array);
        outputlength = 0;
        goto bail_out;
    }

    for( j=0;j<outputlength;j++) 
    {
        // (*t_and_y_and_dydt_out)->data[j] = buffers->data[j];
        for( i=0;i<dimn;i++) 
        {
            (*t_and_y_and_dydt_out)->data[i*outputlength + j] = buffers->data[i*bufferlength + j];
            // print_out("%.16e", buffers->data[i*bufferlength + j]);
            // if (i==dimn-1)
            //     print_out("\n");
            // else
            //     print_out("\t");
        }

    }
#endif
    /* deallocate stuff and return */
    bail_out:
    XLAL_ENDGSL;

    // if (buffers)
    //     XLALDestroyREAL8Array(buffers); /* let's be careful, although all these checks may not be needed */
    STRUCTFREE(buffers, REAL8Array);
    if (temp)
        MYFree(temp);

    if (errnum)
        return -1;

    return outputlength;
}

int XLALAdaptiveRungeKutta4NoInterpolateWithDerivPrec(ARKIntegrator * integrator, 
                void * params, REAL8 * yinit, REAL8 tinit, REAL8 tend, 
                REAL8 deltat_or_h0, REAL8 min_deltat_or_h0,
                REAL8Array **t_and_y_and_dydt_out)
{
    int errnum = CEV_SUCCESS;
    int status; /* used throughout */
    unsigned int i, j;
    /* needed for the integration */
    size_t dim, outputlength=0, bufferlength, retries;
    REAL8 t, tnew, h0, h0old;
    REAL8Array *buffers = NULL;
    REAL8 *temp = NULL, *y, *y0, *dydt_in, *dydt_in0, *dydt_out, *yerr; /* aliases */

    /* note: for speed, this replaces the single CALLGSL wrapper applied before each GSL call */
    XLAL_BEGINGSL;

    /* allocate the buffers!
        * note: REAL8Array has a field dimLength (UINT4Vector) with dimensions, and a field data that points to a single memory block;
        * dimLength itself has fields length and data */
    dim = integrator->sys->dimension;
    bufferlength = (int)((tend - tinit) / deltat_or_h0) + 2;   /* allow for the initial value and possibly a final semi-step */

    UINT dimn;/* Variable for different loop indices below */
    UINT nextra = 3, nderiv = 3;
    dimn = dim + nderiv + 1 + nextra; // SEOBPrec: t, (x,y,z,pTx,pTy,pTz,s1x,s1y,s1z,s2x,s2y,s2z,phiM,phi), (vx,vy,vz), H, flux, prTDot
    // dimn = dim + 1;
    buffers = CreateREAL8Array(2, dimn, bufferlength); 

    temp = MYCalloc(6 * dim, sizeof(REAL8));

    if (!buffers || !temp) 
    {
        errnum = CEV_ENOMEM;
        goto bail_out;
    }

    y = temp;
    y0 = temp + dim;
    dydt_in = temp + 2 * dim;
    dydt_in0 = temp + 3 * dim;
    dydt_out = temp + 4 * dim;
    yerr = temp + 5 * dim;      /* aliases */

    /* set up to get started */
    integrator->sys->params = params;
    SpinEOBParams *seobpms = (SpinEOBParams *) params;

    integrator->returncode = 0;

    retries = integrator->retries;

    t = tinit;
    h0 = deltat_or_h0;
    h0old = h0; /* initialized so that it will not trigger the check h0<h0old at the first step */
    memcpy(y, yinit, dim * sizeof(REAL8));

    /* store the first data point */
    buffers->data[0] = t;
    for ( i = 1; i <= dim; i++)
        buffers->data[i * bufferlength] = y[i - 1];
// print_debug("y0 = (%.5e, %.5e, %.5e), (%.5e, %.5e, %.5e)\n", 
//     y[0], y[1], y[2], y[3], y[4], y[5]);
    /* compute derivatives at the initial time (dydt_in), bail out if impossible */
    if ((status = integrator->dydt(t, y, dydt_in, params)) != GSL_SUCCESS) 
    {
        integrator->returncode = status;
        errnum = CEV_FAILURE;
        goto bail_out;
    }
    for (i = 1; i <= nderiv; i++)
        buffers->data[(dim+i)*bufferlength] = dydt_in[i-1];

    for (i=1; i<=nextra; i++)
        buffers->data[(i+dim+nderiv)*bufferlength] = seobpms->cache[i-1];

    // UINT loop;/*variable for different loop indices below. */
    // if(EOBversion==2) loop = dim;
    // else loop = dim + 3;
    // INT debug_used = 0;
    PRINT_LOG_INFO(LOG_INFO, "Run RK4 integrator...");
    while (1) 
    {
        if (!integrator->stopontestonly && t >= tend)
            break;

        if (integrator->stop) 
        {
            if ((status = integrator->stop(t, y, dydt_in, params)) != GSL_SUCCESS) 
            {
                integrator->returncode = status;
                break;
            }
        }

        /* ready to try stepping! */
try_step:

        /* if we would be stepping beyond the final time, stop there instead... */

#if 0
        // if (!integrator->stopontestonly && t + h0 > tend)
        //     h0 = tend - t;
        if (debug_used == 0){
        PRINT_LOG_INFO(LOG_DEBUG, "INTEGRATE CURRENT STATE (t, y):");
        print_err("\t%g", t);
        for (i=0; i< dim; i++)
        {
            print_err("\t%g", y[i]);
        }
        print_err("\n");
        debug_used = 1;
        }
#endif
        memcpy(y0, y, dim * sizeof(REAL8));     /* save y to y0, dydt_in to dydt_in0 */
        memcpy(dydt_in0, dydt_in, dim * sizeof(REAL8));

// DEBUG_START;
        /* call the GSL stepper function */
        status = gsl_odeiv_step_apply(integrator->step, t, h0, y, yerr, dydt_in, dydt_out, integrator->sys);
        /* note: If the user-supplied functions defined in the system dydt return a status other than GSL_SUCCESS,
            * the step will be aborted. In this case, the elements of y will be restored to their pre-step values,
            * and the error code from the user-supplied function will be returned. */
// DEBUG_END;
        /* did the stepper report a derivative-evaluation error? */
        if (status != GSL_SUCCESS) 
        {
            if (retries--) 
            {
                h0 = h0 / 10.0; /* if we have singularity retries left, reduce the timestep and try again */
                goto try_step;
            } 
            else 
            {
                integrator->returncode = status;
                break;  /* otherwise exit the loop */
            }
        } 
        else 
        {
            retries = integrator->retries;      /* we stepped successfully, reset the singularity retries */
        }

        tnew = t + h0;

        /* call the GSL error-checking function */
        status = gsl_odeiv_control_hadjust(integrator->control, integrator->step, y, yerr, dydt_out, &h0);

        /* Enforce a minimal allowed time step */
        /* To ignore this type of constraint, set min_deltat_or_h0 = 0 */
        if (h0 < min_deltat_or_h0) h0 = min_deltat_or_h0;
        // if (status == GSL_ODEIV_HADJ_INC)
        //     print_debug("h_old = %.16e, h_new = %.16e\n", h0old, h0);

        /* did the error-checker reduce the stepsize?
            * note: previously, was using status == GSL_ODEIV_HADJ_DEC
            * other possible return codes are GSL_ODEIV_HADJ_INC if it was increased,
            * GSL_ODEIV_HADJ_NIL if it was unchanged
            * since we introduced a minimal step size we simply compare to the saved value of h0 */
        /* if (status == GSL_ODEIV_HADJ_DEC) { */
        if (h0 < h0old) 
        {
            h0old = h0;
            // print_debug("[%zu/%zu]h_old = %.16e, h_new = %.16e\n", outputlength, bufferlength, h0old, h0);
            memcpy(y, y0, dim * sizeof(REAL8)); /* if so, undo the step, and try again */
            memcpy(dydt_in, dydt_in0, dim * sizeof(REAL8));
            goto try_step;
        }

        /* update the current time and input derivatives, save the time step */
        t = tnew;
        h0old = h0;
        memcpy(dydt_in, dydt_out, dim * sizeof(REAL8));
        outputlength++;

        /* check if interpolation buffers need to be extended */
#if 1
        if (outputlength >= bufferlength) 
        {
            REAL8Array *rebuffers;

            /* sadly, we cannot use ResizeREAL8Array, because it would only work if we extended the first array dimension,
                * so we need to copy everything over and switch the buffers. Oh well. */
            if (!(rebuffers = CreateREAL8Array(2, dimn, 2 * bufferlength))) 
            {
                errnum = CEV_ENOMEM;   /* ouch, that hurt */
                goto bail_out;
            } 
            else 
            {
                for ( i = 0; i < dimn /* dim */; i++)
                    memcpy(&rebuffers->data[i * 2 * bufferlength], &buffers->data[i * bufferlength], outputlength * sizeof(REAL8));
                STRUCTFREE(buffers, REAL8Array);
                buffers = rebuffers;
                bufferlength *= 2;
            }
        }

        /* copy time and state into output buffers */
        buffers->data[outputlength] = t;
        for ( i = 1; i <= dim; i++)
        {
            buffers->data[i * bufferlength + outputlength] = y[i - 1];   /* y does not have time */
        }
        for (i=1; i<= nderiv; i++)
            buffers->data[(i+dim)*bufferlength + outputlength] = dydt_out[i-1];
        for (i=1; i<=nextra; i++)
            buffers->data[(i+dim+nderiv)*bufferlength + outputlength] = seobpms->cache[i-1];

#endif

    } // End loop

    // some thing wrong happens during the integratioin
    if (integrator->returncode == GSL_ERANGE)
    {
        errnum = CEV_FAILURE;
        goto bail_out;
    }
    /* copy the final state into yinit */
    PRINT_LOG_INFO(LOG_DEBUG, "dim = %zu, OutputLen = %zu", dim, outputlength);
    memcpy(yinit, y, dim * sizeof(REAL8));

    /* if we have completed at least one step, allocate the GSL interpolation object and the output array */
    if (outputlength == 0)
        goto bail_out;

    // outputlength++;
#if 1
    (*t_and_y_and_dydt_out) = CreateREAL8Array(2, dimn, outputlength);//OPTV3: Include derivatives

    if (!(*t_and_y_and_dydt_out)) 
    {
        errnum = CEV_ENOMEM;   /* ouch again, ran out of memory */
        STRUCTFREE(*t_and_y_and_dydt_out, REAL8Array);
        outputlength = 0;
        goto bail_out;
    }

    for( j=0;j<outputlength;j++) 
    {
        // (*t_and_y_and_dydt_out)->data[j] = buffers->data[j];
        for( i=0;i<dimn;i++) 
        {
            (*t_and_y_and_dydt_out)->data[i*outputlength + j] = buffers->data[i*bufferlength + j];
            // print_out("%.16e", buffers->data[i*bufferlength + j]);
            // if (i==dimn-1)
            //     print_out("\n");
            // else
            //     print_out("\t");
        }

    }
#endif
    /* deallocate stuff and return */
    bail_out:
    XLAL_ENDGSL;

    // if (buffers)
    //     XLALDestroyREAL8Array(buffers); /* let's be careful, although all these checks may not be needed */
    STRUCTFREE(buffers, REAL8Array);
    if (temp)
        MYFree(temp);

    if (errnum)
        return -1;

    return outputlength;
}

int XLALAdaptiveRungeKutta4WithDerivPrec(ARKIntegrator * integrator, 
                void * params, REAL8 * yinit, REAL8 tinit, REAL8 tend, 
                REAL8 deltat_or_h0, REAL8 min_deltat_or_h0,
                REAL8Array **t_and_y_and_dydt_out)
{
    int errnum = CEV_SUCCESS;
    int status; /* used throughout */
    unsigned int i, j;
    /* needed for the integration */
    size_t dim, outputlength=0, bufferlength, retries;
    REAL8 t, tnew, h0, h0old;
    REAL8Array *buffers = NULL;
    REAL8 *temp = NULL, *y, *y0, *dydt_in, *dydt_in0, *dydt_out, *yerr; /* aliases */

    /* note: for speed, this replaces the single CALLGSL wrapper applied before each GSL call */

    /* needed for the final interpolation */
    gsl_spline *interp = NULL;
    gsl_interp_accel *accel = NULL;
    int outputlen = 0;
    REAL8Array *output = NULL;
    REAL8 *times, *vector;      /* aliases */

    XLAL_BEGINGSL;

    /* allocate the buffers!
        * note: REAL8Array has a field dimLength (UINT4Vector) with dimensions, and a field data that points to a single memory block;
        * dimLength itself has fields length and data */
    dim = integrator->sys->dimension;
    bufferlength = (int)((tend - tinit) / deltat_or_h0) + 2;   /* allow for the initial value and possibly a final semi-step */

    UINT dimn;/* Variable for different loop indices below */
    UINT nextra = 3, nderiv = 3;
    dimn = dim + nderiv + 1 + nextra; // SEOBPrec: t, (x,y,z,pTx,pTy,pTz,s1x,s1y,s1z,s2x,s2y,s2z,phiM,phi), (vx,vy,vz), H, flux, prTDot
    // dimn = dim + 1;
    buffers = CreateREAL8Array(2, dimn, bufferlength); 

    temp = MYCalloc(6 * dim, sizeof(REAL8));

    if (!buffers || !temp) 
    {
        errnum = CEV_ENOMEM;
        goto bail_out;
    }

    y = temp;
    y0 = temp + dim;
    dydt_in = temp + 2 * dim;
    dydt_in0 = temp + 3 * dim;
    dydt_out = temp + 4 * dim;
    yerr = temp + 5 * dim;      /* aliases */

    /* set up to get started */
    integrator->sys->params = params;
    SpinEOBParams *seobpms = (SpinEOBParams *) params;

    integrator->returncode = 0;

    retries = integrator->retries;

    t = tinit;
    h0 = deltat_or_h0;
    h0old = h0; /* initialized so that it will not trigger the check h0<h0old at the first step */
    memcpy(y, yinit, dim * sizeof(REAL8));

    /* store the first data point */
    buffers->data[0] = t;
    for ( i = 1; i <= dim; i++)
        buffers->data[i * bufferlength] = y[i - 1];

    /* compute derivatives at the initial time (dydt_in), bail out if impossible */
    if ((status = integrator->dydt(t, y, dydt_in, params)) != GSL_SUCCESS) 
    {
        integrator->returncode = status;
        errnum = CEV_FAILURE;
        goto bail_out;
    }

    for (i = 1; i <= nderiv; i++)
        buffers->data[(dim+i)*bufferlength] = dydt_in[i-1];

    for (i=1; i<=nextra; i++)
        buffers->data[(i+dim+nderiv)*bufferlength] = seobpms->cache[i-1];

    // UINT loop;/*variable for different loop indices below. */
    // if(EOBversion==2) loop = dim;
    // else loop = dim + 3;
    PRINT_LOG_INFO(LOG_INFO, "Run RK4 integrator...");
    while (1) 
    {
        if (!integrator->stopontestonly && t >= tend)
            break;

        if (integrator->stop) 
        {
            if ((status = integrator->stop(t, y, dydt_in, params)) != GSL_SUCCESS) 
            {
                integrator->returncode = status;
                break;
            }
        }

        /* ready to try stepping! */
try_step:

        /* if we would be stepping beyond the final time, stop there instead... */

#if LOCAL_USE_DEBUG
        if (!integrator->stopontestonly && t + h0 > tend)
            h0 = tend - t;
        PRINT_LOG_INFO(LOG_DEBUG, "INTEGRATE CURRENT STATE (t, y):");
        print_err("\t%g", t);
        for (i=0; i< dim; i++)
        {
            print_err("\t%g", y[i]);
        }
        print_err("\n");
#endif
        memcpy(y0, y, dim * sizeof(REAL8));     /* save y to y0, dydt_in to dydt_in0 */
        memcpy(dydt_in0, dydt_in, dim * sizeof(REAL8));

        /* call the GSL stepper function */
        status = gsl_odeiv_step_apply(integrator->step, t, h0, y, yerr, dydt_in, dydt_out, integrator->sys);
        /* note: If the user-supplied functions defined in the system dydt return a status other than GSL_SUCCESS,
            * the step will be aborted. In this case, the elements of y will be restored to their pre-step values,
            * and the error code from the user-supplied function will be returned. */

        /* did the stepper report a derivative-evaluation error? */
        if (status != GSL_SUCCESS) 
        {
            if (retries--) 
            {
                h0 = h0 / 10.0; /* if we have singularity retries left, reduce the timestep and try again */
                goto try_step;
            } 
            else 
            {
                integrator->returncode = status;
                break;  /* otherwise exit the loop */
            }
        } 
        else 
        {
            retries = integrator->retries;      /* we stepped successfully, reset the singularity retries */
        }

        tnew = t + h0;

        /* call the GSL error-checking function */
        status = gsl_odeiv_control_hadjust(integrator->control, integrator->step, y, yerr, dydt_out, &h0);

        /* Enforce a minimal allowed time step */
        /* To ignore this type of constraint, set min_deltat_or_h0 = 0 */
        if (h0 < min_deltat_or_h0) h0 = min_deltat_or_h0;
        // if (status == GSL_ODEIV_HADJ_INC)
        //     print_debug("h_old = %.16e, h_new = %.16e\n", h0old, h0);

        /* did the error-checker reduce the stepsize?
            * note: previously, was using status == GSL_ODEIV_HADJ_DEC
            * other possible return codes are GSL_ODEIV_HADJ_INC if it was increased,
            * GSL_ODEIV_HADJ_NIL if it was unchanged
            * since we introduced a minimal step size we simply compare to the saved value of h0 */
        /* if (status == GSL_ODEIV_HADJ_DEC) { */
        if (h0 < h0old) 
        {
            h0old = h0;
            // print_debug("[%zu/%zu]h_old = %.16e, h_new = %.16e\n", outputlength, bufferlength, h0old, h0);
            memcpy(y, y0, dim * sizeof(REAL8)); /* if so, undo the step, and try again */
            memcpy(dydt_in, dydt_in0, dim * sizeof(REAL8));
            goto try_step;
        }

        /* update the current time and input derivatives, save the time step */
        t = tnew;
        h0old = h0;
        memcpy(dydt_in, dydt_out, dim * sizeof(REAL8));
        outputlength++;

        /* check if interpolation buffers need to be extended */
#if 1
        if (outputlength >= bufferlength) 
        {
            REAL8Array *rebuffers;

            /* sadly, we cannot use ResizeREAL8Array, because it would only work if we extended the first array dimension,
                * so we need to copy everything over and switch the buffers. Oh well. */
            if (!(rebuffers = CreateREAL8Array(2, dimn, 2 * bufferlength))) 
            {
                errnum = CEV_ENOMEM;   /* ouch, that hurt */
                goto bail_out;
            } 
            else 
            {
                for ( i = 0; i < dimn /* dim */; i++)
                    memcpy(&rebuffers->data[i * 2 * bufferlength], &buffers->data[i * bufferlength], outputlength * sizeof(REAL8));
                STRUCTFREE(buffers, REAL8Array);
                buffers = rebuffers;
                bufferlength *= 2;
            }
        }

        /* copy time and state into output buffers */
        buffers->data[outputlength] = t;
        for ( i = 1; i <= dim; i++)
        {
            buffers->data[i * bufferlength + outputlength] = y[i - 1];   /* y does not have time */
        }
        for (i=1; i<= nderiv; i++)
            buffers->data[(i+dim)*bufferlength + outputlength] = dydt_out[i-1];
        for (i=1; i<=nextra; i++)
            buffers->data[(i+dim+nderiv)*bufferlength + outputlength] = seobpms->cache[i-1];

#endif

    } // End loop
    /* copy the final state into yinit */
    PRINT_LOG_INFO(LOG_DEBUG, "dim = %zu, OutputLen = %zu", dim, outputlength);
    memcpy(yinit, y, dim * sizeof(REAL8));

    /* if we have completed at least one step, allocate the GSL interpolation object and the output array */
    if (outputlength == 0)
        goto bail_out;
    interp = gsl_spline_alloc(gsl_interp_cspline, outputlength + 1);
    accel = gsl_interp_accel_alloc();
    outputlen = (int)(t / deltat_or_h0) + 1;
    output = CreateREAL8Array(2, dimn, outputlen);
    if (!interp || !accel || !output) {
        errnum = CEV_ENOMEM;   /* ouch again, ran out of memory */
        STRUCTFREE(output, REAL8Array);
        outputlen = 0;
        goto bail_out;
    }
    /* make an array of times */
    times = output->data;
    for ( j = 0; j < outputlen; j++)
        times[j] = deltat_or_h0 * j;

    /* interpolate! */
    for ( i = 1; i <= dimn-1; i++) {
        gsl_spline_init(interp, &buffers->data[0], &buffers->data[bufferlength * i], outputlength + 1);

        vector = output->data + outputlen * i;
        for ( j = 0; j < outputlen; j++) {
            gsl_spline_eval_e(interp, times[j], accel, &(vector[j]));
        }
    }

    // outputlength++;
#if 0
    (*t_and_y_and_dydt_out) = CreateREAL8Array(2, dimn, outputlength);//OPTV3: Include derivatives

    if (!(*t_and_y_and_dydt_out)) 
    {
        errnum = CEV_ENOMEM;   /* ouch again, ran out of memory */
        STRUCTFREE(*t_and_y_and_dydt_out, REAL8Array);
        outputlength = 0;
        goto bail_out;
    }

    for( j=0;j<outputlength;j++) 
    {
        // (*t_and_y_and_dydt_out)->data[j] = buffers->data[j];
        for( i=0;i<dimn;i++) 
        {
            (*t_and_y_and_dydt_out)->data[i*outputlength + j] = buffers->data[i*bufferlength + j];
            // print_out("%.16e", buffers->data[i*bufferlength + j]);
            // if (i==dimn-1)
            //     print_out("\n");
            // else
            //     print_out("\t");
        }

    }
#endif
    /* deallocate stuff and return */
    bail_out:
    XLAL_ENDGSL;

    // if (buffers)
    //     XLALDestroyREAL8Array(buffers); /* let's be careful, although all these checks may not be needed */
    STRUCTFREE(buffers, REAL8Array);
    if (temp)
        MYFree(temp);
    if (interp)
        gsl_spline_free(interp);
    if (accel)
        gsl_interp_accel_free(accel);
    if (errnum)
        return -1;
    *t_and_y_and_dydt_out = output;
    return outputlen;
}


int XLALAdaptiveRungeKutta4(ARKIntegrator * integrator,
        void *params, REAL8 * yinit, REAL8 tinit, REAL8 tend, 
        REAL8 deltat, REAL8Array ** yout)
{
// print_debug("END time = %f\n", tend);
    int errnum = CEV_SUCCESS;
    int status; /* used throughout */
    unsigned int i, j;
    /* needed for the integration */
    size_t dim, bufferlength, cnt, retries;
    REAL8 t, tnew, h0;
    REAL8Array *buffers = NULL;
    REAL8 *temp = NULL, *y, *y0, *dydt_in, *dydt_in0, *dydt_out, *yerr; /* aliases */

    /* needed for the final interpolation */
    gsl_spline *interp = NULL;
    gsl_interp_accel *accel = NULL;
    int outputlen = 0;
    REAL8Array *output = NULL;
    REAL8 *times, *vector;      /* aliases */

    /* note: for speed, this replaces the single CALLGSL wrapper applied before each GSL call */
    XLAL_BEGINGSL;

    /* allocate the buffers!
     * note: REAL8Array has a field dimLength (UINT4Vector) with dimensions, and a field data that points to a single memory block;
     * dimLength itself has fields length and data */
    dim = integrator->sys->dimension;
    bufferlength = (int)((tend - tinit) / deltat) + 2;  /* allow for the initial value and possibly a final semi-step */
    buffers = CreateREAL8Array(2, dim + 1, bufferlength);  /* 2-dimensional array, (dim+1) x bufferlength */
    temp = MYCalloc(6 * dim, sizeof(REAL8));

    if (!buffers || !temp) 
    {
        errnum = CEV_ENOMEM;
        goto bail_out;
    }

    y = temp;
    y0 = temp + dim;
    dydt_in = temp + 2 * dim;
    dydt_in0 = temp + 3 * dim;
    dydt_out = temp + 4 * dim;
    yerr = temp + 5 * dim;      /* aliases */

    /* set up to get started */
    integrator->sys->params = params;

    integrator->returncode = 0;

    cnt = 0;
    retries = integrator->retries;

    t = tinit;
    h0 = deltat;
    memcpy(y, yinit, dim * sizeof(REAL8));

    /* store the first data point */
    buffers->data[0] = t;
    for ( i = 1; i <= dim; i++)
        buffers->data[i * bufferlength] = y[i - 1];

    /* compute derivatives at the initial time (dydt_in), bail out if impossible */
    if ((status = integrator->dydt(t, y, dydt_in, params)) != GSL_SUCCESS) 
    {
        integrator->returncode = status;
        errnum = CEV_FAILURE;
        PRINT_LOG_INFO(LOG_CRITICAL, "First step of integration failed.");
        goto bail_out;
    }

    while (1) 
    {
// print_debug("[%d]prec: %e/%e\n", integrator->stopontestonly, t, tend);
        if (!integrator->stopontestonly && t >= tend) 
        {
            break;
        }

        if (integrator->stop) 
        {
            if ((status = integrator->stop(t, y, dydt_in, params)) != GSL_SUCCESS) {
                integrator->returncode = status;
                break;
            }
        }

        /* ready to try stepping! */
        try_step:
#if LOCAL_USE_DEBUG
        if (!integrator->stopontestonly && t + h0 > tend)
            h0 = tend - t;
        PRINT_LOG_INFO(LOG_DEBUG, "INTEGRATE CURRENT STATE (t, y):");
        print_err("\t%g", t);
        for (i=0; i< dim; i++)
        {
            print_err("\t%g", y[i]);
        }
        print_err("\n");
#endif

        /* if we would be stepping beyond the final time, stop there instead... */
        if (!integrator->stopontestonly && t + h0 > tend)
            h0 = tend - t;
        memcpy(y0, y, dim * sizeof(REAL8));     /* save y to y0, dydt_in to dydt_in0 */
        memcpy(dydt_in0, dydt_in, dim * sizeof(REAL8));

        /* call the GSL stepper function */
        status = gsl_odeiv_step_apply(integrator->step, t, h0, y, yerr, dydt_in, dydt_out, integrator->sys);
        /* note: If the user-supplied functions defined in the system dydt return a status other than GSL_SUCCESS,
         * the step will be aborted. In this case, the elements of y will be restored to their pre-step values,
         * and the error code from the user-supplied function will be returned. */

        /* did the stepper report a derivative-evaluation error? */
        if (status != GSL_SUCCESS) 
        {
            if (retries--) 
            {
                h0 = h0 / 10.0; /* if we have singularity retries left, reduce the timestep and try again */
                goto try_step;
            } 
            else 
            {
                integrator->returncode = status;
                break;  /* otherwise exit the loop */
            }
        } 
        else 
        {
            retries = integrator->retries;      /* we stepped successfully, reset the singularity retries */
        }

        tnew = t + h0;

        /* call the GSL error-checking function */
        status = gsl_odeiv_control_hadjust(integrator->control, integrator->step, y, yerr, dydt_out, &h0);

        /* did the error-checker reduce the stepsize?
         * note: other possible return codes are GSL_ODEIV_HADJ_INC if it was increased,
         * GSL_ODEIV_HADJ_NIL if it was unchanged */
        if (status == GSL_ODEIV_HADJ_DEC) 
        {
            memcpy(y, y0, dim * sizeof(REAL8)); /* if so, undo the step, and try again */
            memcpy(dydt_in, dydt_in0, dim * sizeof(REAL8));
            goto try_step;
        }

        /* update the current time and input derivatives */
        t = tnew;
        memcpy(dydt_in, dydt_out, dim * sizeof(REAL8));
        cnt++;

        /* check if interpolation buffers need to be extended */
        if (cnt >= bufferlength) 
        {
            REAL8Array *rebuffers;

            /* sadly, we cannot use ResizeREAL8Array, because it would only work if we extended the first array dimension,
             * so we need to copy everything over and switch the buffers. Oh well. */
            if (!(rebuffers = CreateREAL8Array(2, dim + 1, 2 * bufferlength))) 
            {
                errnum = CEV_ENOMEM;   /* ouch, that hurt */
                PRINT_LOG_INFO(LOG_CRITICAL, "Could not resize REAL8Array");
                goto bail_out;
            } 
            else 
            {
                for (i = 0; i <= dim; i++)
                    memcpy(&rebuffers->data[i * 2 * bufferlength], &buffers->data[i * bufferlength], cnt * sizeof(REAL8));
                STRUCTFREE(buffers, REAL8Array);
                buffers = rebuffers;
                bufferlength *= 2;
            }
        }

        /* copy time and state into interpolation buffers */
        buffers->data[cnt] = t;
        for (i = 1; i <= dim; i++)
            buffers->data[i * bufferlength + cnt] = y[i - 1];   /* y does not have time */
    }

    /* copy the final state into yinit */

    memcpy(yinit, y, dim * sizeof(REAL8));

    /* if we have completed at least one step, allocate the GSL interpolation object and the output array */
    if (cnt == 0)
        goto bail_out;

    interp = gsl_spline_alloc(gsl_interp_cspline, cnt + 1);
    accel = gsl_interp_accel_alloc();

    outputlen = (int)(t / deltat) + 1;
    output = CreateREAL8Array(2, dim + 1, outputlen);

    if (!interp || !accel || !output) {
        errnum = CEV_ENOMEM;   /* ouch again, ran out of memory */
        STRUCTFREE(output, REAL8Array);
        outputlen = 0;
        goto bail_out;
    }

    /* make an array of times */
    times = output->data;
    for ( j = 0; j < outputlen; j++)
        times[j] = tinit + deltat * j;

    /* interpolate! */
    for ( i = 1; i <= dim; i++) {
        gsl_spline_init(interp, &buffers->data[0], &buffers->data[bufferlength * i], cnt + 1);

        vector = output->data + outputlen * i;
        for ( j = 0; j < outputlen; j++) {
            gsl_spline_eval_e(interp, times[j], accel, &(vector[j]));
        }
    }

    /* deallocate stuff and return */
bail_out:

    XLAL_ENDGSL;

    // if (buffers)
    //     XLALDestroyREAL8Array(buffers); /* let's be careful, although all these checks may not be needed */
    STRUCTFREE(buffers, REAL8Array);
    if (temp)
        MYFree(temp);

    if (interp)
        gsl_spline_free(interp);
    if (accel)
        gsl_interp_accel_free(accel);

    if (errnum)
        return -1;
    
    *yout = output;
    return outputlen;
}


int XLALAdaptiveRungeKutta4WithDeriv(ARKIntegrator * integrator,
        void *params, REAL8 * yinit, REAL8 tinit, REAL8 tend, 
        REAL8 deltat, REAL8Array ** yout)
{
// print_debug("END time = %f\n", tend);
    int errnum = CEV_SUCCESS;
    int status; /* used throughout */
    unsigned int i, j;
    /* needed for the integration */
    size_t dim, nextra, bufferlength, cnt, retries;
    REAL8 t, tnew, h0;
    REAL8Array *buffers = NULL;
    REAL8 *temp = NULL, *y, *y0, *dydt_in, *dydt_in0, *dydt_out, *yerr; /* aliases */

    /* needed for the final interpolation */
    gsl_spline *interp = NULL;
    gsl_interp_accel *accel = NULL;
    int outputlen = 0;
    REAL8Array *output = NULL;
    REAL8 *times, *vector;      /* aliases */

    /* note: for speed, this replaces the single CALLGSL wrapper applied before each GSL call */
    XLAL_BEGINGSL;

    /* allocate the buffers!
     * note: REAL8Array has a field dimLength (UINT4Vector) with dimensions, and a field data that points to a single memory block;
     * dimLength itself has fields length and data */
    dim = integrator->sys->dimension;
    nextra = 3;
    size_t dimn = 1 + dim + 2 + nextra;
    bufferlength = (int)((tend - tinit) / deltat) + 2;  /* allow for the initial value and possibly a final semi-step */
    buffers = CreateREAL8Array(2, dimn, bufferlength);  /* 2-dimensional array, (dim+1) x bufferlength */
    temp = MYCalloc(6 * dim, sizeof(REAL8));

    if (!buffers || !temp) 
    {
        errnum = CEV_ENOMEM;
        goto bail_out;
    }

    y = temp;
    y0 = temp + dim;
    dydt_in = temp + 2 * dim;
    dydt_in0 = temp + 3 * dim;
    dydt_out = temp + 4 * dim;
    yerr = temp + 5 * dim;      /* aliases */

    /* set up to get started */
    integrator->sys->params = params;
    SpinEOBParams *seobpms = (SpinEOBParams *)params;
    integrator->returncode = 0;

    cnt = 0;
    retries = integrator->retries;

    t = tinit;
    h0 = deltat;
    memcpy(y, yinit, dim * sizeof(REAL8));

    /* store the first data point */
    buffers->data[0] = t;
    for ( i = 1; i <= dim; i++)
        buffers->data[i * bufferlength] = y[i - 1];

    /* compute derivatives at the initial time (dydt_in), bail out if impossible */
    if ((status = integrator->dydt(t, y, dydt_in, params)) != GSL_SUCCESS) 
    {
        integrator->returncode = status;
        errnum = CEV_FAILURE;
        PRINT_LOG_INFO(LOG_CRITICAL, "First step of integration failed.");
        goto bail_out;
    }

    for ( i = 1; i <= 2; i++)
        buffers->data[(i+dim) * bufferlength] = dydt_in[i - 1];

    for ( i=1; i<=nextra; i++)
        buffers->data[(i+dim+2) * bufferlength] = seobpms->cache[i-1];

    while (1) 
    {
// print_debug("[%d]prec: %e/%e\n", integrator->stopontestonly, t, tend);
        if (!integrator->stopontestonly && t >= tend) 
        {
            break;
        }

        if (integrator->stop) 
        {
            if ((status = integrator->stop(t, y, dydt_in, params)) != GSL_SUCCESS) {
                integrator->returncode = status;
                break;
            }
        }

        /* ready to try stepping! */
        try_step:
#if LOCAL_USE_DEBUG
        if (!integrator->stopontestonly && t + h0 > tend)
            h0 = tend - t;
        PRINT_LOG_INFO(LOG_DEBUG, "INTEGRATE CURRENT STATE (t, y):");
        print_err("\t%g", t);
        for (i=0; i< dim; i++)
        {
            print_err("\t%g", y[i]);
        }
        print_err("\n");
#endif

        /* if we would be stepping beyond the final time, stop there instead... */
        if (!integrator->stopontestonly && t + h0 > tend)
            h0 = tend - t;
        memcpy(y0, y, dim * sizeof(REAL8));     /* save y to y0, dydt_in to dydt_in0 */
        memcpy(dydt_in0, dydt_in, dim * sizeof(REAL8));

        /* call the GSL stepper function */
        status = gsl_odeiv_step_apply(integrator->step, t, h0, y, yerr, dydt_in, dydt_out, integrator->sys);
        /* note: If the user-supplied functions defined in the system dydt return a status other than GSL_SUCCESS,
         * the step will be aborted. In this case, the elements of y will be restored to their pre-step values,
         * and the error code from the user-supplied function will be returned. */

        /* did the stepper report a derivative-evaluation error? */
        if (status != GSL_SUCCESS) 
        {
            if (retries--) 
            {
                h0 = h0 / 10.0; /* if we have singularity retries left, reduce the timestep and try again */
                goto try_step;
            } 
            else 
            {
                integrator->returncode = status;
                break;  /* otherwise exit the loop */
            }
        } 
        else 
        {
            retries = integrator->retries;      /* we stepped successfully, reset the singularity retries */
        }

        tnew = t + h0;

        /* call the GSL error-checking function */
        status = gsl_odeiv_control_hadjust(integrator->control, integrator->step, y, yerr, dydt_out, &h0);

        /* did the error-checker reduce the stepsize?
         * note: other possible return codes are GSL_ODEIV_HADJ_INC if it was increased,
         * GSL_ODEIV_HADJ_NIL if it was unchanged */
        if (status == GSL_ODEIV_HADJ_DEC) 
        {
            memcpy(y, y0, dim * sizeof(REAL8)); /* if so, undo the step, and try again */
            memcpy(dydt_in, dydt_in0, dim * sizeof(REAL8));
            goto try_step;
        }

        /* update the current time and input derivatives */
        t = tnew;
        memcpy(dydt_in, dydt_out, dim * sizeof(REAL8));
        cnt++;

        /* check if interpolation buffers need to be extended */
        if (cnt >= bufferlength) 
        {
            REAL8Array *rebuffers;

            /* sadly, we cannot use ResizeREAL8Array, because it would only work if we extended the first array dimension,
             * so we need to copy everything over and switch the buffers. Oh well. */
            if (!(rebuffers = CreateREAL8Array(2, dimn, 2 * bufferlength))) 
            {
                errnum = CEV_ENOMEM;   /* ouch, that hurt */
                PRINT_LOG_INFO(LOG_CRITICAL, "Could not resize REAL8Array");
                goto bail_out;
            } 
            else 
            {
                for (i = 0; i < dimn; i++)
                    memcpy(&rebuffers->data[i * 2 * bufferlength], &buffers->data[i * bufferlength], cnt * sizeof(REAL8));

                STRUCTFREE(buffers, REAL8Array);
                buffers = rebuffers;
                bufferlength *= 2;
            }
        }

        /* copy time and state into interpolation buffers */
        buffers->data[cnt] = t;
        for (i = 1; i <= dim; i++)
            buffers->data[i * bufferlength + cnt] = y[i - 1];   /* y does not have time */
        for ( i = 1; i <= 2; i++)
            buffers->data[(i+dim) * bufferlength + cnt] = dydt_out[i - 1];
        for ( i=1; i<=nextra; i++)
            buffers->data[(i+dim+2) * bufferlength + cnt] = seobpms->cache[i-1];
    }

    /* copy the final state into yinit */

    memcpy(yinit, y, dim * sizeof(REAL8));

    /* if we have completed at least one step, allocate the GSL interpolation object and the output array */
    if (cnt == 0)
        goto bail_out;

    interp = gsl_spline_alloc(gsl_interp_cspline, cnt + 1);
    accel = gsl_interp_accel_alloc();

    outputlen = (int)(t / deltat) + 1;
    output = CreateREAL8Array(2, dimn, outputlen);

    if (!interp || !accel || !output) {
        errnum = CEV_ENOMEM;   /* ouch again, ran out of memory */
        STRUCTFREE(output, REAL8Array);
        outputlen = 0;
        goto bail_out;
    }

    /* make an array of times */
    times = output->data;
    for ( j = 0; j < outputlen; j++)
        times[j] = deltat * j;

    /* interpolate! */
    for ( i = 1; i <= dimn-1; i++) {
        gsl_spline_init(interp, &buffers->data[0], &buffers->data[bufferlength * i], cnt + 1);

        vector = output->data + outputlen * i;
        for ( j = 0; j < outputlen; j++) {
            gsl_spline_eval_e(interp, times[j], accel, &(vector[j]));
        }
    }

    // for (i=0; i<outputlen; i++)
    // {
    //     for (j=0; j<dimn; j++)
    //     {
    //         vector = output->data + outputlen * j;
    //         print_out("%.16e", vector[i]);
    //         if (j==dimn-1)
    //             print_out("\n");
    //         else
    //             print_out("\t");
    //     }
    // }

    /* deallocate stuff and return */
bail_out:

    XLAL_ENDGSL;

    // if (buffers)
    //     XLALDestroyREAL8Array(buffers); /* let's be careful, although all these checks may not be needed */
    STRUCTFREE(buffers, REAL8Array);
    if (temp)
        MYFree(temp);

    if (interp)
        gsl_spline_free(interp);
    if (accel)
        gsl_interp_accel_free(accel);

    if (errnum)
        return -1;
    
    *yout = output;
    return outputlen;
}
