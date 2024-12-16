/**
 * Writer: Xiaolin.liu
 * xiaolin.liu@mail.bnu.edu.cn
 *
 * This head file define basic utils for the math lib.
 *
 * 2019.05.21, UWA
 **/

#ifndef __INCLUDE_MY_UTILS__
#define __INCLUDE_MY_UTILS__
#include "myAlloc.h"
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <gsl/gsl_math.h>

/** < MACRO DEF > **/
/* Internal */
#define ___PASTE(a,b) a##b
#define __PASTE(a,b) ___PASTE(a,b)
#define __UNIQUE_ID(prefix) __PASTE(__PASTE(__UNIQUE_ID_, prefix), __COUNTER__)

/* Color */
//LOG \033[1;34;47mLOG\033[0m
//WARNING \033[4;31mWARNING\033[0m

/* Math */
#undef	MAX
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))

#undef	MIN
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))

#undef	ABS
#define ABS(a)	   (((a) < 0) ? -(a) : (a))

#undef	CLAMP
#define CLAMP(x, low, high)  (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x)))

#define GET_MAX(a,b) MAX(a,b)
#define GET_MIN(a,b) MIN(a,b)
#define GET_ABS(a) ABS(a)
#define GET_POW(a,b) pow((a), (b))
#define GET_CBRT(a) pow((a), 1.0/3.0)
#define GET_EXP(a) exp(a)
#define GET_SIN(a) sin(a)
#define GET_COS(a) cos(a)
#define GET_SQRT(a) sqrt(a)
#define GET_LN(a) log(a)
#define GET_LOG10(a) log10(a)
#define GET_LOG2(a) GET_LN(a)/CST_LN2
#define GET_CEIL(a) ceil(a)
#define GET_FLOOR(a) floor(a)
#define GET_FMA(a,b,c) fma((a),(b),(c))
#define GET_ATAN(a) atan(a)
#define GET_ATAN2(a,b) atan2((a), (b))

#define C_ABS(a) cabs(a)
#define C_CONJUGATE(a) cong(a)
#define C_ANGLE(a) carg(a)
#define C_SQRT(a) csqrt(a)
#define C_EXP(a) cexp(a)
#define C_POW(a, b) cpow((a), (b))
#define C_SIN(a) csin(a)
#define C_COS(a) ccos(a)
#define C_REAL(a) creal(a)
#define C_IMAG(a) cimag(a)

#define IS_ODD(a) GSL_IS_ODD(a)
#define IS_EVEN(a) GSL_IS_EVEN(a)

#define _SWAP(a,b,tmp) \
do{typeof(a) tmp = (a);(a)=(b);(b)=tmp;} while(0)

#define SWAP(a,b) \
_SWAP(a,b,__UNIQUE_ID(tmp))

#define IS_REAL8_FAIL_NAN(val) IsREAL8FailNaN(val)
/** < Vector print > **/
#define _PRINT_LEN(x, len, fmt, itr) \
UINT itr; \
print_err("["); \
for(itr=0; itr < len; ++itr)\
{\
print_err(fmt, x[itr]); \
}\
print_err("]\n");

#define PRINT_LEN(x, len, fmt) \
_PRINT_LEN(x, len, fmt, __UNIQUE_ID(_itr))

#define PRINT_VEC(x, fmt) PRINT_LEN((x)->data, (x)->length, fmt)
#define PRINT_REAL8VEC(x) PRINT_VEC(x, " %f ")
#define PRINT_UINTVEC(x) PRINT_VEC(x, " %d ")


/** < FUNCTIONS DECLARATION > **/

/* Used for standart print */
INT print_err(const CHAR *fmt, ...);
INT print_log(CHAR *fmt, ...);
INT print_warning(CHAR *fmt, ...);
INT print_out(const CHAR *fmt, ...);
INT print_debug(CHAR *fmt, ...);
INT print_error(CHAR *fmt, ...);

/* Used for creating & destroying basic data struct */
REAL8Vector* CreateREAL8Vector(UINT length);
void DestroyREAL8Vector(REAL8Vector* vector);

CHARVector* CreateCHARVector(UINT length, UINT STR_LEN);
void DestroyCHARVector(CHARVector* vector);

COMPLEX16Vector* CreateCOMPLEX16Vector(UINT length);
void DestroyCOMPLEX16Vector(COMPLEX16Vector* vector);
void CalculateAmpPhaseFromCOMPLEX16Vector(COMPLEX16Vector *vector,
        REAL8Vector **ampVector, REAL8Vector **phaseVector);
void CalculateRIFromCOMPLEX16Vector(COMPLEX16Vector *vector,
        REAL8Vector **rVector, REAL8Vector **iVector);

UINTVector* CreateUINTVector(UINT length);
void DestroyUINTVector(UINTVector* vector);

#define STRUCTFREE(vec, struc) \
do{\
if((vec)){\
__PASTE(Destroy,struc)((vec));\
(vec)=NULL;}\
}while(0)


REAL8Array *CreateREAL8Array_comm(UINT *dims, UINT ndim);
REAL8Array *CreateREAL8Array(UINT ndim, ...);
void DestroyREAL8Array(REAL8Array* array);

REAL8TimeSeries *CreateREAL8TimeSeries (const REAL8 epoch, REAL8 deltaT, UINT length);
void DestroyREAL8TimeSeries( REAL8TimeSeries * series );

COMPLEX16TimeSeries *CreateCOMPLEX16TimeSeries (const REAL8 epoch, REAL8 deltaT, UINT length);
void DestroyCOMPLEX16TimeSeries( COMPLEX16TimeSeries * series );

REAL8FrequencySeries *CreateREAL8FrequencySeries(const REAL8 epoch, REAL8 f0, REAL8 deltaF, UINT length);
void DestroyREAL8FrequencySeries( REAL8FrequencySeries * series );

COMPLEX16FrequencySeries *CreateCOMPLEX16FrequencySeries(const REAL8 epoch, REAL8 f0, REAL8 deltaF, UINT length);
void DestroyCOMPLEX16FrequencySeries( COMPLEX16FrequencySeries * series );

VectorArray *CreateVectorArray_comm(UINT *dims, UINT ndim);
VectorArray *CreateVectorArray(UINT ndim, ...);
void DestroyVectorArray(VectorArray *varr);

MatrixArray *CreateMatrixArray_comm(UINT *dims, UINT ndim);
MatrixArray *CreateMatrixArray(UINT ndim, ...);
void DestroyMatrixArray(MatrixArray *marr);

/* REAL8Array */
INT getArrayValue_comm(REAL8 *val, REAL8Array *arr, UINT *cood);
INT getArrayValue(REAL8 *val, REAL8Array *arr, UINT ndim,...);
void setArrayValue_comm(REAL8 val, REAL8Array *arr, UINT *cood);
INT setArrayValue(REAL8 val, REAL8Array *arr, UINT ndim,...);

/* VectorArray */
void setVectorArrayValue_comm(REAL8 val, VectorArray *varr, UINT *cood);
INT setVectorArrayValue(REAL8 val, VectorArray *varr, UINT ndim, ...);

INT getVectorArrayValue_comm(REAL8 *val, VectorArray *varr, UINT *cood);
INT getVectorArrayValue(REAL8 *val, VectorArray *varr, UINT ndim, ...);

/* MatrixArray */

INT Array2gslMatrixArray(REAL8Array *arr, MatrixArray **marr);
INT gslMatrixArray2Array(MatrixArray *marr, REAL8Array **arr);
INT MatrixArrayTranspose(MatrixArray *marr);
INT getMatrixArrayValue_comm(REAL8 *val, MatrixArray *marr, UINT *cood);
INT getMatrixArrayValue(REAL8 *val, MatrixArray *marr, UINT ndim, ...);
void setMatrixArrayValue_comm(REAL8 val, MatrixArray *marr, UINT *cood);
INT setMatrixArrayValue(REAL8 val, MatrixArray *marr, UINT ndim, ...);

/* Math */
COMPLEX16 CX16polar(REAL8 r, REAL8 phi);

/* String */
char * StringCaseSubstring(const char *haystack, const char *needle);


#define TP_INT8_C(c) INT64_C(c)
#define REAL8_FAIL_NAN_INT TP_INT8_C(0x7ff80000000001a1)
#define IS_REAL8_FAIL_NAN(val) IsREAL8FailNaN(val)
#define REAL8_FAIL_NAN ( REAL8FailNaN() )
static inline REAL8 REAL8FailNaN(void)
{
    volatile const union {
        INT8 i;
        REAL8 x;
    } val = {
        REAL8_FAIL_NAN_INT};
    return val.x;
}

static inline int IsREAL8FailNaN(REAL8 val)
{
    volatile const union {
        INT8 i;
        unsigned char s[8];
    } a = {
        REAL8_FAIL_NAN_INT};
    volatile union {
        REAL8 x;
        unsigned char s[8];
    } b;
    size_t n;
    b.x = val;
    for (n = 0; n < sizeof(val); ++n)
        if (a.s[n] != b.s[n])
            return 0;
    return 1;
}

int SphericalToCartesian(REAL8 qCart[],
                         REAL8 pCart[],
                         const REAL8 qSph[],
                         const REAL8 pSph[]);

int CartesianToSpherical(REAL8 qSph[],
                         REAL8 pSph[],
                         const REAL8 qCart[],
                         const REAL8 pCart[]);
UINT argmax(REAL8Vector *vec);
REAL8Vector *get_slice(REAL8Vector *vec, UINT lo, UINT hi);
UINT FindClosestIndex(
    REAL8Vector *vec, /**<< Input: monotonically increasing vector */
    REAL8 value       /**<< Input: value to look for */
);
int EulerAnglesZYZFromRotationMatrixActive(REAL8* alpha, REAL8* beta, REAL8* gamma, REAL8Array* R);
int RotationMatrixActiveFromBasisVectors(REAL8Array* R, const REAL8 e1p[], const REAL8 e2p[], const REAL8 e3p[]);
int XLALREAL8VectorUnwrapAngle( REAL8Vector *out, const REAL8Vector *in );

COMPLEX16Vector *CutCOMPLEX16Vector (
	COMPLEX16Vector *sequence,
	size_t first,
	size_t length
);

COMPLEX16TimeSeries  *CutCOMPLEX16TimeSeries(
	const COMPLEX16TimeSeries *series,
	size_t first,
	size_t length
);


#endif




