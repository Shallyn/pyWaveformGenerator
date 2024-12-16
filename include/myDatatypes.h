/**
 * Writer: Xiaolin.liu
 * xiaolin.liu@mail.bnu.edu.cn
 *
 * This head file define basic datatype in the math lib.
 * Referenced from LAL
 * Based on glib.h and comples.h
 *
 * 2019.05.21, UWA
 **/

#ifndef __INCLUDE_MY_DATA_TYPE__
#define __INCLUDE_MY_DATA_TYPE__

#include <math.h>
// #include <glib.h>
#include <complex.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <stdio.h>
#include <string.h>

#define MAX_STRLENGTH 256

#define DEFAULT_m1 10
#define DEFAULT_m2 10
#define DEFAULT_f_min 40
#define DEFAULT_sx 0.0
#define DEFAULT_sy 0.0
#define DEFAULT_sz 0.0
#define DEFAULT_deltat 16384
#define DEFAULT_distance 100
#define DEFAULT_inclination 0.0
#define DEFAULT_PHIREF 0.0
#define DEFAULT_ecc 0.0
#define DEFAULT_semilatus -1 // unused

/** < DATA TYPE DEF > **/
typedef int16_t INT2; /**< Two-byte signed integer */
typedef int32_t INT4; /**< Four-byte signed integer. */
typedef int64_t INT8; /**< Eight-byte signed integer; on some platforms this is
                         equivalent to <tt>long int</tt> instead. */
typedef signed short INT16;
typedef signed int INT32;

typedef unsigned char UINT8;
typedef unsigned short UINT16;
typedef unsigned int UINT32;

typedef short SINT;
typedef long LINT;
typedef int INT;
typedef unsigned int UINT;

typedef char CHAR;

typedef float REAL4;
typedef double REAL8;
typedef long double REAL16;

typedef double complex COMPLEX16;
typedef long double complex COMPLEX32;

typedef int BOOLEAN;

typedef void *POINTER;

#ifndef INT64_C
#define INT64_C(c) (c##LL)
#define UINT64_C(c) (c##ULL)
#endif

/** < CONSTANT DEF (CST) > **/
#define CST_E 2.7182818284590452353602874713526625           /** e */
#define CST_PI 3.141592653589793238462643383279502884        /** pi **/
#define CST_2PI 6.283185307179586476925286766559005768       /** 2pi **/
#define CST_PI_2 1.570796326794896619231321691639751442      /**< pi/2 */
#define CST_1_PI 0.318309886183790671537767526745028724      /** 1/pi **/
#define CST_PI_180 1.745329251994329576923690768488612713e-2 /** pi/180 **/
#define CST_180_PI 57.295779513082320876798154814105170332   /** 180/pi **/
#define CST_LN2 0.6931471805599453094172321214581766         /** ln2 **/
#define CST_GAMMA 0.5772156649015328606065120900824024       /**< gamma */
#define CST_INV_LOG2E 0.69314718055994530941723212145817656807550013436026
#define CST_SQ3 1.7320508075688772935274463415059
#define MAX_STR_LEN 256

#define CST_MSUN_SI 1.988409870698050731911960804878414216e30  /**< Solar mass, kg */
#define CST_PC_SI 3.085677581491367278913937957796471611e16    /**< Parsec, m */
#define CST_MRSUN_SI 1.476625038050124729627979840144936351e3  /**< Geometrized solar mass, m */
#define CST_MTSUN_SI 4.925490947641266978197229498498379006e-6 /**< Geometrized solar mass, s */

#define CST_G_SI 6.67384e-11 /**< Gravitational constant, N m^2 kg^-2 */
#define CST_C_SI 299792458e0 /**< Speed of light in vacuo, m s^-1 */

#ifndef FALSE
#define FALSE (0)
#endif

#ifndef TRUE
#define TRUE (!FALSE)
#endif

/** < STRUCTURE DEF > **/
typedef struct tagREAL8Vector
{
    UINT length; /** Number of elements **/
    REAL8 *data; /** Pointer to the data array **/
} REAL8Vector;

typedef struct tagCOMPLEX16Vector
{
    UINT length;
    COMPLEX16 *data;
} COMPLEX16Vector;

typedef struct tagUINTVector
{
    UINT length;
    UINT *data;
} UINTVector;

typedef struct tagREAL8Array
{
    UINTVector *dimLength; /** Vector of array dimensions **/
    REAL8 *data;           /** Pointer to the data array **/
    UINT size;             /** Number of data **/
} REAL8Array;

typedef struct tagCOMPLEX16Array
{
    UINTVector *dimLength; /** Vector of array dimensions **/
    COMPLEX16 *data;       /** Pointer to the data array **/
    UINT size;             /** Number of data **/
} COMPLEX16Array;

typedef struct tagREAL8TimeSeries
{
    REAL8 epoch;       /** The start time **/
    REAL8 deltaT;      /** The time step **/
    REAL8Vector *data; /** The sequence of sampled data. **/
} REAL8TimeSeries;

typedef struct tagCOMPLEX16TimeSeries
{
    REAL8 epoch;
    REAL8 deltaT;
    COMPLEX16Vector *data;
} COMPLEX16TimeSeries;

typedef struct tagREAL8FrequencySeries
{
    REAL8 epoch; /** The start time **/
    REAL8 f0;    /** Initial frequency **/
    REAL8 deltaF;
    REAL8Vector *data;
} REAL8FrequencySeries;

typedef struct tagCOMPLEX16FrequencySeries
{
    REAL8 epoch;
    REAL8 f0;
    REAL8 deltaF;
    COMPLEX16Vector *data;
} COMPLEX16FrequencySeries;

/** < gsl datatype negotiation > **/
typedef struct tagMatrixArray
{
    UINTVector *dimLength;
    gsl_matrix **data;
    UINT *mdim;
    UINT size;
} MatrixArray;

typedef struct tagVectorArray
{
    UINTVector *dimLength;
    gsl_vector **data;
    UINT vlength;
    UINT size;
} VectorArray;

typedef struct tagCHARVecotor
{
    UINT length;
    CHAR **data;
    UINT STR_LEN;
} CHARVector;

/** < CONTROL ERROR VALUE DEF (CEV) > **/
enum ControlErrorValue
{
    CEV_SUCCESS = GSL_SUCCESS, /** Success return value **/
    CEV_FAILURE = GSL_FAILURE, /** Failure return **/
    CEV_EFUNC = 1024,          /** Internal function call failed bit. **/
    CEV_ENOMEM = 12,           /** Memory allocation error **/
    CEV_EFAULT = 14,           /** Invalid pointer **/
    CEV_EINVAL = 22            /** Invalid argument **/
};
#endif
