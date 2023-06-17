/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 2019.07.03, Sydney
**/

#ifndef __INCLUDE_MYOPTPARSER__
#define __INCLUDE_MYOPTPARSER__
#include "myDatatypes.h"
#include <string.h>

#undef no_argument
#undef required_argument
#undef optional_argument

#define opt_no_argument        0
#define opt_required_argument    1
#define opt_optional_argument    2


typedef struct tagOPTION
{
    const CHAR *name;
    INT has_arg;
    INT *flag;
    INT val;
} OPTION;

INT getopt_long_only(INT argc, CHAR *const *argv, const CHAR *options,
                     const OPTION *opt, INT *opt_index);

#endif

