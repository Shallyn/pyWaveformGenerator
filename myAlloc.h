/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_MYALLOC__
#define __INCLUDE_MYALLOC__
#include <stddef.h>
#include "myDatatypes.h"
#define CHECK_MEM_LEAK 1

void *myMallocLong(size_t n, const char *file, int line);
void *myCallocLong(size_t m, size_t n, const char *file, int line);
void myFree(void *p);
void CheckMemoryLeak();

#if CHECK_MEM_LEAK
#define MYMalloc(n)     myMallocLong(n, __FILE__, __LINE__)
#define MYCalloc(m, n)     myCallocLong(m, n, __FILE__, __LINE__)
#define MYFree(p)       myFree(p)
#else
#define MYMalloc(n)     malloc(n)
#define MYCalloc(m, n)     calloc(m, n)
#define MYFree(p)       myFree(p)
#endif
#endif

