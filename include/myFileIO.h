
#ifndef __INCLUDE_MYFIO__
#define __INCLUDE_MYFIO__
#include "myUtils.h"
// #include <hdf5.h>

INT get_REAL8TimeSeries_waveform(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, FILE *fp);
INT get_REAL8TimeSeries(REAL8TimeSeries **Series, FILE *fp);
INT get_COMPLEX16TimeSeries_waveform(COMPLEX16TimeSeries **hout, FILE *fp);

INT cmd_mkdir(CHAR *folderName);
INT read_waveform(REAL8Vector **time, REAL8Vector **hreal, REAL8Vector **himag, FILE *file);

// INT dump_to_extendible_H5_2D(const char* FILENAME, char *dname, REAL8Array
// *arr); INT writeREAL8Tohdf5(CHAR *fname, CHAR *dname, REAL8 val, INT
// is_delete); INT writeINTTohdf5(CHAR *fname, CHAR *dname, INT val, INT
// is_delete); INT DumpREAL8ArrayTohdf5(CHAR *fname, CHAR *dname, REAL8Array
// *array, INT is_delete); INT DumpREAL8VectorTohdf5(CHAR *fname, CHAR *dname,
// REAL8Vector *vec, INT is_delete);

#endif
