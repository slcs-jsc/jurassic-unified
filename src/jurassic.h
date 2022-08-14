#ifndef JURASSIC_UNIFIED_H
#define JURASSIC_UNIFIED_H

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_statistics.h>

#ifdef MPI
#include <mpi.h>
#endif


/* ------------------------------------------------------------
   Macros...
   ------------------------------------------------------------ */

/* Allocate memory. */
#define ALLOC(ptr, type, n)       \
  if((ptr=malloc((size_t)(n)*sizeof(type)))==NULL)  \
  ERRMSG("Out of memory!");

/* Compute angle between two vectors. */
#define ANGLE(a, b)           \
  acos(GSL_MIN(GSL_MAX(DOTP(a, b)/NORM(a)/NORM(b), -1), 1))

/* Compute Cartesian distance between two vectors. */
#define DIST(a, b) sqrt(DIST2(a, b))

/* Compute squared distance between two vectors. */
#define DIST2(a, b)             \
  ((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]))

/* Compute dot product of two vectors. */
#define DOTP(a, b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])

/* Print error message and quit program. */
#define ERRMSG(msg) {                                                   \
  printf("\nError (%s, %s, l%d): %s\n\n",                             \
      __FILE__, __func__, __LINE__, msg);                      \
  exit(EXIT_FAILURE);                                                 \
}

/* Compute exponential interpolation. */
#define EXP(x0, y0, x1, y1, x)          \
  (((y0)>0 && (y1)>0)           \
   ? ((y0)*exp(log((y1)/(y0))/((x1)-(x0))*((x)-(x0))))          \
   : LIN(x0, y0, x1, y1, x))

/* Read binary data. */
#define FREAD(ptr, type, nmemb, stream) {                               \
  if(fread(ptr, sizeof(type), (size_t)nmemb, stream)!=(size_t)nmemb)  \
  ERRMSG("Error while reading!");                                   \
}

/* Write binary data. */
#define FWRITE(ptr, type, nmemb, stream) {        \
  if(fwrite(ptr, sizeof(type), (size_t)nmemb, stream)!=(size_t)nmemb) \
  ERRMSG("Error while writing!");         \
}

/* Compute linear interpolation. */
#define LIN(x0, y0, x1, y1, x)      \
  ((y0)+((y1)-(y0))/((x1)-(x0))*((x)-(x0)))

/* Execute netCDF library command and check result. */
#define NC(cmd) {            \
  if((cmd)!=NC_NOERR)            \
  ERRMSG(nc_strerror(cmd));          \
}

/* Compute norm of a vector. */
#define NORM(a) sqrt(DOTP(a, a))

/* Print macro for debugging. */
#define PRINT(format, var)                                              \
  printf("Print (%s, %s, l%d): %s= "format"\n",                         \
      __FILE__, __func__, __LINE__, #var, var);

//Added:

/*! Start or stop a timer. */
#define TIMER(name, mode) jur_timer(name, __FILE__, __func__, __LINE__, mode)

#define __deprecated__ __attribute__((deprecated))

/* stringify the value of a macro, two expansion levels needed */
#define xstr(a) str(a)
#define str(b) #b

/*! Read string tokens. */
// in scatter version there were 4 parameters
#define TOK_FIVE_ARGS(line, tok, format, var, saveptr) {      \
  if(((tok)=strtok_r((line), " \t",saveptr))) {     \
    if(sscanf(tok, format, &(var))!=1) continue;  \
  } else ERRMSG("Error while reading!");    \
}

/* Read string tokens. */
#define TOK(line, tok, format, var) {     \
  if(((tok)=strtok((line), " \t"))) {     \
    if(sscanf(tok, format, &(var))!=1) continue;  \
  } else ERRMSG("Error while reading!");    \
}

/* ------------------------------------------------------------
   Constants...
   ------------------------------------------------------------ */

/* First spectroscopic constant (c_1 = 2 h c^2) [W/(m^2 sr cm^-4)]. */
#define C1 1.19104259e-8

/* Second spectroscopic constant (c_2 = h c / k) [K/cm^-1]. */
#define C2 1.43877506

/* Standard pressure [hPa]. */
#define P0 1013.25

/* Mean radius of Earth [km]. */
#define RE 6367.421

/* Mass of Earth [kg]. */
#define ME 5.976e24

/* Temperature of the Sun [K]. */
#define TSUN 5780.

//Added:
/* Standard gravity [m/s^2]. */
#define G0 9.80665

/* Standard temperature [K]. */
#define T0 273.15

/* ------------------------------------------------------------
   Quantity indices...
   ------------------------------------------------------------ */

/* Pressure. */
#define IDXP 0

/* Temperature. */
#define IDXT 1

/* Volume mixing ratios. */
#define IDXQ(ig) (2+ig)

/* Extinction. */
#define IDXK(iw) (2+ctl->ng+iw)

/* Particle concentration. */
#define IDXNN  (2+ctl->ng+ctl->nw+1)

/* Particle size (radius). */
#define IDXRR  (2+ctl->ng+ctl->nw+2)

/* Particle size distribution width. */
#define IDXSS  (2+ctl->ng+ctl->nw+3)

#include "jurassic_dimensions.h"
#include "jurassic_structs.h"
#include "jurassic_functions.h"

#endif
