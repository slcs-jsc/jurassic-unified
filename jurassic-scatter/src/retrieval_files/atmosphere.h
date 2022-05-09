#ifndef ATMOSPHERE_H
#define ATMOSPHERE_H

/* Module for atmosphere related functions. */

#include "jurassic.h"

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/* Compose state vector or parameter vector. */
size_t atm2x(ctl_t *ctl,
	     atm_t *atm,
	     aero_t *aero,
	     gsl_vector *x,
	     int *iqa,
	     int *ipa);

/* Add elements to state vector. */
void atm2x_help(atm_t *atm,
		double zmin,
		double zmax,
		double *value,
		int val_iqa,
		gsl_vector *x,
		int *iqa,
		int *ipa,
		size_t *n);

/* Write atmospheric data. */
void write_atm(const char *dirname,
	       const char *filename,
	       ctl_t *ctl,
	       atm_t *atm);

#endif