#ifndef FORWARDMODEL_H
#define FORWARDMODEL_H

#include "jurassic.h"
#include "scatter.h"
#include "atmosphere.h"

/* Compute brightness temperature. */
double brightness(double rad,
		  double nu);

void copy_obs_row(obs_t const *source, int rs, obs_t *dest, int rd);

void advanced_execute(ctl_t *ctl, atm_t *atm, aero_t *aero, queue_t *qs, int nr);

/* Determine ray paths and compute radiative transfer. */
void formod(ctl_t *ctl,
	    atm_t *atm,
	    obs_t *obs,
	    aero_t *aero);

/* Apply field of view convolution. */
void formod_fov(ctl_t *ctl,
		obs_t *obs);

/* Compute radiative transfer for a pencil beam. */
void formod_pencil(ctl_t *ctl,
		   atm_t *atm,
 		   obs_t *obs,
 		   aero_t *aero,
 		   int scattering,
 		   int ir,
       queue_t *q);

/* Compute Planck function. */
double planck(double t,
	      double nu);

/* Read shape function. */
void read_shape(const char *filename,
		double *x,
		double *y,
		int *n);

#endif
