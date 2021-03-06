#ifndef FORWARDMODEL_H
#define FORWARDMODEL_H

#include "jurassic.h"
#include "sca_gpu_interface.h"
#include "sca_scatter.h"

void jur_sca_copy_obs_row(obs_t const *source, int rs, obs_t *dest, int rd);

void jur_sca_advanced_execute(ctl_t *ctl, atm_t *atm, aero_t *aero, queue_t *qs, int nr);

/* Determine ray paths and compute radiative transfer. */
void jur_sca_formod(ctl_t *ctl,
	    atm_t *atm,
	    obs_t *obs,
	    aero_t *aero);

/* Apply field of view convolution. */
void jur_sca_formod_fov(ctl_t *ctl,
		obs_t *obs);

/* Compute radiative transfer for a pencil beam. */
void jur_sca_formod_pencil(ctl_t *ctl,
		   atm_t *atm,
 		   obs_t *obs,
 		   aero_t *aero,
 		   int scattering,
 		   int ir,
       queue_t *q);

/* Read shape function. */
void jur_sca_read_shape(const char *filename,
		double *x,
		double *y,
		int *n);

#endif
