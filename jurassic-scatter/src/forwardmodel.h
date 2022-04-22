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

/* Compute absorption coefficient of continua. */
//void formod_continua(ctl_t *ctl,
//         pos_t los[],
//		     int ip,
//		     double *beta);

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

/* Get transmittance from look-up tables. */
//void intpol_tbl(ctl_t *ctl,
//		tbl_t *tbl,
//    pos_t los[],
//		int ip,
//		double tau_path[NGMAX][NDMAX],
//		double tau_seg[NDMAX]);

/* Interpolate emissivity from look-up tables. */
//double intpol_tbl_eps(tbl_t *tbl,
//		      int ig,
//		      int id,
//		      int ip,
//		      int it,
//		      double u);

/* Interpolate column density from look-up tables. */
//double intpol_tbl_u(tbl_t *tbl,
//		    int ig,
//		    int id,
//		    int ip,
//		    int it,
//		    double eps);

/* Find array index in float array. */
//int locate_tbl(float *xx,
//	       int n,
//	       double x);

// needed for srcfunc_sca_sun from scatter.c
/* Compute Planck function. */
double planck(double t,
	      double nu);

// needed for formod_fov
/* Read shape function. */
void read_shape(const char *filename,
		double *x,
		double *y,
		int *n);

/* Read emissivity look-up tables. */
//void read_tbl(ctl_t *ctl,
//	      tbl_t *tbl);

/* Compute Planck source function. */
//void srcfunc_planck(ctl_t *ctl,
//		    double t,
//		    double *src);

#endif
