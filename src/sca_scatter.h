#ifndef SCATTER_H
#define SCATTER_H

/* Module containing all variables and functions related to scattering simulations */

#include "jurassic.h"
#include "sca_gpu_interface.h"
#include "sca_forwardmodel.h"

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/* Compute orthonormal basis. */
void jur_sca_bascoord(double *dz,
	      double *dy,
	      double *ex,
	      double *ey,
	      double *ez);

/* Compute Mie parameters. */
void jur_sca_bhmie(double x,
	   double n_real,
	   double n_imag,
	   double *phase,
	   double *qext,
	   double *qsca);

/* Gauss Hermite abcissas and weights. */
void jur_sca_gauher(double *x,
	    double *w);

/* Copy and initialize aerosol variables. */
void jur_sca_copy_aero(ctl_t *ctl,
	       aero_t *aero_dest,
	       aero_t *aero_src,
	       int init);

/* Get aerosol/cloud optical properties (1D). */
void jur_sca_get_opt_prop(ctl_t *ctl,
		  aero_t *aero);

/* Calculate optical properties with Mie theory for a log-normal mode. */
void jur_sca_opt_prop_mie_log(ctl_t *ctl,
		    aero_t *aero,
		    int count,
		    double *beta_ext,
		    double *beta_sca,
		    double phase[NDMAX][NTHETAMAX]);

/* Get optical properties from external database. - New */
void jur_sca_opt_prop_external(ctl_t *ctl,
		      aero_t *aero,
		      int count,
		      double *beta_ext,
		      double *beta_sca,
		      double phase[NDMAX][NTHETAMAX]);

/* Read aerosol/cloud data. */
void jur_sca_read_aero(const char *dirname,
	       const char *filename,
	       ctl_t *ctl,
	       aero_t *aero);

/* Compute scattering source. */
void jur_sca_srcfunc_sca(ctl_t *ctl,
		 atm_t *atm,
		 aero_t *aero,
		 double sec,
		 double *x,
		 double *dx,
		 int il,
		 double *src_sca,
		 int scattering,
     queue_t *q);

/* Compute scattering source (thermal emissions). */
void jur_sca_srcfunc_sca_1d(ctl_t *ctl,
		    atm_t *atm,
		    aero_t *aero,
        double sec,
		    double *x,
		    double *dx,
		    int il,
		    double *src_sca,
		    int scattering,
        queue_t *q);

/* Compute scattering source (thermal emissions). */
void jur_sca_srcfunc_sca_3d(ctl_t *ctl,
		    atm_t *atm,
		    aero_t *aero,
        double sec,
		    double *x,
		    double *dx,
		    int il,
		    double *src_sca,
		    int scattering,
        queue_t *q);

/* Add solar radiation to scattering source. */
void jur_sca_srcfunc_sca_sun(ctl_t *ctl,
		     atm_t *atm,
		     aero_t *aero,
		     double sec,
		     double *x,
		     double *dx,
		     int il,
		     double *src_sun,
         queue_t *q);

/* Compute Sun's angular coordinates. */
void jur_sca_suncoord(double sec,
	      double lon,
	      double lat,
	      double *azi,
	      double *sza);

/* Write particle data. */
void jur_sca_write_aero(const char *dirname,
		const char *filename,
		aero_t *aero);

#endif
