#ifndef SONSTIGE_H
#define SONSTIGE_H

#include "jurassic.h"
#include "scatter.h"
#include "forwardmodel.h"

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/* Initialize look-up tables. */
// void init_tbl(ctl_t *ctl,
//	      tbl_t *tbl);

/* Find array index. */
int locate(double *xx,  /* array */
	   int n,       /* array size */ 
	   double x);   /* value */

/* Read observation data. */
/* Reads observations e.g for retrieval */
void read_obs(const char *dirname,
	      const char *filename,
	      ctl_t *ctl,
	      obs_t *obs);

/* Write observation data. */
void write_obs(const char *dirname,
	       const char *filename,
	       ctl_t *ctl,
	       obs_t *obs);

#endif
