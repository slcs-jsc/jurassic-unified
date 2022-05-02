#include "atmosphere.h"

#define __host__
#include "interface.h"

/*****************************************************************************/

size_t atm2x(ctl_t *ctl,
	     atm_t *atm,
	     aero_t *aero,
	     gsl_vector *x,
	     int *iqa,
	     int *ipa) {
  
  int ig, iw;
  
  size_t n=0;
  
  /* Add pressure... */
  atm2x_help(atm, ctl->retp_zmin, ctl->retp_zmax,
	     atm->p, IDXP, x, iqa, ipa, &n);
  
  /* Add temperature... */
  atm2x_help(atm, ctl->rett_zmin, ctl->rett_zmax,
	     atm->t, IDXT, x, iqa, ipa, &n);
  
  /* Add volume mixing ratios... */
  for(ig=0; ig<ctl->ng; ig++)
    atm2x_help(atm, ctl->retq_zmin[ig], ctl->retq_zmax[ig],
	       atm->q[ig], IDXQ(ig), x, iqa, ipa, &n);
  
  /* Add extinction... */
  for(iw=0; iw<ctl->nw; iw++)
    atm2x_help(atm, ctl->retk_zmin[iw], ctl->retk_zmax[iw],
	       atm->k[iw], IDXK(iw), x, iqa, ipa, &n);
 
  /* Add particle concentration... */
  if(ctl->retnn) {
    if(x!=NULL)
      gsl_vector_set(x, n, aero->nn[0]);
    if(iqa!=NULL)
      iqa[n]=IDXNN;
    if(ipa!=NULL)
      ipa[n]=-1;
    (n)++;
  }

  /* Add particle size... */
  if(ctl->retrr) {
    if(x!=NULL)
      gsl_vector_set(x, n, aero->rr[0]);
    if(iqa!=NULL)
      iqa[n]=IDXRR;
    if(ipa!=NULL)
      ipa[n]=-1;
    (n)++;
  }

  /* Add particle size distribution width... */
  if(ctl->retss) {
    if(x!=NULL)
      gsl_vector_set(x, n, aero->ss[0]);
    if(iqa!=NULL)
      iqa[n]=IDXSS;
    if(ipa!=NULL)
      ipa[n]=-1;
    (n)++;
  }
  
  return n;
}

/*****************************************************************************/

void atm2x_help(atm_t *atm,
		double zmin,
		double zmax,
		double *value,
		int val_iqa,
		gsl_vector *x,
		int *iqa,
		int *ipa,
		size_t *n) {
  
  int ip;

  /* Add elements to state vector... */
  for(ip=0; ip<atm->np; ip++)
    if(atm->z[ip]>=zmin && atm->z[ip]<=zmax) {
      if(x!=NULL)
	gsl_vector_set(x, *n, value[ip]);
      if(iqa!=NULL)
	iqa[*n]=val_iqa;
      if(ipa!=NULL)
	ipa[*n]=ip;
      (*n)++;
    }
}

/*****************************************************************************/

void write_atm(const char *dirname,
	       const char *filename,
	       ctl_t *ctl,
	       atm_t *atm) {
  
  FILE *out;
  
  char file[LEN];
  
  int ig, ip, iw, n=6;
  
  /* Set filename... */
  if(dirname!=NULL)
    sprintf(file, "%s/%s", dirname, filename);
  else
    sprintf(file, "%s", filename);
  
  /* Write info... */
  printf("Write atmospheric data: %s\n", file);
  
  /* Create file... */
  if(!(out=fopen(file, "w")))
    ERRMSG("Cannot create file!");
  
  /* Write header... */
  fprintf(out,
	  "# $1 = time (seconds since 2000-01-01T00:00Z)\n"
	  "# $2 = altitude [km]\n"
	  "# $3 = longitude [deg]\n"
	  "# $4 = latitude [deg]\n"
	  "# $5 = pressure [hPa]\n"
	  "# $6 = temperature [K]\n");
  for(ig=0; ig<ctl->ng; ig++)
    fprintf(out, "# $%d = %s volume mixing ratio\n",
	    ++n, ctl->emitter[ig]);
  for(iw=0; iw<ctl->nw; iw++)
    fprintf(out, "# $%d = window %d: extinction [1/km]\n", ++n, iw);
  
  /* Write data... */
  for(ip=0; ip<atm->np; ip++) {
    if(ip==0 || atm->lat[ip]!=atm->lat[ip-1] || atm->lon[ip]!=atm->lon[ip-1])
      fprintf(out, "\n");
    fprintf(out, "%.2f %g %g %g %g %g", atm->time[ip], atm->z[ip],
	    atm->lon[ip], atm->lat[ip], atm->p[ip], atm->t[ip]);
    for(ig=0; ig<ctl->ng; ig++)
      fprintf(out, " %g", atm->q[ig][ip]);
    for(iw=0; iw<ctl->nw; iw++)
      fprintf(out, " %g", atm->k[iw][ip]);
    fprintf(out, "\n");
  }
  
  /* Close file... */
  fclose(out);
}

/*****************************************************************************/
