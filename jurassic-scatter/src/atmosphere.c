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

void copy_atm(ctl_t *ctl,
	      atm_t *atm_dest,
	      atm_t *atm_src,
	      int init) {

  int ig, ip, iw;
  
  size_t s;

  /* Data size... */
  s=(size_t)atm_src->np*sizeof(double);
  
  /* Copy data... */
  atm_dest->np=atm_src->np;
  memcpy(atm_dest->time, atm_src->time, s);
  memcpy(atm_dest->z, atm_src->z, s);
  memcpy(atm_dest->lon, atm_src->lon, s);
  memcpy(atm_dest->lat, atm_src->lat, s);
  memcpy(atm_dest->p, atm_src->p, s);
  memcpy(atm_dest->t, atm_src->t, s);
  for(ig=0; ig<ctl->ng; ig++)
    memcpy(atm_dest->q[ig], atm_src->q[ig], s);
  for(iw=0; iw<ctl->nw; iw++)
    memcpy(atm_dest->k[iw], atm_src->k[iw], s);
  atm_dest->init=atm_src->init;
  
  /* Initialize... */
  if(init)
    for(ip=0; ip<atm_dest->np; ip++) {
      atm_dest->p[ip]=0;
      atm_dest->t[ip]=0;
      for(ig=0; ig<ctl->ng; ig++)
        atm_dest->q[ig][ip]=0;
      for(iw=0; iw<ctl->nw; iw++)
	atm_dest->k[iw][ip]=0;
    }
}

/*****************************************************************************/

int find_emitter(ctl_t *ctl,
		 const char *emitter) {
  
  int ig;

  for(ig=0; ig<ctl->ng; ig++)
    if(strcasecmp(ctl->emitter[ig], emitter)==0)
      return ig;
  
  return -1;
}

/*****************************************************************************/

double gravity(double z, 
	       double lat) {
  
  /* Compute gravity according to 1967 Geodetic Reference System... */
  return 9.780318*(1+0.0053024*gsl_pow_2(sin(lat/180*M_PI))
		   -0.0000058*gsl_pow_2(sin(2*lat/180*M_PI)))-3.086e-3*z;
}

/*****************************************************************************/

void intpol_atm(ctl_t *ctl,
		atm_t *atm_dest,
		atm_t *atm_src) {
  
  double k[NWMAX], q[NGMAX];
  
  int ig, ip, iw;
  
  /* Interpolate atmospheric data... */
  for(ip=0; ip<atm_dest->np; ip++) {
    intpol_atm_geo(ctl, atm_src, atm_dest->z[ip], atm_dest->lon[ip],
		   atm_dest->lat[ip], &atm_dest->p[ip], &atm_dest->t[ip],
		   q, k);
    for(ig=0; ig<ctl->ng; ig++)
      atm_dest->q[ig][ip]=q[ig];
    for(iw=0; iw<ctl->nw; iw++)
      atm_dest->k[iw][ip]=k[iw];
  }
}

/*****************************************************************************/

void intpol_atm_geo(ctl_t *ctl,
		    atm_t *atm,
		    double z0,
		    double lon0,
		    double lat0,
		    double *p,
		    double *t,
		    double *q,
		    double *k) {
  
  /* 1D interpolation (vertical profile)... */
  if(ctl->ip==1)
    intpol_atm_1d(ctl, atm, 0, atm->np, z0, p, t, q, k);
  
  /* 2D interpolation (satellite track)... */
  else if(ctl->ip==2)
    intpol_atm_2d(ctl, atm, z0, lon0, lat0, p, t, q, k);
  
  /* 3D interpolation (Lagrangian grid)... */
  else if(ctl->ip==3)
    intpol_atm_3d(ctl, atm, z0, lon0, lat0, p, t, q, k);
  
  /* Wrong parameter... */
  else
    ERRMSG("Unknown interpolation method, check IP!");
}

/*****************************************************************************/

void intpol_atm_1d(ctl_t *ctl,
		   atm_t *atm,
		   int idx0,
		   int n,
		   double z0,
		   double *p,
		   double *t,
		   double *q,
		   double *k) {
  
  int ig, ip, iw;
  
  /* Get array index... */
  ip=idx0+jur_locate(&atm->z[idx0], n, z0);
  
  /* Interpolate... */
  *p=EXP(atm->z[ip], atm->p[ip], atm->z[ip+1], atm->p[ip+1], z0);
  *t=LIN(atm->z[ip], atm->t[ip], atm->z[ip+1], atm->t[ip+1], z0);
  for(ig=0; ig<ctl->ng; ig++)
    q[ig]=LIN(atm->z[ip], atm->q[ig][ip], atm->z[ip+1], atm->q[ig][ip+1], z0);
  for(iw=0; iw<ctl->nw; iw++)
    k[iw]=LIN(atm->z[ip], atm->k[iw][ip], atm->z[ip+1], atm->k[iw][ip+1], z0);
}

/*****************************************************************************/

void intpol_atm_2d(ctl_t *ctl,
		   atm_t *atm,
		   double z0,
		   double lon0,
		   double lat0,
		   double *p,
		   double *t,
		   double *q,
		   double *k) {
  
  static double x1[NPMAX][3];
  
  static int idx[NPMAX], nx, nz[NPMAX];
  
  double dh, dhmin0=1e99, dhmin1=1e99, k0[NWMAX], k1[NWMAX], lat1=-999,
    lon1=-999, p0, p1, q0[NGMAX], q1[NGMAX], r, r0, r1, t0, t1, x0[3], x, x2;
  
  int ig, ip, ix, iw, ix0=0, ix1=0;
  
  /* Initialize... */
  if(!atm->init) {
    atm->init=1;
    
    /* Determine grid dimensions... */
    nx=0;
    for(ip=0; ip<atm->np; ip++) {
      if(atm->lon[ip]!=lon1 || atm->lat[ip]!=lat1) {
	if((++nx)>NPMAX)
	  ERRMSG("Too many profiles!");
	nz[nx-1]=0;
	lon1=atm->lon[ip];
	lat1=atm->lat[ip];
	jur_geo2cart(0, lon1, lat1, x1[nx-1]);
	idx[nx-1]=ip;
      }
      nz[nx-1]++;
    }

    /* Check profiles... */
    for(ix=0; ix<nx; ix++)
      if(nz[ix]<=1)
	ERRMSG("Cannot identify profiles. Check ordering of data points!");
  }
  
  /* Get Cartesian coordinates... */
  jur_geo2cart(0, lon0, lat0, x0);
  
  /* Find next neighbours... */
  for(ix=0; ix<nx; ix++) {
    
    /* Get squared horizontal distance... */
    dh=DIST2(x0, x1[ix]);
    
    /* Find neighbours... */
    if(dh<=dhmin0) {
      dhmin1=dhmin0;
      ix1=ix0;
      dhmin0=dh;
      ix0=ix;
    } else if(dh<=dhmin1) {
      dhmin1=dh;
      ix1=ix;
    }
  }
  
  /* Interpolate vertically... */
  intpol_atm_1d(ctl, atm, idx[ix0], nz[ix0], z0, &p0, &t0, q0, k0);
  intpol_atm_1d(ctl, atm, idx[ix1], nz[ix1], z0, &p1, &t1, q1, k1);
  
  /* Interpolate horizontally... */
  x2=DIST2(x1[ix0], x1[ix1]);
  x=sqrt(x2);
  r0=(dhmin0-dhmin1+x2)/(2*x);
  r1=x-r0;
  if(r0<=0)
    r=0;
  else if(r1<=0)
    r=1;
  else
    r=r0/(r0+r1);
  
  *p=(1-r)*p0+r*p1;
  *t=(1-r)*t0+r*t1;
  for(ig=0; ig<ctl->ng; ig++)
    q[ig]=(1-r)*q0[ig]+r*q1[ig];
  for(iw=0; iw<ctl->nw; iw++)
    k[iw]=(1-r)*k0[iw]+r*k1[iw];
}

/*****************************************************************************/

void intpol_atm_3d(ctl_t *ctl,
		   atm_t *atm,
		   double z0,
		   double lon0,
		   double lat0,
		   double *p,
		   double *t,
		   double *q,
		   double *k) {
  
  static double x1[NPMAX][3];
  
  double dx2, dz, w, wsum, x0[3];
  
  int ig, ip, iw;
  
  /* Get Cartesian coordinates... */
  if(!atm->init) {
    atm->init=1;
    for(ip=0; ip<atm->np; ip++)
      jur_geo2cart(0, atm->lon[ip], atm->lat[ip], x1[ip]);
  }
  
  /* Initialize.. */
  *p=*t=wsum=0;
  for(ig=0; ig<ctl->ng; ig++)
    q[ig]=0;
  for(iw=0; iw<ctl->nw; iw++)
    k[iw]=0;
  
  /* Loop over grid points... */
  for(ip=0; ip<atm->np; ip++) {
    
    /* Get vertical distance... */
    dz=fabs(atm->z[ip]-z0);
    if(dz>ctl->cz)
      continue;
    
    /* Check latitude distance... */
    if(fabs(atm->lat[ip]-lat0)*111.13>ctl->cx)
      continue;
    
    /* Get horizontal distance... */
    jur_geo2cart(0, lon0, lat0, x0);
    dx2=DIST2(x0, x1[ip]);
    if(dx2>gsl_pow_2(ctl->cx))
      continue;
    
    /* Distance-based weighting... */
    w=(1-dz/ctl->cz)*(1-sqrt(dx2)/ctl->cx);
    
    /* Average data... */
    wsum+=w;
    *p+=w*atm->p[ip];
    *t+=w*atm->t[ip];
    for(ig=0; ig<ctl->ng; ig++)
      q[ig]+=w*atm->q[ig][ip];
    for(iw=0; iw<ctl->nw; iw++)
      k[iw]+=w*atm->k[iw][ip];
  }
  
  /* Normalize... */
  *p/=wsum;
  *t/=wsum;
  for(ig=0; ig<ctl->ng; ig++)
    q[ig]/=wsum;
  for(iw=0; iw<ctl->nw; iw++)
    k[iw]/=wsum;
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
