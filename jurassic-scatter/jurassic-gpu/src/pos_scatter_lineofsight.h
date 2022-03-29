#pragma once

#include "jr_common.h"
#include "assert.h"

/* The functions of this file have been copied from lineofsight.h and 
 * lineofsight.c from jurassic-scatter project. Also, jur_ preffix
 * (and in some cases jur_sca_) is added to the functions to avoid
 * function name collisions.
 * 
 * All functions are host functions.
 */

int pos_scatter_traceray(ctl_t *ctl, 
    atm_t *atm, 
    obs_t *obs, 
    aero_t *aero,
    int const ir, 
    pos_t los[], 
    double *tsurf,
    int const ignore_scattering); 

/* Find ground or TOA intersection point of a LOS. */
void pos_scatter_jur_intersection_point(ctl_t *ctl,
    atm_t *atm,
    double *znew,
    pos_t los[],
    int ip,
    pos_t los_aero[],
    int jp);

/* Add points to LOS for fine sampling of aerosol/cloud layers. */
int pos_scatter_jur_add_aerosol_layers(ctl_t *ctl,
    atm_t *atm,
    pos_t los[],
    aero_t *aero,
    int np);

int pos_scatter_jur_add_aerosol_layers_with_los_aero(ctl_t *ctl,
    atm_t *atm,
    pos_t los[],
    aero_t *aero,
    int np);


int pos_scatter_traceray(ctl_t *ctl, atm_t *atm, obs_t *obs, aero_t *aero, int const ir, 
    pos_t los[], double *tsurf, int const ignore_scattering) {
  double ex0[3], ex1[3], q[NG], k[NW], lat, lon, p, t, x[3], xobs[3], xvp[3], z = 1e99, z_low=z, zmax, zmin, zrefrac = 60;
  // Initialize
  *tsurf = -999;
  for(int ig = 0; ig < NG; ig++) q[ig] = 0;
  for(int iw = 0; iw < NW; iw++) k[iw] = 0;
  obs->tpz[ir]   = obs->vpz[ir];
  obs->tplon[ir] = obs->vplon[ir];
  obs->tplat[ir] = obs->vplat[ir];
  size_t atmIdx=0; int atmNp=0;
  locate_atm(atm, obs->time[ir], &atmIdx, &atmNp);
  altitude_range_nn(atm, atmIdx, atmNp, &zmin, &zmax);
  if(obs->obsz[ir] < zmin)        return 0;																		// Check observer altitude
  if(obs->vpz[ir] > zmax - 0.001) return 0;																		// Check view point altitude
  jur_geo2cart(obs->obsz[ir], obs->obslon[ir], obs->obslat[ir], xobs);				// Cart. coordinates of observer
  jur_geo2cart(obs->vpz[ir],	obs->vplon[ir],  obs->vplat[ir],	xvp);					// and view point
  UNROLL
    for(int i = 0; i < 3; i++) ex0[i] = xvp[i] - xobs[i];											// Determine initial tangent vector
  double const norm = NORM(ex0);
  UNROLL
    for(int i = 0; i < 3; i++) {
      ex0[i] /= norm;
      x[i] = xobs[i];																													// Observer within atmosphere
    } // i
  if(obs->obsz[ir] > zmax) {																									// Above atmosphere, search entry point
    double dmax = norm, dmin = 0.;
    while(fabs(dmin - dmax) > 0.001) {
      double const d = 0.5*(dmax + dmin);
      UNROLL
        for(int i = 0; i < 3; i++) x[i] = xobs[i] + d*ex0[i];
      z = cart2alt(x); // no need to compute lat and lon here
      if((z <= zmax) && (z > zmax - 0.001)) break;
      if(z < zmax - 0.0005) dmax = d;
      else									dmin = d;
    } // while
  }

  int np = 0, z_low_idx=-1;
  for(int stop = 0; np < NLOS; ++np) {																				// Ray-tracing
    double ds = ctl->rayds, dz = ctl->raydz;																	// Set step length
    if(dz > 0.) {
      double const norm_x = 1.0/NORM(x);
      double dot = 0.;
      UNROLL
        for(int i = 0; i < 3; i++) {
          dot += ex0[i]*x[i]*norm_x;
        }
      double const cosa = fabs(dot);
      if(cosa != 0.) ds = fmin(ds, dz/cosa);
    }
    jur_cart2geo(x, &z, &lon, &lat);																					// Determine geolocation
    double EPS = 0.001;
    if(ignore_scattering) EPS = 0.0;
    if((z < zmin + EPS) || (z > zmax + EPS)) {																						// LOS escaped
      double xh[3];
      stop = (z < zmin + EPS) ? 2 : 1;
      jur_geo2cart(los[np - 1].z, los[np - 1].lon, los[np - 1].lat, xh);
      double const zfrac = (z < zmin + EPS) ? zmin : zmax;
      double const frac = (zfrac - los[np - 1].z)/(z - los[np - 1].z);
      UNROLL
        for(int i = 0; i < 3; i++) x[i] = xh[i] + frac*(x[i] - xh[i]);
      jur_cart2geo(x, &z, &lon, &lat);
      if(ignore_scattering) {
        los[np - 1].ds = ds*frac;
        ds = 0.;
      }
      else {
        intpol_atm_geo_pt(ctl, atm, (int) atmIdx, atmNp, z, lon, lat, &los[np].p, &los[np].t);		// Interpolate atmospheric data
        intpol_atm_geo_qk(ctl, atm, (int) atmIdx, atmNp, z, lon, lat, los[np].q, los[np].k);			// Interpolate atmospheric data
        los[np].z = z;
        los[np].lon = lon;
        los[np].lat = lat;
        los[np].ds=0.;
      }
    }
    if(ignore_scattering || (!ignore_scattering && stop == 0)) {
      intpol_atm_geo_pt(ctl, atm, (int) atmIdx, atmNp, z, lon, lat, &p, &t);		// Interpolate atmospheric data
      intpol_atm_geo_qk(ctl, atm, (int) atmIdx, atmNp, z, lon, lat, q, k);			// Interpolate atmospheric data
      write_pos_point(los + np, lon, lat, z, p, t, q, k, ds);
    }
#ifdef GPUDEBUG
      los[np].ir = ir; los[np].ip = np; // for DEBUGging
#endif
    if(z < z_low) {
      z_low = z;
      z_low_idx = np; // store the index of the point where the altitude is lowest, used to compute the tangent point later
    }

    if(stop) { *tsurf = (stop == 2 ? t : -999); break; }											// Hit ground or space?

    double n = 1., ng[] = {0., 0., 0.};
    if(ctl->refrac && z <= zrefrac) {																					// Compute gradient of refractivity
      n += jur_refractivity(p, t);
      double xh[3];
      UNROLL
        for(int i = 0; i < 3; i++) xh[i] = x[i] + 0.5*ds*ex0[i];
      jur_cart2geo(xh, &z, &lon, &lat);
      intpol_atm_geo_pt(ctl, atm, (int) atmIdx, atmNp, z, lon, lat, &p, &t);
      double const n2 = jur_refractivity(p, t);
      for(int i = 0; i < 3; i++) {
        double const h = 0.02;
        xh[i] += h;
        jur_cart2geo(xh, &z, &lon, &lat);
        intpol_atm_geo_pt(ctl, atm, (int) atmIdx, atmNp, z, lon, lat, &p, &t);
        ng[i] = (jur_refractivity(p, t) - n2)/h;
        xh[i] -= h;
      } // i
    }
    UNROLL
      for(int i = 0; i < 3; i++) ex1[i] = ex0[i]*n + ds*ng[i];								// Construct new tangent vector

    double const norm_ex1 = NORM(ex1);
    for(int i = 0; i < 3; i++) {
      ex1[i] /= norm_ex1;
      x[i] += 0.5*ds*(ex0[i] + ex1[i]);																				// Determine next point of LOS
      ex0[i] = ex1[i];																												// Copy tangent vector
    } // i
  } // np
  ++np;
#ifndef __NVCC__
  if (NLOS <= np) ERRMSG("Too many LOS points!");
#endif

  // FIXME: added..
  if(!ignore_scattering) {
    /* Check length of last segment... */
    if(los[np - 2].ds < 1e-3 && np - 1 > 1)
      np--;
  }

  // Get tangent point (before changing segment lengths!)
  jur_tangent_point(los, np, z_low_idx, &obs->tpz[ir], &obs->tplon[ir], &obs->tplat[ir]);
  trapezoid_rule_pos(np, los);
  column_density(ctl->ng, los, np);
#ifdef CURTIS_GODSON
  if(ctl->formod == 1) curtis_godson(ctl, los, np);
  // this could be done during the while loop using the aux-variables:
  //					double cgpxu[NG];	/*! Curtis-Godson pressure times column density */
  //					double cgtxu[NG];	/*! Curtis-Godson temperature times column density */
#else
  assert(1 != ctl->formod);
#endif

  if(!ignore_scattering) {
    if (ctl->sca_n > 0) {
      np = pos_scatter_jur_add_aerosol_layers(ctl, atm, los, aero, np);
    }
  }
  return np;
} // traceray

/*****************************************************************************/


// FIXME: these 2 functions are added to remove gsl min and 
//        gsl_stats_min and gsl_stats_max functions
double min_value_in_array(double *arr, int n) {
  double ret = arr[0];
  for(int i = 1; i < n; i++)
    if(arr[i] < ret)
      ret = arr[i];
  return ret;
}

double max_value_in_array(double *arr, int n) {
  double ret = arr[0];
  for(int i = 1; i < n; i++)
    if(arr[i] > ret)
      ret = arr[i];
  return ret;
}

// comparator for qsort (reversed)
int cmp(const void *a, const void *b) {
   if (*(double*) a > *(double*) b) return -1;
   if (*(double*) a < *(double*) b) return 1;
   return 0;
} 

// copy pos
void copy_pos(const pos_t *source, pos_t *dest) {
  dest -> z   = source -> z;
  dest -> lat = source -> lat;
  dest -> lon = source -> lon;
  dest -> ds  = source -> ds;  
  dest -> t   = source -> t;
  dest -> p   = source -> p;
  for(int ig = 0; ig < NG; ig++)
    dest -> q[ig] = source -> q[ig];
  for(int iw = 0; iw < NW; iw++)
    dest -> k[iw] = source -> k[iw];
  
}

int pos_scatter_jur_add_aerosol_layers(ctl_t *ctl,
    atm_t *atm,
    pos_t los[],
    aero_t *aero,
    int np) {

  double alti[8*NLMAX], altimax, altimin, x1[3], x2[3], x3[3], tt=0., epsilon=0.005; 
  /* deltatop=10., deltabot=10., */

  int il, ig, jl=0, ip, it;

  /* Create altitudes to sample aerosol edges */
  for (il=0; il<aero->nl;il++){
    alti[jl] = aero->top[il] + epsilon;
    alti[jl+1] = aero->top[il] - epsilon;
    alti[jl+2] = aero->bottom[il] + epsilon;
    alti[jl+3] = aero->bottom[il] - epsilon;
    jl = jl+4;

    /* Create altitudes to sample transition layers */
    if (aero->trans[il] > ctl->transs) {
      tt = aero->trans[il] / ctl->transs;

      alti[jl] = aero->top[il] + aero->trans[il] + epsilon;
      alti[jl+1] = aero->top[il] + aero->trans[il] - epsilon;
      alti[jl+2] = aero->bottom[il] - aero->trans[il] + epsilon;
      alti[jl+3] = aero->bottom[il] - aero->trans[il] - epsilon;
      jl = jl+4;
      for (it=1; it<(int)tt; it++){
        alti[jl] = aero->top[il] + aero->trans[il] - epsilon - it*ctl->transs;
        jl++;
        alti[jl] = aero->bottom[il] - aero->trans[il] + epsilon + it*ctl->transs;
        jl++;
      }
    }  
  }

  if(jl >= 8 * NLMAX) ERRMSG("You should increase NLMAX!");

  qsort(alti, jl, sizeof(double), cmp); 

  altimax = max_value_in_array(alti, jl);
  altimin = min_value_in_array(alti, jl);

  int los_aero_np = 0;
  for(int i = 0; i < np; i++) {
    copy_pos(&los[i], &los[NLOS - np + i]);
  }
  copy_pos(&los[NLOS - np + 0], &los[0]);
  los_aero_np = 1;

  int prev_pos_index = 0;

  for (ip=1; ip<np; ip++){

    /* add new los points around cloud edges */
    if ( (los[prev_pos_index].z < altimax || los[NLOS - np + ip].z < altimax) &&
        (los[prev_pos_index].z > altimin || los[NLOS - np + ip].z > altimin) ) { 
      for (il=0; il<jl;il++){ /* loop over cloud edges */
        /* von oben */
        if(los[prev_pos_index].z > alti[il] && los[NLOS - np + ip].z < alti[il]){
          if(los_aero_np >= NLOS - np + ip) ERRMSG("Too many LOS points!");
          pos_scatter_jur_intersection_point(ctl, atm, &alti[il], los, NLOS - np + ip, los, los_aero_np);
          los_aero_np++; 
        }
        /* von unten */
        if(los[prev_pos_index].z < alti[jl-il-1] && los[NLOS - np + ip].z > alti[jl-il-1]){
          if(los_aero_np >= NLOS - np + ip) ERRMSG("Too many LOS points!");
          pos_scatter_jur_intersection_point(ctl, atm, &alti[jl-il-1], los, NLOS - np + ip, los, los_aero_np);
          los_aero_np++;
        }
      }
    }
    copy_pos(&los[NLOS - np + ip], &los[los_aero_np]);
    prev_pos_index = los_aero_np;

    if(los_aero_np >= NLOS - np + ip)
      ERRMSG("Too many LOS points!");

    /* Increment and check number of new LOS points */
    if((los_aero_np++)>NLOS)
      ERRMSG("Too many LOS points!");
  }

  /* Compute segment length following trapezoidal rule */
  jur_geo2cart(los[0].z, los[0].lon,los[0].lat, x1);
  jur_geo2cart(los[1].z, los[1].lon,los[1].lat, x2);
  los[0].ds= 0.5 * (DIST(x1,x2));
  for(ip=1; ip<los_aero_np-1; ip++){
    jur_geo2cart(los[ip-1].z, los[ip-1].lon, los[ip-1].lat, x1);
    jur_geo2cart(los[ip].z, los[ip].lon, los[ip].lat, x2);
    jur_geo2cart(los[ip+1].z, los[ip+1].lon, los[ip+1].lat, x3);
    los[ip].ds = 0.5 * (DIST(x1,x2) + DIST(x2,x3));
  }
  jur_geo2cart(los[los_aero_np-1].z, los[los_aero_np-1].lon, 
      los[los_aero_np-1].lat, x1);
  jur_geo2cart(los[los_aero_np-2].z, los[los_aero_np-2].lon,
      los[los_aero_np-2].lat, x2);
  los[los_aero_np-1].ds = 0.5 * (DIST(x1,x2));

  /* add aerosol/cloud information and column density u to new los  */
  for (ip=0; ip<los_aero_np; ip++){

    /* Compute column density... */
    for(ig=0; ig<ctl->ng; ig++)
      los[ip].u[ig]=10*los[ip].q[ig]*los[ip].p
        /(GSL_CONST_MKSA_BOLTZMANN*los[ip].t)*los[ip].ds;

    /* Get aerosol/cloud layer id and factor */
    los[ip].aeroi = -999;
    los[ip].aerofac = 0.;
    // FIXME: added ip > 0 2 times..
    if ( ((ip > 0 && los[ip-1].z < altimax) || los[ip].z < altimax) &&
        ((ip > 0 && los[ip-1].z > altimin) || los[ip].z > altimin) ) { 
      for (il=0; il<aero->nl;il++){
        /* Aerosol info within layer centre */
        if (los[ip].z <= aero->top[il] && 
            los[ip].z >= aero->bottom[il]){
          los[ip].aeroi = il;
          los[ip].aerofac = 1.;	  
        }
        /* Aerosol info in transition region */
        if (aero->trans[il] > ctl->transs &&
            los[ip].z <= (aero->top[il] + aero->trans[il]) &&
            los[ip].z > aero->top[il]){
          los[ip].aeroi = il;
          los[ip].aerofac = (aero->top[il] + aero->trans[il] - los[ip].z)/
            aero->trans[il];
        }
        if (aero->trans[il] > ctl->transs &&
            los[ip].z < aero->bottom[il] &&
            los[ip].z >= (aero->bottom[il] - aero->trans[il])){
          los[ip].aeroi = il;
          los[ip].aerofac = fabs(aero->bottom[il]-aero->trans[il]-los[ip].z)/
            aero->trans[il];
        }
      }
    }
  }

  np = los_aero_np;
  return np;
}

int pos_scatter_jur_add_aerosol_layers_with_los_aero(ctl_t *ctl,
    atm_t *atm,
    pos_t los[],
    aero_t *aero,
    int np) {

  /* Allocate extended los... */
  pos_t *los_aero = (pos_t*) malloc((NLOS) * sizeof(pos_t));
  int los_aero_np = 0;

  double alti[8*NLMAX], altimax, altimin, x1[3], x2[3], x3[3], tt=0., epsilon=0.005; 
  /* deltatop=10., deltabot=10., */

  int il, ig, iw, jl=0, ip, it;

  size_t s;

  /* Create altitudes to sample aerosol edges */
  for (il=0; il<aero->nl;il++){
    alti[jl] = aero->top[il] + epsilon;
    alti[jl+1] = aero->top[il] - epsilon;
    alti[jl+2] = aero->bottom[il] + epsilon;
    alti[jl+3] = aero->bottom[il] - epsilon;
    jl = jl+4;

    /* Create altitudes to sample transition layers */
    if (aero->trans[il] > ctl->transs) {
      tt = aero->trans[il] / ctl->transs;

      alti[jl] = aero->top[il] + aero->trans[il] + epsilon;
      alti[jl+1] = aero->top[il] + aero->trans[il] - epsilon;
      alti[jl+2] = aero->bottom[il] - aero->trans[il] + epsilon;
      alti[jl+3] = aero->bottom[il] - aero->trans[il] - epsilon;
      jl = jl+4;
      for (it=1; it<(int)tt; it++){
        alti[jl] = aero->top[il] + aero->trans[il] - epsilon - it*ctl->transs;
        jl++;
        alti[jl] = aero->bottom[il] - aero->trans[il] + epsilon + it*ctl->transs;
        jl++;
      }
    }  
  }

  if(jl >= 8 * NLMAX) ERRMSG("You should increase NLMAX!");

  /* Sort all altitudes from top-down */
  for (il=0; il<jl;il++)
    alti[il]=alti[il]*(-1.);
  gsl_sort(alti,1,(size_t)jl);
  for (il=0; il<jl;il++)
    alti[il]=alti[il]*(-1.);

  altimax = gsl_stats_max(alti, 1, (size_t)jl);
  altimin = gsl_stats_min(alti, 1, (size_t)jl);

  /* Copy los to new los and add additional points */
  los_aero[0].z   = los[0].z;
  los_aero[0].lat = los[0].lat;
  los_aero[0].lon = los[0].lon;
  los_aero[0].ds  = los[0].ds;  
  los_aero[0].t   = los[0].t;
  los_aero[0].p   = los[0].p;
  for(ig=0; ig<ctl->ng; ig++)
    los_aero[0].q[ig] = los[0].q[ig];
  for(iw=0; iw<ctl->nw; iw++)
    los_aero[0].k[iw] = los[0].k[iw];
  los_aero_np = 1;

  for (ip=1; ip<np; ip++){

    /* add new los points around cloud edges */
    if ( (los[ip-1].z < altimax || los[ip].z < altimax) &&
        (los[ip-1].z > altimin || los[ip].z > altimin) ) { 
      for (il=0; il<jl;il++){ /* loop over cloud edges */
        /* von oben */
        if(los[ip-1].z > alti[il] && los[ip].z < alti[il]){
          pos_scatter_jur_intersection_point(ctl, atm, &alti[il], los, ip, los_aero, los_aero_np);
          los_aero_np++; 
        }
        /* von unten */
        if(los[ip-1].z < alti[jl-il-1] && los[ip].z > alti[jl-il-1]){
          pos_scatter_jur_intersection_point(ctl, atm, &alti[jl-il-1], los, ip, los_aero, los_aero_np);
          los_aero_np++;
        }
      }
    }

    /* /\* check if current altitude is closer than 2m *\/ */
    /* /\* to any cloud top or bottom *\/ */
    /* deltatop = 10; */
    /* deltabot = 10; */
    /* for (il=0; il<aero->nl;il++){ */
    /*   deltatop = fabs(los[ip].z - aero->top[il]); */
    /*   deltabot = fabs(los[ip].z - aero->bottom[il]); */
    /*   if ( deltatop < epsilon*2. || deltabot < epsilon*2. ){ */
    /* 	continue; */
    /*   } */
    /* } */

    /* only copy old los points, if they are outside top||bottom +-2m */ 
    /* if ( deltatop >= epsilon*2. && deltabot >= epsilon*2. ) {  */
    /* copy old los points */
    los_aero[los_aero_np].z = los[ip].z;
    los_aero[los_aero_np].lat = los[ip].lat;
    los_aero[los_aero_np].lon = los[ip].lon;
    los_aero[los_aero_np].t = los[ip].t;
    los_aero[los_aero_np].p = los[ip].p;
    for(ig=0; ig<ctl->ng; ig++)
      los_aero[los_aero_np].q[ig] = los[ip].q[ig];
    for(iw=0; iw<ctl->nw; iw++)
      los_aero[los_aero_np].k[iw] = los[ip].k[iw];

    /* Increment and check number of new LOS points */
    if((los_aero_np++)>NLOS)
      ERRMSG("Too many LOS points!");
    /* } */
  }

  /* Compute segment length following trapezoidal rule */
  jur_geo2cart(los_aero[0].z, los_aero[0].lon,los_aero[0].lat, x1);
  jur_geo2cart(los_aero[1].z, los_aero[1].lon,los_aero[1].lat, x2);
  los_aero[0].ds= 0.5 * (DIST(x1,x2));
  for(ip=1; ip<los_aero_np-1; ip++){
    jur_geo2cart(los_aero[ip-1].z, los_aero[ip-1].lon, los_aero[ip-1].lat, x1);
    jur_geo2cart(los_aero[ip].z, los_aero[ip].lon, los_aero[ip].lat, x2);
    jur_geo2cart(los_aero[ip+1].z, los_aero[ip+1].lon, los_aero[ip+1].lat, x3);
    los_aero[ip].ds = 0.5 * (DIST(x1,x2) + DIST(x2,x3));
  }
  jur_geo2cart(los_aero[los_aero_np-1].z, los_aero[los_aero_np-1].lon, 
      los_aero[los_aero_np-1].lat, x1);
  jur_geo2cart(los_aero[los_aero_np-2].z, los_aero[los_aero_np-2].lon,
      los_aero[los_aero_np-2].lat, x2);
  los_aero[los_aero_np-1].ds = 0.5 * (DIST(x1,x2));

  /* add aerosol/cloud information and column density u to new los  */
  for (ip=0; ip<los_aero_np; ip++){

    /* Compute column density... */
    for(ig=0; ig<ctl->ng; ig++)
      los_aero[ip].u[ig]=10*los_aero[ip].q[ig]*los_aero[ip].p
        /(GSL_CONST_MKSA_BOLTZMANN*los_aero[ip].t)*los_aero[ip].ds;

    /* Get aerosol/cloud layer id and factor */
    los_aero[ip].aeroi = -999;
    los_aero[ip].aerofac = 0.;
    // FIXME: added ip > 0 2 times..
    if ( ((ip > 0 && los_aero[ip-1].z < altimax) || los_aero[ip].z < altimax) &&
        ((ip > 0 && los_aero[ip-1].z > altimin) || los_aero[ip].z > altimin) ) { 
      for (il=0; il<aero->nl;il++){
        /* Aerosol info within layer centre */
        if (los_aero[ip].z <= aero->top[il] && 
            los_aero[ip].z >= aero->bottom[il]){
          los_aero[ip].aeroi = il;
          los_aero[ip].aerofac = 1.;	  
        }
        /* Aerosol info in transition region */
        if (aero->trans[il] > ctl->transs &&
            los_aero[ip].z <= (aero->top[il] + aero->trans[il]) &&
            los_aero[ip].z > aero->top[il]){
          los_aero[ip].aeroi = il;
          los_aero[ip].aerofac = (aero->top[il] + aero->trans[il] - los_aero[ip].z)/
            aero->trans[il];
        }
        if (aero->trans[il] > ctl->transs &&
            los_aero[ip].z < aero->bottom[il] &&
            los_aero[ip].z >= (aero->bottom[il] - aero->trans[il])){
          los_aero[ip].aeroi = il;
          los_aero[ip].aerofac = fabs(aero->bottom[il]-aero->trans[il]-los_aero[ip].z)/
            aero->trans[il];
        }
      }
    }
  }

  /* Copy los */
  /*   *los = *los_aero; */
  s=(NLOS) * sizeof(pos_t);
  memcpy(los, los_aero, s);

  np = los_aero_np;

  /* Free help los... */
  free(los_aero);

  return np;
}

/*****************************************************************************/

void pos_scatter_jur_intersection_point(ctl_t *ctl,
    atm_t *atm,
    double *znew,
    pos_t los[],
    int ip,
    pos_t los_aero[],
    int jp){

  double frac, x1[3], x2[3];
  int i ;

  frac = (los[ip-1].z - *znew) / (los[ip-1].z - los[ip].z);
  jur_geo2cart(los[ip-1].z, los[ip-1].lon, los[ip-1].lat, x1);
  jur_geo2cart(los[ip].z, los[ip].lon, los[ip].lat, x2);

  for(i=0; i<3; i++)
    x2[i]=x1[i]+frac*(x2[i]-x1[i]);
  /* get new coordinates */
  jur_cart2geo(x2, &los_aero[jp].z, &los_aero[jp].lon, &los_aero[jp].lat);
  /* get atmosphere parameters */
  jur_intpol_atm_geo(ctl, atm, los_aero[jp].z, los_aero[jp].lon, los_aero[jp].lat,
      &los_aero[jp].p, &los_aero[jp].t, los_aero[jp].q, los_aero[jp].k);
}
