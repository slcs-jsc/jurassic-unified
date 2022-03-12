#pragma once

#include "assert.h"

#include "jr_common.h" // contains jurassic-gpu raytracer 
#include "scatter_lineofsight.h" // containts jurassic-scatter raytracer
// #include "plain_lineofsight.h" // contains jurassic-plain raytracer

/*
// jurassic-scatter los type:
typedef struct { /// Line-of-sight data. ///////////////////////////////////////////////
  int np;                 /// Number of LOS points.
  double z[NLOS];         /// Altitude [km].
  double lon[NLOS];       /// Longitude [deg].
  double lat[NLOS];       /// Latitude [deg].
  double p[NLOS];         /// Pressure [hPa].
  double t[NLOS];         /// Temperature [K].
  double q[NLOS][NGMAX];  /// Volume mixing ratio.
  double k[NLOS][NWMAX];  /// Extinction [1/km].
  int aeroi [NLOS];       /// Aerosol/cloud layer index
  double aerofac[NLOS];   /// Aerosol/cloud layer scaling factor for transition layer
  double tsurf;           /// Surface temperature [K].
  double ds[NLOS];        /// Segment length [km].
  double u[NLOS][NGMAX];  /// Column density [molecules/cm^2].
} los_t; ///////////////////////////////////////////////////////////////////////////////

// jurassic-gpu los type: 
// FIXME: (acctually, this is jurassic-scatter-gpu, because it contains aeroi and aerofac)
typedef struct {    /// Point on the Line-of-sight data without storing //////////
  double z;		      /// Altitude [km].
  double lon;	      /// Longitude [deg].
  double lat;	      /// Latitude [deg].
  double p;		      /// Pressure [hPa].
  double t;		      /// Temperature [K].
  double q[NG];	    /// Volume mixing ratio.
  double k[NW];	    /// Extinction [1/km].
  int aeroi;        /// Aerosol/cloud layer index
  double aerofac;   /// Aerosol/cloud layer scaling factor for transition layer
  double ds;	      /// Segment length [km].
  double u[NG];	    /// Column density [molecules/cm^2].
#ifdef CURTIS_GODSON
  double cgp[NG];	  /// Curtis-Godson pressure [hPa].
  double cgt[NG];	  /// Curtis-Godson temperature [K].
  double cgu[NG];	  /// Curtis-Godson column density [molecules/cm^2].
#endif
#ifdef GPUDEBUG
  int ip, ir;       /// debug helpers
#endif
} pos_t; //////////////////////////////////////////////////////////////////////

// jurassic-plain los type:

typedef struct {
  int np;                 // Number of LOS points.
  double z[NLOS];         // Altitude [km].
  double lon[NLOS];       // Longitude [deg].
  double lat[NLOS];       // Latitude [deg].
  double p[NLOS];         // Pressure [hPa].
  double t[NLOS];         // Temperature [K].
  double q[NLOS][NG];     // Volume mixing ratio [ppv].
  double k[NLOS][ND];     // Extinction [1/km].
  double sft;             // Surface temperature [K].                     FIXME: "t_surf" in gpu and scatter versions
  double sfeps[NSF];      // Surface emissivity.                          FIXME: not in gpu and scatter versions
  double ds[NLOS];        // Segment length [km].
  double u[NLOS][NG];     // Column density [molecules/cm^2].
  double eps[NLOS][ND];   // Segment emissivity.                          FIXME: not in gpu and scatter versions
  double src[NLOS][ND];   // Segment source function [W/(m^2 sr cm^-1)].  FIXME: not in gpu and scatter versions
} los_t;
                                                                          FIXME: it doesn't contain aeroi
                                                                          FIXME: it doesn't contain aerofac
*/

// Additional functions for development and debugging
void convert_los_to_pos(pos_t *p, int *np, double *tsurf, los_t const *l) { 
  *np =     l -> np;
  *tsurf =  l -> tsurf;
  p = (pos_t*) malloc((*np) * sizeof(pos_t));
  for(int ip = 0; ip < l->np; ip++) {
    p->z       = l->z[ip];
    p->lon     = l->lon[ip];
    p->lat     = l->lat[ip];
    p->p       = l->p[ip];
    p->t       = l->t[ip];
    p->aerofac = l->aerofac[ip];
    p->ds      = l->ds[ip];
    for(int j = 0; j < NG; j++) {
      p->q[j] = l->q[ip][j];
      p->u[j] = l->u[ip][j];
    }
    for(int j = 0; j < NW; j++)
      p->k[j] = l->k[ip][j];
    // special case: not exactly the same in the both structures
    p->aeroi = (-999 == l->aeroi[ip] ? 0 : l->aeroi[ip]);
  }
}

void convert_pos_to_los(los_t *l, pos_t const *p, int const np, double const tsurf) {
  l -> np = np;
  l -> tsurf = tsurf;
  for(int ip = 0; ip < np; ip++) {
    l->z[ip]       = p[ip].z;
    l->lon[ip]     = p[ip].lon;
    l->lat[ip]     = p[ip].lat;
    l->p[ip]       = p[ip].p;
    l->t[ip]       = p[ip].t;
    l->aerofac[ip] = p[ip].aerofac;
    l->ds[ip]      = p[ip].ds;  
    for(int j = 0; j < NGMAX; j++) {
      l->q[ip][j] = p[ip].q[j];
      l->u[ip][j] = p[ip].u[j];
    }
    for(int j = 0; j < NWMAX; j++)
      l->k[ip][j] = p[ip].k[j];
    // special case: not exactly the same in the both structures
    l->aeroi[ip]   = (p[ip].aeroi == 0 ? -999 : p[ip].aeroi);
  }
}

void combined_raytrace(ctl_t *ctl,
    atm_t *atm,
    obs_t *obs,
    aero_t *aero,
    los_t *los,
    int ir) {

  pos_t *p = (pos_t*) malloc((NLOS) * sizeof(pos_t));
  double tsurf;
  int n = traceray(ctl, atm, obs, ir, p, &tsurf);

  los_t *l = (los_t*) malloc(sizeof(los_t));
  convert_pos_to_los(l, p, n, tsurf); 


  /* Add additional los points for aerosol layers and add aerosol data */
  if (ctl->sca_n > 0)
    jur_add_aerosol_layers(ctl,atm,l,aero);

  printf("num of points in ray: %d\n", n);

  memcpy(los, l, sizeof(los_t));

}

void update(int *ret_index, double *max_diff, int id, double diff) {
  if(diff < 0) diff *= -1;
  if(diff > *max_diff) {
    *max_diff = diff;
    *ret_index = id;
  }
}

int los_diff(los_t *a, los_t *b, ctl_t *ctl, int *ret_index, double *max_diff) {
  *ret_index = 0;
  *max_diff = 0;

  if(a->np != b->np)
    printf("DEBUG #%d number of los points differ: %d vs. %d\n", ctl->MPIglobrank, a->np, b->np);

  update(ret_index, max_diff, 1, a->np - b->np);

  for(int i = 0; i < a->np; i++) {
    update(ret_index, max_diff, 2, a->z[i] - b->z[i]);
    update(ret_index, max_diff, 3, a->lon[i] - b->lon[i]);
    update(ret_index, max_diff, 4, a->lat[i] - b->lat[i]);
    update(ret_index, max_diff, 5, a->p[i] - b->p[i]);
    update(ret_index, max_diff, 6, a->t[i] - b->t[i]);
    // TODO: YET NOT PART OF jurassic-gpu raytracer
    //update(ret_index, max_diff, 7, a->aerofac[i] - b->aerofac[i]);

    // TODO: Segment length [km]. WHY?
    // update(ret_index, max_diff, 8, a->ds[i] - b->ds[i]);

    for(int j = 0; j < ctl->ng; j++) {
      update(ret_index, max_diff, 9, a->q[i][j] - b->q[i][j]);
      // TODO: Column density [molecules/cm^2]. WHY? (also, very large values..)
      /* update(ret_index, max_diff, 10, a->u[i][j] - b->u[i][j]);
         if(a->u[i][j] - b->u[i][j] > 10000 || b->u[i][j] - a->u[i][j] > 10000)
         printf("DEBUG #%d %d-th column density: %lf vs. %lf\n", ctl->MPIglobrank, j, a->u[i][j], b->u[i][j]); */
    }

    for(int j = 0; j < ctl->nw; j++)
      update(ret_index, max_diff, 11, a->k[i][j] - b->k[i][j]);
  }

  // TODO: YET NOT PART OF jurassic-gpu raytracer
  /* for(int i = 0; i < a->np; i++) {
     update(ret_index, max_diff, 12, a->aeroi[i] - b->aeroi[i]);
     if(a->aeroi[i] - b->aeroi[i] == 999 || b->aeroi[i] - a->aeroi[i] == 999)
     printf("DEBUG #%d i-th point cloud layer index: %d vs. %d\n", ctl->MPIglobrank, i, a->aeroi[i], b->aeroi[i]);
     } */

  // TODO: WHY?
  // update(ret_index, max_diff, 13, a->tsurf - b->tsurf);
  return 0;
}

void compare_raytracers(ctl_t *ctl,
    atm_t *atm,
    obs_t *obs,
    aero_t *aero) {

  // compare computed lines of sight one by one
  for(int ir = 0; ir < obs-> nr; ir++) { 
    // jurassic-scatter (copied from scatter) raytracer
    // without calling jur_add_aerosol_layers() function
    los_t *jr_scatter_los = (los_t*) malloc(sizeof(los_t));
    jur_raytrace(ctl, atm, obs, aero, jr_scatter_los, ir, 1); // with ignoring scattering!

    // jurassic-gpu raytracer
    pos_t *p = (pos_t*) malloc((NLOS) * sizeof(pos_t));
    double tsurf;
    int n = traceray(ctl, atm, obs, ir, p, &tsurf);
    los_t *jr_gpu_los = (los_t*) malloc(sizeof(los_t));
    convert_pos_to_los(jr_gpu_los, p, n, tsurf); 

    // jurassic-plain raytracer
    // los_t *jr_plain_los = (los_t*) malloc(sizeof(los_t));
    // plain_raytrace(ctl, atm, obs, jr_plain_los, ir);

    // comparison of jurassic-scatter and jurassic-gpu raytracers
    int ret_index;
    double max_diff;
    int diff = los_diff(jr_scatter_los, jr_gpu_los, ctl, &ret_index, &max_diff);
    printf("DEBUG #%d ray #%d, max_diff = %lf, ret_index = %d\n", ctl->MPIglobrank, ir, max_diff, ret_index);
  }
}

