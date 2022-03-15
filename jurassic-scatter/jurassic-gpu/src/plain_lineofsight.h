#pragma once

#include "jr_common.h"

int plain_locate_irr(
    double *xx,
    int n,
    double x) {

  int ilo = 0;
  int ihi = n - 1;
  int i = (ihi + ilo) >> 1;

  if (xx[i] < xx[i + 1])
    while (ihi > ilo + 1) {
      i = (ihi + ilo) >> 1;
      if (xx[i] > x)
        ihi = i;
      else
        ilo = i;
    } else
      while (ihi > ilo + 1) {
        i = (ihi + ilo) >> 1;
        if (xx[i] <= x)
          ihi = i;
        else
          ilo = i;
      }

    return ilo;
}

void plain_intpol_atm(
    ctl_t * ctl,
    atm_t * atm,
    double z,
    double *p,
    double *t,
    double *q,
    double *k) {

  /* Get array index... */
  int ip = plain_locate_irr(atm->z, atm->np, z);

  /* Interpolate... */
  *p = EXP(atm->z[ip], atm->p[ip], atm->z[ip + 1], atm->p[ip + 1], z);
  *t = LIN(atm->z[ip], atm->t[ip], atm->z[ip + 1], atm->t[ip + 1], z);
  for (int ig = 0; ig < ctl->ng; ig++)
    q[ig] =
      LIN(atm->z[ip], atm->q[ig][ip], atm->z[ip + 1], atm->q[ig][ip + 1], z);
  for (int iw = 0; iw < ctl->nw; iw++)
    k[iw] =
      LIN(atm->z[ip], atm->k[iw][ip], atm->z[ip + 1], atm->k[iw][ip + 1], z);
}

void plain_raytrace(
    ctl_t * ctl,
    atm_t * atm,
    obs_t * obs,
    los_t * los,
    int ir) {

  const double h = 0.02, zrefrac = 60;

  double ds, ex0[3], ex1[3], k[NW], lat, lon, n, ng[3], norm,
         p, q[NG], t, x[3], xh[3], xobs[3], xvp[3], z = 1e99, zmax, zmin;

  int stop = 0;

  /* Initialize... */
  los->np = 0;
  los->tsurf = -999;
  obs->tpz[ir] = obs->vpz[ir];
  obs->tplon[ir] = obs->vplon[ir];
  obs->tplat[ir] = obs->vplat[ir];

  /* Get altitude range of atmospheric data... */
  gsl_stats_minmax(&zmin, &zmax, atm->z, 1, (size_t) atm->np);

  /*(if (ctl->nsf > 0) { // ignore this because ctl->nsf is not in jurassic-scatter-gpu ctl_t struct
    zmin = GSL_MAX(atm->sfz, zmin);
    if (atm->sfp > 0) {
    int ip = locate_irr(atm->p, atm->np, atm->sfp);
    double zip = LIN(log(atm->p[ip]), atm->z[ip],
    log(atm->p[ip + 1]), atm->z[ip + 1], log(atm->sfp));
    zmin = GSL_MAX(zip, zmin);
    }
    }*/

  /* Check observer altitude... */
  if (obs->obsz[ir] < zmin)
    ERRMSG("Observer below surface!");

  /* Check view point altitude... */
  if (obs->vpz[ir] > zmax)
    return;

  /* Determine Cartesian coordinates for observer and view point... */
  jur_geo2cart(obs->obsz[ir], obs->obslon[ir], obs->obslat[ir], xobs);
  jur_geo2cart(obs->vpz[ir], obs->vplon[ir], obs->vplat[ir], xvp);

  /* Determine initial tangent vector... */
  for (int i = 0; i < 3; i++)
    ex0[i] = xvp[i] - xobs[i];
  norm = NORM(ex0);
  for (int i = 0; i < 3; i++)
    ex0[i] /= norm;

  /* Observer within atmosphere... */
  for (int i = 0; i < 3; i++)
    x[i] = xobs[i];

  /* Observer above atmosphere (search entry point)... */
  if (obs->obsz[ir] > zmax) {
    double dmax = norm, dmin = 0;
    while (fabs(dmin - dmax) > 0.001) {
      double d = (dmax + dmin) / 2;
      for (int i = 0; i < 3; i++)
        x[i] = xobs[i] + d * ex0[i];
      jur_cart2geo(x, &z, &lon, &lat);
      if (z <= zmax && z > zmax - 0.001)
        break;
      if (z < zmax - 0.0005)
        dmax = d;
      else
        dmin = d;
    }
  }

  /* Ray-tracing... */
  while (1) {

    /* Set step length... */
    ds = ctl->rayds;
    if (ctl->raydz > 0) {
      norm = NORM(x);
      for (int i = 0; i < 3; i++)
        xh[i] = x[i] / norm;
      double cosa = fabs(DOTP(ex0, xh));
      if (cosa != 0)
        ds = GSL_MIN(ctl->rayds, ctl->raydz / cosa);
    }

    /* Determine geolocation... */
    jur_cart2geo(x, &z, &lon, &lat);

    /* Check if LOS hits the ground or has left atmosphere... */
    if (z < zmin || z > zmax) {
      stop = (z < zmin ? 2 : 1);
      double frac =
        ((z <
          zmin ? zmin : zmax) - los->z[los->np - 1]) / (z - los->z[los->np -
          1]);
      jur_geo2cart(los->z[los->np - 1], los->lon[los->np - 1],
          los->lat[los->np - 1], xh);
      for (int i = 0; i < 3; i++)
        x[i] = xh[i] + frac * (x[i] - xh[i]);
      jur_cart2geo(x, &z, &lon, &lat);
      los->ds[los->np - 1] = ds * frac;
      ds = 0;
    }

    /* Interpolate atmospheric data... */
    plain_intpol_atm(ctl, atm, z, &p, &t, q, k);

    /* Save data... */
    los->lon[los->np] = lon;
    los->lat[los->np] = lat;
    los->z[los->np] = z;
    los->p[los->np] = p;
    los->t[los->np] = t;
    for (int ig = 0; ig < ctl->ng; ig++)
      los->q[los->np][ig] = q[ig];
    for (int id = 0; id < ctl->nd; id++)
      los->k[los->np][id] = k[ctl->window[id]];
    los->ds[los->np] = ds;
    /* Add cloud extinction... */ // ignore this because ctl->ncl is not in jurassic-scatter-gpu ctl_t struct
    /*
       if (ctl->ncl > 0 && atm->cldz > 0) {
       double aux = exp(-0.5 * POW2((z - atm->clz) / atm->cldz));
       for (int id = 0; id < ctl->nd; id++) {
       int icl = locate_irr(ctl->clnu, ctl->ncl, ctl->nu[id]);
       los->k[los->np][id]
       += aux * LIN(ctl->clnu[icl], atm->clk[icl],
       ctl->clnu[icl + 1], atm->clk[icl + 1], ctl->nu[id]);
       }
       } */

    /* Increment and check number of LOS points... */
    if ((++los->np) > NLOS)
      ERRMSG("Too many LOS points!");

    /* Check stop flag... */
    if (stop) {

      /* Set surface temperature...         // ignore this because ctl->ncl is not in jurassic-scatter-gpu ctl_t struct
         if (ctl->nsf > 0 && atm->t_surf > 0)  // maybe just if statemet should be ignored
         t = atm->t_surf; */
      los->tsurf = (stop == 2 ? t : -999);

      /* Set surface emissivity... */       // ignore this because los->sfeps is not in jurassic-scatter-gpu ctl_t struct
      /*for (int id = 0; id < ctl->nd; id++) {
        los->sfeps[id] = 1.0;
        if (ctl->nsf > 0) {
        int isf = locate_irr(ctl->sfnu, ctl->nsf, ctl->nu[id]);
        los->sfeps[id] = LIN(ctl->sfnu[isf], atm->sfeps[isf],
        ctl->sfnu[isf + 1], atm->sfeps[isf + 1],
        ctl->nu[id]);
        }
        }*/

      /* Leave raytracer... */
      break;
    }

    /* Determine refractivity... */
    if (ctl->refrac && z <= zrefrac)
      n = 1 + jur_refractivity(p, t);
    else
      n = 1;

    /* Construct new tangent vector (first term)... */
    for (int i = 0; i < 3; i++)
      ex1[i] = ex0[i] * n;

    /* Compute gradient of refractivity... */
    if (ctl->refrac && z <= zrefrac) {
      for (int i = 0; i < 3; i++)
        xh[i] = x[i] + 0.5 * ds * ex0[i];
      jur_cart2geo(xh, &z, &lon, &lat);
      plain_intpol_atm(ctl, atm, z, &p, &t, q, k);
      n = jur_refractivity(p, t);
      for (int i = 0; i < 3; i++) {
        xh[i] += h;
        jur_cart2geo(xh, &z, &lon, &lat);
        plain_intpol_atm(ctl, atm, z, &p, &t, q, k);
        ng[i] = (jur_refractivity(p, t) - n) / h;
        xh[i] -= h;
      }
    } else
      for (int i = 0; i < 3; i++)
        ng[i] = 0;

    /* Construct new tangent vector (second term)... */
    for (int i = 0; i < 3; i++)
      ex1[i] += ds * ng[i];

    /* Normalize new tangent vector... */
    norm = NORM(ex1);
    for (int i = 0; i < 3; i++)
      ex1[i] /= norm;

    /* Determine next point of LOS... */
    for (int i = 0; i < 3; i++)
      x[i] += 0.5 * ds * (ex0[i] + ex1[i]);

    /* Copy tangent vector... */
    for (int i = 0; i < 3; i++)
      ex0[i] = ex1[i];
  }

  /* Get tangent point (to be done before changing segment lengths!)... */
  jur_sca_tangent_point(los, &obs->tpz[ir], &obs->tplon[ir], &obs->tplat[ir]);

  /* Change segment lengths according to trapezoid rule... */
  for (int ip = los->np - 1; ip >= 1; ip--)
    los->ds[ip] = 0.5 * (los->ds[ip - 1] + los->ds[ip]);
  los->ds[0] *= 0.5;

  /* Compute column density... */
  for (int ip = 0; ip < los->np; ip++)
    for (int ig = 0; ig < ctl->ng; ig++)
      los->u[ip][ig] = 10 * los->q[ip][ig] * los->p[ip]
        / (GSL_CONST_MKSA_BOLTZMANN * los->t[ip]) * los->ds[ip];
}

/*

// compile errors: 
// the differences between jurassic and jurassic-scatter-gpu struct variables

plain_lineofsight.h(21): error: class "los_t" has no member "sft"         -> tsurf
plain_lineofsight.h(28): error: class "ctl_t" has no member "nsf"         -> ignore (comment)
plain_lineofsight.h(29): error: class "atm_t" has no member "sfz"         -> ignore
plain_lineofsight.h(30): error: class "atm_t" has no member "sfp"         -> ignore
plain_lineofsight.h(31): error: identifier "locate_irr" is undefined      -> ignore
plain_lineofsight.h(47): error: identifier "geo2cart" is undefined        -> jur_geo2cart
plain_lineofsight.h(68): error: identifier "cart2geo" is undefined        -> jur_cart2geo
plain_lineofsight.h(112): error: identifier "intpol_atm" is undefined     -> jur_intpol_atm  <- copied plain_intpol_atm & plain_locate_irr from the plain  from the plain versionversion
plain_lineofsight.h(127): error: class "ctl_t" has no member "ncl"        -> ignore
plain_lineofsight.h(127): error: class "atm_t" has no member "cldz"       -> ignore
plain_lineofsight.h(128): error: class "atm_t" has no member "clz"        -> ignore
plain_lineofsight.h(128): error: identifier "POW2" is undefined           -> ignore
plain_lineofsight.h(130): error: class "ctl_t" has no member "clnu"       -> ignore
plain_lineofsight.h(132): error: class "atm_t" has no member "clk"        -> ignore
plain_lineofsight.h(145): error: class "ctl_t" has no member "nsf"        -> ignore, maybe just if statement should be ignored
plain_lineofsight.h(151): error: class "los_t" has no member "sfeps"      -> ignore
plain_lineofsight.h(153): error: class "ctl_t" has no member "sfnu"       -> ignore
plain_lineofsight.h(166): error: identifier "refractivity" is undefined   -> jur_refractivity
plain_lineofsight.h(211): error: identifier "tangent_point" is undefined  -> jur_sca_tangent_point
plain_lineofsight.h(222): error: identifier "KB" is undefined             -> GSL_CONST_MKSA_BOLTZMANN
*/
