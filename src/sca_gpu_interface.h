#ifndef _INTERFACE_H_
#define _INTERFACE_H_

#include "jurassic.h"

#define ptr const restrict

void jur_formod_multiple_packages(ctl_t const *ctl, atm_t *atm, obs_t *obs, int n, int32_t const *atm_id, aero_t const *aero);

void jur_table_initialization(ctl_t const *ctl);

int jur_traceray(ctl_t *ctl, atm_t *atm, obs_t *obs, int const ir,
    pos_t los[], double *tsurf, aero_t *aero, int scattering_included);

trans_table_t* jur_get_tbl(ctl_t const *ctl);

double jur_src_planck_core(trans_table_t const *tbl, double const t, int const id);

double jur_continua_core_CPU(ctl_t const *ctl, pos_t const *los, int const id);

double jur_apply_ega_core(trans_table_t const *tbl, pos_t const *los, double (*ptr tau_path), int const ng, int const id);

void jur_read_atm(char const *dirname, char const *filename, ctl_t *ctl, atm_t *atm);

double jur_scan_ctl(int argc, char *argv[], char const *varname, int arridx, char const *defvalue, char *value);

void jur_cart2geo(double const x[], double *alt, double *lon, double *lat);

void jur_copy_obs(ctl_t const *const ctl, obs_t *obs_dest, obs_t const *const obs_src, int const init);

void jur_geo2cart(double const alt, double const lon, double const lat, double x[]);

int jur_locate(double const *ptr xx, int const n, double const x);

void jur_read_obs(char const *dirname, char const *filename, ctl_t *ctl, obs_t *obs);

void jur_write_obs(char const *dirname, char const *filename, ctl_t *ctl, obs_t *obs);

void jur_hydrostatic1d_CPU(ctl_t const *ctl, atm_t *atm, int const nr, int const ig_h2o);

int jur_find_emitter(ctl_t const *ctl, char const *emitter);

double jur_brightness_core(double const rad, double const nu);

double jur_planck(double const t, double const nu);

void jur_copy_atm(ctl_t const *ctl, atm_t *atm_dest, atm_t const *atm_src, int const init);

void jur_read_ctl(int argc, char *argv[], ctl_t *ctl);

#endif
