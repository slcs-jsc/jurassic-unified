#ifndef _INTERFACE_H_
#define _INTERFACE_H_

#include "jr_scatter_gpu.h"

#define ptr const restrict

void jur_formod_multiple_packages(ctl_t const *ctl, atm_t *atm, obs_t *obs, int n, aero_t const *aero);

void jur_table_initialization(ctl_t *ctl);

int pos_scatter_traceray(ctl_t *ctl, atm_t *atm, obs_t *obs, aero_t *aero, int const ir,
                            pos_t los[], double *tsurf, int const ignore_scattering);

trans_table_t* get_tbl(ctl_t const *ctl);

double src_planck_core(trans_table_t const *tbl, double const t, int const id);

double continua_core_CPU(ctl_t const *ctl, pos_t const *los, int const id);

double apply_ega_core(trans_table_t const *tbl, pos_t const *los, double (*ptr tau_path), int const ng, int const id);

void jur_read_atm(char const *dirname, char const *filename, ctl_t *ctl, atm_t *atm);

double jur_scan_ctl(int argc, char *argv[], char const *varname, int arridx, char const *defvalue, char *value);

void jur_cart2geo(double const x[], double *alt, double *lon, double *lat);

void jur_copy_obs(ctl_t const *const ctl, obs_t *obs_dest, obs_t const *const obs_src, int const init);

void jur_geo2cart(double const alt, double const lon, double const lat, double x[]);

int jur_locate(double const *ptr xx, int const n, double const x);

void jur_read_obs(char const *dirname, char const *filename, ctl_t *ctl, obs_t *obs);

void jur_write_obs(char const *dirname, char const *filename, ctl_t *ctl, obs_t *obs);

void hydrostatic1d_CPU(ctl_t const *ctl, atm_t *atm, int const nr, int const ig_h2o);

int jur_find_emitter(ctl_t const *ctl, char const *emitter);

double brightness_core(double const rad, double const nu);

double jur_planck(double const t, double const nu);

#endif
