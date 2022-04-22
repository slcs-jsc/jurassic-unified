#ifndef _INTERFACE_H_
#define _INTERFACE_H_

#include "jr_scatter_gpu.h"

#define ptr const restrict

void jur_formod(ctl_t const *ctl, atm_t *atm, obs_t *obs, aero_t const *aero, int n); 

void jur_table_initialization(ctl_t *ctl); 

int pos_scatter_traceray(ctl_t *ctl, atm_t *atm, obs_t *obs, aero_t *aero, int const ir, 
                            pos_t los[], double *tsurf, int const ignore_scattering);

trans_table_t* get_tbl(ctl_t const *ctl);

double src_planck_core(trans_table_t const *tbl, double const t, int const id); 

double continua_core_CPU(ctl_t const *ctl, pos_t const *los, int const id);

double apply_ega_core(trans_table_t const *tbl, pos_t const *los, double (*ptr tau_path), int const ng, int const id); 

void jur_read_atm(char const *dirname, char const *filename, ctl_t *ctl, atm_t *atm); 

double jur_scan_ctl(int argc, char *argv[], char const *varname, int arridx, char const *defvalue, char *value);

#endif
