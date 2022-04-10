#ifndef _INTERFACE_H_
#define _INTERFACE_H_

#include "jr_scatter_gpu.h"

void formod_multiple_packages(ctl_t *ctl, atm_t *atm, aero_t *aero, int n, obs_t *obs_packages); 
void initialize_jurassic_gpu_table(ctl_t *ctl);

int call_raytrace(ctl_t *ctl, atm_t *atm, obs_t *obs, aero_t *aero, int const ir, 
                            pos_t los[], double *tsurf, int const ignore_scattering);

int raytrace_from_jr_common(ctl_t *ctl, atm_t *atm, obs_t *obs, aero_t *aero, int const ir, 
                            pos_t los[], double *tsurf, int const ignore_scattering);

trans_table_t* call_get_tbl(ctl_t const *ctl);
trans_table_t* get_tbl_from_jr_common(ctl_t const *ctl);

double call_src_planck_core(trans_table_t const *tbl, double const t, int const id); 
double src_planck_core_from_jr_common(trans_table_t const *tbl, double const t, int const id); 


double continua_core_CPU(ctl_t const *ctl, pos_t const *los, int const id);
double continua_core_CPU_from_CPUdrivers(ctl_t const *ctl, pos_t const *los, int const id);

#endif
