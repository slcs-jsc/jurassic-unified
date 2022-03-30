#ifndef _INTERFACE_H_
#define _INTERFACE_H_

#include "jr_scatter_gpu.h"

void formod_multiple_packages(ctl_t *ctl, atm_t *atm, aero_t *aero, int n, obs_t *obs_packages); 
void initialize_jurassic_gpu_table(ctl_t *ctl);
int raytrace_from_jr_common(ctl_t *ctl, atm_t *atm, obs_t *obs, aero_t *aero, int const ir, 
                            pos_t los[], double *tsurf, int const ignore_scattering);

#endif
