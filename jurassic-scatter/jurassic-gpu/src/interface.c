#include "inter_folder/interface.h"

#include "jurassic.h"
#include <stdio.h>

void formod_multiple_packages(ctl_t *ctl, atm_t *atm,  aero_t *aero, int n, obs_t *obs_packages) { 
  printf("DEBUG #%d number of packages.. %d\n", ctl->MPIglobrank, n);
  printf("DEBUG #%d call jur_formod..\n", ctl->MPIglobrank);
	jur_formod(ctl, atm, obs_packages, aero, n);
}

void initialize_jurassic_gpu_table(ctl_t *ctl) {
  printf("DEBUG #%d call table_initializaiton..\n", ctl->MPIglobrank);
  jur_table_initialization(ctl); 
}

int raytrace_from_jr_common(ctl_t *ctl, atm_t *atm, obs_t *obs, aero_t *aero, int const ir, 
                            pos_t los[], double *tsurf, int const ignore_scattering) {
  return call_traceray(ctl, atm, obs, aero, ir, los, tsurf, ignore_scattering);
}

trans_table_t* get_tbl_from_jr_common(ctl_t const *ctl) {
  return call_get_tbl(ctl);
}

double src_planck_core_from_jr_common(trans_table_t const *tbl, double const t, int const id) {
  return call_src_planck_core(tbl, t, id);
}

