#include "jurassic_unified_library.h"
#include <assert.h>

// declarations of functions from jurassic-unified:
void jur_get_dimensions(int32_t *dimensions);
int jur_get_num_of_atms(int const nr, int32_t const *atm_id);
void jur_formod_multiple_packages(jur_ctl_t const *ctl, jur_atm_t *atm, jur_obs_t *obs, int n, int32_t const *atm_id, jur_aero_t const *aero); 
void jur_read_ctl(int argc, char *argv[], jur_ctl_t *ctl);
void jur_table_initialization(jur_ctl_t *ctl);

// declarations of functions from this file, which are not visible to Lars:
void convert_atm_from_reference_to_unified(atm_t const *atm, jur_atm_t *jur_atm); 
void convert_obs_from_reference_to_unified(obs_t const *obs, jur_obs_t *jur_obs);
void convert_obs_from_unified_to_reference(jur_obs_t const *jur_obs, obs_t *obs);

jur_ctl_t *jur_unified_init(int argc, char *argv[]) {
  static jur_ctl_t *jur_ctl = NULL;
  if(NULL == jur_ctl) {
    assert(argc && "Please call jur_unified_init(argc, argv) at the beginning of the program");
    int32_t dimensions[23];
    jur_get_dimensions(dimensions);
    assert(LENMAX     == dimensions[0]  && "jurassic_dimensions.h is inconsistent with the compiled libjurassic_unified.a");
    assert(MMAX       == dimensions[1]  && "jurassic_dimensions.h is inconsistent with the compiled libjurassic_unified.a");
    assert(NMAX       == dimensions[2]  && "jurassic_dimensions.h is inconsistent with the compiled libjurassic_unified.a");
    assert(NQMAX      == dimensions[3]  && "jurassic_dimensions.h is inconsistent with the compiled libjurassic_unified.a");
    assert(NWMAX      == dimensions[4]  && "jurassic_dimensions.h is inconsistent with the compiled libjurassic_unified.a");
    assert(NDMAX      == dimensions[5]  && "jurassic_dimensions.h is inconsistent with the compiled libjurassic_unified.a");
    assert(NGMAX      == dimensions[6]  && "jurassic_dimensions.h is inconsistent with the compiled libjurassic_unified.a");
    assert(NLOSMAX    == dimensions[7]  && "jurassic_dimensions.h is inconsistent with the compiled libjurassic_unified.a");
    assert(NPMAX      == dimensions[8]  && "jurassic_dimensions.h is inconsistent with the compiled libjurassic_unified.a");
    assert(NRMAX      == dimensions[9]  && "jurassic_dimensions.h is inconsistent with the compiled libjurassic_unified.a");
    assert(NSHAPEMAX  == dimensions[10] && "jurassic_dimensions.h is inconsistent with the compiled libjurassic_unified.a");
    assert(NFOVMAX    == dimensions[11] && "jurassic_dimensions.h is inconsistent with the compiled libjurassic_unified.a");
    assert(TBLNPMAX   == dimensions[12] && "jurassic_dimensions.h is inconsistent with the compiled libjurassic_unified.a");
    assert(TBLNSMAX   == dimensions[13] && "jurassic_dimensions.h is inconsistent with the compiled libjurassic_unified.a");
    assert(TBLNTMAX   == dimensions[14] && "jurassic_dimensions.h is inconsistent with the compiled libjurassic_unified.a");
    assert(TBLNUMAX   == dimensions[15] && "jurassic_dimensions.h is inconsistent with the compiled libjurassic_unified.a");
    assert(SCAMODMAX  == dimensions[16] && "jurassic_dimensions.h is inconsistent with the compiled libjurassic_unified.a");
    assert(NLMAX      == dimensions[17] && "jurassic_dimensions.h is inconsistent with the compiled libjurassic_unified.a");
    assert(NTHETAMAX  == dimensions[18] && "jurassic_dimensions.h is inconsistent with the compiled libjurassic_unified.a");
    assert(NRADMAX    == dimensions[19] && "jurassic_dimensions.h is inconsistent with the compiled libjurassic_unified.a");
    assert(REFMAX     == dimensions[20] && "jurassic_dimensions.h is inconsistent with the compiled libjurassic_unified.a");
    assert(RFMNPTSMAX == dimensions[21] && "jurassic_dimensions.h is inconsistent with the compiled libjurassic_unified.a");
    assert(RFMLINEMAX == dimensions[22] && "jurassic_dimensions.h is inconsistent with the compiled libjurassic_unified.a");

    jur_ctl = (jur_ctl_t*) malloc(sizeof(jur_ctl_t));
    jur_read_ctl(argc, argv, jur_ctl); 
    jur_table_initialization(jur_ctl); 
  }
  return jur_ctl;
}

void convert_atm_from_reference_to_unified(atm_t const *atm, jur_atm_t *jur_atm) {
  printf("[INFO] Warning: unified atm_t doesn't contain clz, cldz, clk, sfz, sfp, sft, ...\n");  
  jur_atm->np = atm->np;
  jur_atm->init = 0;
  for(int i = 0; i < atm->np; i++) {
    jur_atm->time[i] = atm->time[i];
    jur_atm->z[i]    = atm->z[i];
    jur_atm->lon[i]  = atm->lon[i];
    jur_atm->lat[i]  = atm->lat[i];
    jur_atm->p[i]    = atm->p[i];
    jur_atm->t[i]    = atm->t[i];
    for(int j = 0; j < NG; j++)
      jur_atm->q[j][i] = atm->q[j][i];
    for(int j = 0; j < NW; j++)
      jur_atm->k[j][i] = atm->k[j][i];
  }
}

void convert_obs_from_reference_to_unified(obs_t const *obs, jur_obs_t *jur_obs) {
  jur_obs->nr = obs->nr;
  for(int i = 0; i < obs->nr; i++) {
    jur_obs->time[i]   = obs->time[i];
    jur_obs->obsz[i]   = obs->obsz[i];
    jur_obs->obslon[i] = obs->obslon[i];
    jur_obs->obslat[i] = obs->obslat[i];
    jur_obs->vpz[i]    = obs->vpz[i];
    jur_obs->vplon[i]  = obs->vplon[i];
    jur_obs->vplat[i]  = obs->vplat[i];
    jur_obs->tpz[i]    = obs->tpz[i];
    jur_obs->tplon[i]  = obs->tplon[i];
    jur_obs->tplat[i]  = obs->tplat[i];
    for(int j = 0; j < ND; j++) {
      jur_obs->tau[i][j] = obs->tau[j][i];
      jur_obs->rad[i][j] = obs->rad[j][i];
    }
  }
}

void convert_obs_from_unified_to_reference(jur_obs_t const *jur_obs, obs_t *obs) {
  obs->nr = jur_obs->nr;
  for(int i = 0; i < obs->nr; i++) {
    obs->time[i]   = jur_obs->time[i];
    obs->obsz[i]   = jur_obs->obsz[i];
    obs->obslon[i] = jur_obs->obslon[i];
    obs->obslat[i] = jur_obs->obslat[i];
    obs->vpz[i]    = jur_obs->vpz[i];
    obs->vplon[i]  = jur_obs->vplon[i];
    obs->vplat[i]  = jur_obs->vplat[i];
    obs->tpz[i]    = jur_obs->tpz[i];
    obs->tplon[i]  = jur_obs->tplon[i];
    obs->tplat[i]  = jur_obs->tplat[i];
    for(int j = 0; j < ND; j++) {
      obs->tau[j][i] = jur_obs->tau[i][j];
      obs->rad[j][i] = jur_obs->rad[i][j];
    }
  }
}

void jur_unified_formod_multiple_packages(atm_t const *atm, obs_t *obs, int num_of_obs_packages, int32_t const *atm_id) {
  assert(omp_get_num_threads() == 1);
  jur_ctl_t *jur_ctl = jur_unified_init(0, NULL);
    
  int total_number_of_rays = 0;
  for(int i = 0; i < num_of_obs_packages; i++)
    total_number_of_rays += obs[i].nr;
  int const num_of_atms = jur_get_num_of_atms(total_number_of_rays, atm_id);
  jur_atm_t *jur_atm = (jur_atm_t*) malloc((size_t) num_of_atms * sizeof(jur_atm_t));
  for(int i = 0; i < num_of_atms; i++)
    convert_atm_from_reference_to_unified(&atm[i], &jur_atm[i]);
  
  jur_obs_t *jur_obs = (jur_obs_t*) malloc((size_t) num_of_obs_packages * sizeof(jur_obs_t));
  for(int i = 0; i < num_of_obs_packages; i++)
    convert_obs_from_reference_to_unified(&obs[i], &jur_obs[i]);

	jur_formod_multiple_packages(jur_ctl, jur_atm, jur_obs, num_of_obs_packages, atm_id, NULL);
  free(jur_atm);

  for(int i = 0; i < num_of_obs_packages; i++)
    convert_obs_from_unified_to_reference(&jur_obs[i], &obs[i]);
  free(jur_obs);
}
