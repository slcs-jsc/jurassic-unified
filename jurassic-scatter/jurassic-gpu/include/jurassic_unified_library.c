#include "jurassic_unified_library.h"

// declarations of functions from unified..
void jur_get_dimensions(int32_t *dimensions);

jur_ctl_t *jur_unified_init(char const *ctl_name, int argc, char *argv[]) {
  static jur_ctl_t *jur_ctl = NULL;
  if(NULL == jur_ctl) {
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

    jur_read_ctl(0, NULL, jur_ctl); 
    jur_table_initialization(&jur_ctl); 
  }
}

void convert_atm_from_reference_to_unified(atm_t const *atm, jur_atm_t *jur_atm) {
  //...

}

void jur_unified_formod_multiple_packages(char const *ctl_name, atm_t const *atm, obs_t *obs, int num_of_obs, int32_t const *atm_id) {
  assert(omp_num_threads() == 1);
  jur_ctl_t *jur_ctl = jur_unified_init(ctl_name, 0, NULL);
    
  // TODO:
  num_of_atms = jur_get_num_of_atms(int const nr, atm_id);
  jur_atm_t *jur_atm = (jur_atm_t*) malloc(num_of_atms * sizeof(jur_atm_t));
  for(int i = 0; i < n; i++)
    for(int j = 0; j < num_of_atms_i; j++)
      convert_atm_from_reference_to_unified(atm, jur_atm);
  
  jur_obs_t *jur_obs = (jur_obs_t*) malloc(n * sizeof(jur_obs_t));
  convert_obs_from_reference_to_unified(obs, jur_obs);

	jur_formod_multiple_packages(jur_ctl, jur_atm, jur_obs, n, atm_id);
  free(jur_atm);

  convert_obs_from_unified_to_reference(jur_obs, obs);
  free(jur_obs);
}
