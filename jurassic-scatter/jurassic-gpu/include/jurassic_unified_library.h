#ifndef JURASSIC_UNIFIED_LIBRARY_H
#define JURASSIC_UNIFIED_LIBRARY_H
  #include "jurassic.h"
  #include "../unified_library/jurassic_dimensions.h"
  #include "../unified_library/jurassic_structs.h"

  #if ND > NDMAX
    #error "Please compile jurassic-unified with larger NDMAX"
  #endif

  #if ND > NDMAX
    #error "Please compile jurassic-unified with larger NDMAX"
  #endif

  #if ND > NDMAX
    #error "Please compile jurassic-unified with larger NDMAX"
  #endif

  ...

  jur_ctl_t *jur_unified_init(char const *ctl_name, int argc, char *argv[]);

	void jur_unified_formod_multiple_packages(char const *ctl_name, atm_t const *atm, obs_t *obs, int num_of_obs, int32_t const *atm_id); // without aero

#endif
