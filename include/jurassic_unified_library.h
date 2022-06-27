#ifndef JURASSIC_UNIFIED_LIBRARY_H
#define JURASSIC_UNIFIED_LIBRARY_H

  #include "jurassic.h"
  #include "../unified_library/jurassic_dimensions.h"
  #include "../unified_library/jurassic_structs.h" 

  #if LEN > LENMAX
    #error "Please compile jurassic-unified with larger LENMAX"
  #endif

  #if M > MMAX
    #error "Please compile jurassic-unified with larger MMAX"
  #endif

  #if N > NMAX
    #error "Please compile jurassic-unified with larger NMAX"
  #endif

  #if NQ > NQMAX
    #error "Please compile jurassic-unified with larger NQMAX"
  #endif

  #if NW > NWMAX
    #error "Please compile jurassic-unified with larger NWMAX"
  #endif

  #if ND > NDMAX
    #error "Please compile jurassic-unified with larger NDMAX"
  #endif

  #if NG > NGMAX
    #error "Please compile jurassic-unified with larger NGMAX"
  #endif

  #if NLOSMAX > NLOSMAX
    #error "Please compile jurassic-unified with larger NLOSMAX"
  #endif

  #if NP > NPMAX
    #error "Please compile jurassic-unified with larger NPMAX"
  #endif

  #if NR > NRMAX
    #error "Please compile jurassic-unified with larger NRMAX"
  #endif

  #if NSHAPE > NSHAPEMAX
    #error "Please compile jurassic-unified with larger NSHAPEMAX"
  #endif

  #if NFOV > NFOVMAX
    #error "Please compile jurassic-unified with larger NFOVMAX"
  #endif

  #if TBLNP > TBLNPMAX
    #error "Please compile jurassic-unified with larger TBLNPMAX"
  #endif

  #if TBLNS > TBLNSMAX
    #error "Please compile jurassic-unified with larger TBLNSMAX"
  #endif

  #if TBLNT > TBLNTMAX
    #error "Please compile jurassic-unified with larger TBLNTMAX"
  #endif

  #if TBLNU > TBLNUMAX
    #error "Please compile jurassic-unified with larger TBLNUMAX"
  #endif

  jur_ctl_t *jur_unified_init(int argc, char *argv[]);
	void jur_unified_formod_multiple_packages(atm_t const *atm, obs_t *obs, int num_of_obs_packages, int32_t const *atm_id);

#endif
