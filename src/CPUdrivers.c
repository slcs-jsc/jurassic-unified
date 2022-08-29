#include "jr_common.h" // inline definitions of common functions for CPU and GPU code

// ################ CPU driver routines - keep consistent with GPUdrivers.cu ##############

__host__
void jur_radiance_to_brightness_CPU(ctl_t const *ctl, obs_t *obs) {
#pragma omp parallel for
  for(int ir = 0; ir < obs->nr; ir++) { // loop over rays
    for(int id = 0; id < ctl->nd; id++) { // loop over detectors
      // convert in-place
      obs->rad[ir][id] = jur_brightness_core(obs->rad[ir][id], ctl->nu[id]);
    } // id
  } // ir
} // jur_radiance_to_brightness

__host__
void jur_surface_terms_CPU(trans_table_t const *tbl, obs_t *obs, double const tsurf[], int const nd) {
#pragma omp parallel for
  for(int ir = 0; ir < obs->nr; ir++) { // loop over rays
    for(int id = 0; id < nd; id++) { // loop over detectors
      jur_add_surface_core(obs, tbl, tsurf[ir], ir, id);
    } // id
  } // ir
} // jur_surface_terms_CPU

__host__
double jur_continua_core_CPU(ctl_t const *ctl, pos_t const *los, int const id) {
  static int init = 0;
  static int CO2, H2O, N2, O2, ig_co2, ig_h2o;
  if(!init) {
#ifdef _OPENMP
#pragma omp critical
#endif
    if(!init) {
      ig_co2 = -999, ig_h2o = -999;
      if((ctl->ctm_h2o) && (-999 == ig_h2o)) ig_h2o = jur_find_emitter(ctl, "H2O");
      if((ctl->ctm_co2) && (-999 == ig_co2)) ig_co2 = jur_find_emitter(ctl, "CO2");
      CO2 = (1 == ctl->ctm_co2) && (ig_co2 >= 0);
      H2O = (1 == ctl->ctm_h2o) && (ig_h2o >= 0);
      N2 = (1 == ctl->ctm_n2);
      O2 = (1 == ctl->ctm_o2);
      init = 1;
    }
  }

  double const p = los->p;
  double const t = los->t;
  double const ds = los->ds;
  double beta_ds = los->k[ctl->window[id]]*ds;  // extinction
  // make sure that ig_co2 and ig_h2o are both >= 0
  if(CO2) beta_ds += jur_continua_ctmco2(ctl->nu[id], p, t, los->u[ig_co2]);  // co2 continuum
  if(H2O) beta_ds += jur_continua_ctmh2o(ctl->nu[id], p, t, los->q[ig_h2o], los->u[ig_h2o]);  // h2o continuum
  if(N2)  beta_ds += jur_continua_ctmn2(ctl->nu[id], p, t)*ds;  // n2 continuum
  if(O2)  beta_ds += jur_continua_ctmo2(ctl->nu[id], p, t)*ds;  // o2 continuum
  return     beta_ds;
} // jur_continua_core_CPU

__host__
void jur_apply_kernels_CPU(trans_table_t const *tbl, ctl_t const *ctl, obs_t *obs,
    pos_t (*restrict los)[NLOSMAX], int const np[],
    double const (*restrict aero_beta)[NDMAX]) { // aero_beta is added

#pragma omp parallel for
  for(int ir = 0; ir < obs->nr; ir++) { // loop over independent rays
    double tau_path[NDMAX][NGMAX]; // private for each ray
    for(int id = 0; id < NDMAX; id++) { // loop over detectors
      obs->rad[ir][id] = 0.0;
      obs->tau[ir][id] = 1.0;
      for(int ig = 0; ig < NGMAX; ig++) { // loop over gases
        tau_path[id][ig] = 1.0;
      } // ig
    } //  id

    for(int ip = 0; ip < np[ir]; ++ip) { // loop over line-of-sight points
      for(int id = 0; id < ctl->nd; id++) { // loop over detector channels

        // compute extinction coefficient
        double const beta_ds = jur_continua_core_CPU(ctl, &(los[ir][ip]), id);

        // Bug found here:
        double aero_ds = 0;
        if(NULL != aero_beta) // only if scattering is included
          aero_ds = los[ir][ip].aerofac * aero_beta[los[ir][ip].aeroi][id] * los[ir][ip].ds;

        // compute transmission with the EGA method
        double const tau_gas = jur_apply_ega_core(tbl, &(los[ir][ip]), tau_path[id], ctl->ng, id);
        // compute the source term
        double const planck = jur_src_planck_core(tbl, los[ir][ip].t, id);
        // perform integration
        jur_new_obs_core(obs, ir, id, beta_ds + aero_ds, planck, tau_gas);

      } // id --> could be vectorized over detector channels
#ifdef GPUDEBUG
      assert(los[ir][ip].ip == ip); // sanity check
      assert(los[ir][ip].ir == ir); // sanity check
#endif
    } // ip --> non-parallelisable due to loup carried dependency
  } // ir --> OpenMP parallel over rays

} // jur_apply_kernels_CPU

__host__
void jur_raytrace_rays_CPU(ctl_t const *ctl, atm_t const *atm, obs_t *obs,
    pos_t los[NRMAX][NLOSMAX], double tsurf[], int np[],
    int32_t const *atm_id, aero_t const *aero, int scattering_included) {
#pragma omp parallel for
  for(int ir = 0; ir < obs->nr; ir++) { // loop over rays
    np[ir] = jur_traceray(ctl, &atm[(NULL == atm_id ? 0 : atm_id[ir])], obs, ir, los[ir], &tsurf[ir], aero, scattering_included);
  } // ir
} // jur_raytrace_rays_CPU

__host__
void jur_hydrostatic1d_CPU(ctl_t const *ctl, atm_t *atm, int const num_of_atms, int const ig_h2o) {
  if(ctl->hydz < 0) return; // Check reference height
  for(int i = 0; i < num_of_atms; i++) { // Apply hydrostatic equation to individual profiles
    jur_hydrostatic_1d_h2o(ctl, &atm[i], 0, atm[i].np, ig_h2o);
  } // i
} // jur_hydrostatic1d_CPU

// ################ end of CPU driver routines ##############

// The full forward model on the CPU ////////////////////////////////////////////
__host__
void jur_formod_one_package_CPU(ctl_t const *ctl, atm_t *atm, obs_t *obs,
    int32_t const *atm_id, aero_t const *aero) { // NULL == atm_id if all observations use the same atm
  // otherwise |atm_id| == obs -> nr
  printf("DEBUG #%d jur_formod_one_package_CPU was called!\n", ctl->MPIglobrank);

  if (ctl->checkmode) {
    printf("# %s: checkmode = %d, no actual computation is performed!\n", __func__, ctl->checkmode);
    return; // do nothing here
  } // checkmode

  assert(obs);

  char mask[NRMAX][NDMAX];
  jur_save_mask(mask, obs, ctl);

  trans_table_t const *tbl = jur_get_tbl(ctl);
  double *t_surf = (double*)malloc((size_t) obs->nr * sizeof(double));
  int *np = (int*)malloc((size_t) obs->nr * sizeof(int));
  pos_t (*los)[NLOSMAX] = (pos_t (*)[NLOSMAX])malloc((size_t) (obs->nr * NLOSMAX) * sizeof(pos_t));

  static int ig_h2o = -999;
  if((ctl->ctm_h2o) && (-999 == ig_h2o)) ig_h2o = jur_find_emitter(ctl, "H2O");

  jur_hydrostatic1d_CPU(ctl, atm, jur_get_num_of_atms(obs->nr, atm_id), ig_h2o); // in this call atm might get modified
  // if formod function was NOT called from jurassic-scatter project
  if(NULL != aero && ctl->sca_n > 0) { // only if scattering is included
    jur_raytrace_rays_CPU(ctl, atm, obs, los, t_surf, np, atm_id, aero, 1);
  } else {
    jur_raytrace_rays_CPU(ctl, atm, obs, los, t_surf, np, atm_id, NULL, 0);
  }

  // "beta_a" -> 'a', "beta_e" -> 'e'
  char const beta_type = ctl->sca_ext[5];

  jur_apply_kernels_CPU(tbl, ctl, obs, los, np,
      NULL == aero ? NULL : ('a' == beta_type ? aero->beta_a : aero->beta_e));
  jur_surface_terms_CPU(tbl, obs, t_surf, ctl->nd);

  free(los);
  free(np);
  free(t_surf);

  if(ctl->write_bbt && ctl->leaf_nr == -1)  // convert radiance to brightness (in-place)
    jur_radiance_to_brightness_CPU(ctl, obs);

  jur_apply_mask(mask, obs, ctl);
} // jur_formod_one_package_CPU

__host__
void jur_formod_multiple_packages_CPU(ctl_t const *ctl, atm_t *atm, obs_t *obs, int n, int32_t const *atm_id, aero_t const *aero) {
  if(NULL == atm_id || 1 == n) {
#pragma omp parallel for
    for(int i = 0; i < n; i++) {
      jur_formod_one_package_CPU(ctl, atm, &obs[i], atm_id, aero);
    }
  }
  else {
    atm_t **divided_atms = (atm_t **) malloc((size_t) n * sizeof(atm_t *));
    int32_t **divided_atm_ids = (int32_t **) malloc((size_t) n * sizeof(int32_t *));

    jur_divide_atm_data_into_packages(atm, obs, n, atm_id, divided_atms, divided_atm_ids);

#pragma omp parallel for
    for(int i = 0; i < n; i++) {
      jur_formod_one_package_CPU(ctl, divided_atms[i], &obs[i], divided_atm_ids[i], aero);
    }

    for(int i = 0; i < n; i++) {
      free(divided_atms[i]);
      free(divided_atm_ids[i]);
    }
    free(divided_atms);
    free(divided_atm_ids);
  }
} // jur_formod_multiple_packages_CPU

  __host__
void jur_formod_multiple_packages_GPU(ctl_t const *ctl, atm_t *atm, obs_t *obs,
    int n, int32_t const *atm_id, aero_t const *aero)
#ifdef hasGPU
  ; // declaration only, will be provided by GPUdrivers.o at link time
#else
{ // definition here
  static int warnGPU = 1;
  if (ctl->useGPU > 0) { // USEGPU > 0 means use-GPU-always
    fprintf(stdout, "CUDA not found during compilation, cannot run on GPUs. ABORT\n\n");
    fprintf(stdout, "USEGPU = 1 (use-GPU-always) found in controls\n"
        "USEGPU = -1 (use-GPU-if-possible) could help\n\n");
    fflush(stdout);
    fprintf(stderr, "CUDA not found during compilation, cannot run on GPUs. ABORT\n\n");
    exit(EXIT_FAILURE);
  } else {               // USEGPU < 0 means use-GPU-if-possible
    assert(ctl->useGPU < 0);
    // automatic decision: fallback solution is CPU
    if (warnGPU) { // this warning appears only once per process
      printf("CUDA not found during compilation, continue on CPUs instead!\n");
      warnGPU = 0; // switch this warning off
    } // warnGPU
    jur_formod_multiple_packages_CPU(ctl, atm, obs, n, atm_id, aero);
  } //
} // jur_formod_multiple_packages_GPU
#endif

__host__
void jur_formod_multiple_packages(ctl_t const *ctl, atm_t *atm, obs_t *obs, int n, int32_t const *atm_id, aero_t const *aero) {
  printf("DEBUG #%d jur_formod: number of packages.. %d\n", ctl->MPIglobrank, n);
  if (ctl->checkmode) {
    static int nr_last_time = -999;
    if (obs->nr != nr_last_time) {
      printf("# %s: %d max %d rays , %d of max %d gases, %d of max %d channels\n",
          __func__, obs->nr, NRMAX, ctl->ng, NGMAX, ctl->nd, NDMAX);
      // the number of rays can be zero if we skipped to read obs.tab, checkmode=1
      nr_last_time = obs->nr;
    } // only report if nr changed
  } // checkmode
  if (ctl->useGPU) { //has to be changed
    jur_formod_multiple_packages_GPU(ctl, atm, obs, n, atm_id, aero);
  } else { // USEGPU = 0 means use-GPU-never
    jur_formod_multiple_packages_CPU(ctl, atm, obs, n, atm_id, aero);
  } // useGPU
} // jur_formod_multiple_packages

__host__
void jur_formod(ctl_t const *ctl, // function with the original parameters, without aero, one atm, one obs package
    atm_t *atm,
    obs_t *obs) {
  jur_formod_multiple_packages(ctl, atm, obs, 1, NULL, NULL);
} // jur_formod

//we could use the same trick as above but it's not necessary
  __host__
trans_table_t* jur_get_tbl_on_GPU(ctl_t const *ctl)
#ifdef hasGPU
  ; // declaration only, will be provided by GPUdrivers.o at link time
#else
{ // definition here
  static int warnGPU = 1;
  if (ctl->useGPU > 0) { // USEGPU > 0 means use-GPU-always
    fprintf(stdout, "CUDA not found during compilation, cannot run on GPUs. ABORT\n\n");
    fprintf(stdout, "USEGPU = 1 (use-GPU-always) found in controls\n"
        "USEGPU = -1 (use-GPU-if-possible) could help\n\n");
    fflush(stdout);
    fprintf(stderr, "CUDA not found during compilation, cannot run on GPUs. ABORT\n\n");
    exit(EXIT_FAILURE);
  } else {               // USEGPU < 0 means use-GPU-if-possible
    assert(ctl->useGPU < 0);
    // automatic decision: fallback solution is CPU
    if (warnGPU) { // this warning appears only once per process
      printf("CUDA not found during compilation, continue on CPUs instead!\n");
      warnGPU = 0; // switch this warning off
    } // warnGPU
    printf("DEBUG #%d call initilaze CPU..\n", ctl->MPIglobrank);
    return jur_get_tbl_core(ctl);
  } //
} // jur_get_tbl_on_GPU
#endif

__host__
void jur_table_initialization(ctl_t const *ctl) {
  jur_get_tbl(ctl);
} // jur_table_initialization

__host__
trans_table_t* jur_get_tbl(ctl_t const *ctl) {
  trans_table_t *ret = NULL;
  if(ctl->useGPU)
    ret = jur_get_tbl_on_GPU(ctl);
  else
    ret = jur_get_tbl_core(ctl);
  return ret;
} // jur_table_initialization
