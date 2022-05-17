#include "sca_forwardmodel.h"
#include "sca_workqueue.h" /* Queue_Inactive, Queue_Prepare, Queue_Execute, Queue_Execute */
#include <assert.h> /* assert */
#include <omp.h>
#include <math.h>

#define __host__
#include "interface.h"


// declaration of the functions from CPUdrivers.c
double jur_continua_core_CPU(ctl_t const *ctl, pos_t const *los, int const id);
void jur_hydrostatic1d_CPU(ctl_t const *ctl, atm_t *atm, int const nr, int const ig_h2o);

/*****************************************************************************/

void jur_sca_copy_obs_row(obs_t const *source, int rs, obs_t *dest, int rd) {
  dest->time[rd] = source->time[rs];
  dest->obsz[rd] = source->obsz[rs];
  dest->obslon[rd] = source->obslon[rs];
  dest->obslat[rd] = source->obslat[rs];
  dest->vpz[rd] = source->vpz[rs];
  dest->vplon[rd] = source->vplon[rs];
  dest->vplat[rd] = source->vplat[rs];
  dest->tpz[rd] = source->tpz[rs];
  dest->tplon[rd] = source->tplon[rs];
  dest->tplat[rd] = source->tplat[rs];
  for(int i=0; i<ND; i++) {
    dest->tau[rd][i] = source->tau[rs][i]; //CHANGED
    dest->rad[rd][i] = source->rad[rs][i]; //CHANGED
  }
}

/*****************************************************************************/

void jur_sca_advanced_execute(ctl_t *ctl, atm_t *atm, aero_t *aero, queue_t *qs, int nr) {

  int *pref_sizes;
  ALLOC(pref_sizes, int, nr);
  pref_sizes[0] = 0;
  for(int i = 1; i < nr; i++)
    pref_sizes[i] = pref_sizes[i - 1] + qs[i - 1].end - qs[i - 1].begin;

  int sum_sizes = pref_sizes[nr - 1] + qs[nr - 1].end - qs[nr - 1].begin;
  
  int number_of_packages = (sum_sizes + NR - 1) / NR;
 
  double tic, toc;
  tic = omp_get_wtime();
  
  obs_t *obs_packages;
  ALLOC(obs_packages, obs_t, number_of_packages);
 
  int last_package_size = sum_sizes % NR;
  for(int i = 0; i < number_of_packages; i++)
    if(i == number_of_packages - 1 && last_package_size > 0)
      obs_packages[i].nr = last_package_size;
    else
      obs_packages[i].nr = NR;

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int i = 0; i < nr; i++) {
    for(int j = 0; j < qs[i].end - qs[i].begin; j++) {
      int ind;
      obs_t *obs;
      int package_id = (pref_sizes[i] + j) / NR;
      int obs_row = (pref_sizes[i] + j) % NR;
      jur_sca_get_queue_item(&qs[i], (void*)&obs, &ind, qs[i].begin + j);
      jur_sca_copy_obs_row(obs, ind, &obs_packages[package_id], obs_row);
    }
  }

  toc = omp_get_wtime();
  printf("TIMER #%d Execute: copy from queue to packages time: %lf\n", ctl->MPIglobrank, toc - tic);
  
  for(int i = 0; i < number_of_packages; i++)
    printf("%d ", obs_packages[i].nr);
  printf("\n");
  tic = omp_get_wtime();
	jur_formod_multiple_packages(ctl, atm, obs_packages, number_of_packages, NULL, aero);
  toc = omp_get_wtime();
  printf("TIMER #%d Execute: jur_formod_multiple_packages time: %lf\n", ctl->MPIglobrank, toc - tic);

  tic = omp_get_wtime();
 
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int i = 0; i < nr; i++)
    for(int j = 0; j < qs[i].end - qs[i].begin; j++) {
      int ind;
      obs_t *obs; 
      int package_id = (pref_sizes[i] + j) / NR;
      int obs_row = (pref_sizes[i] + j) % NR;
      jur_sca_get_queue_item(&qs[i], (void*)&obs, &ind, qs[i].begin + j);
      jur_sca_copy_obs_row(&obs_packages[package_id], obs_row, obs, ind);
    }

  toc = omp_get_wtime();
  printf("TIMER #%d Execute: copy from packages queue time: %lf\n", ctl->MPIglobrank, toc - tic);
  free(pref_sizes);
  free(obs_packages);
}

/*****************************************************************************/

void jur_sca_formod(ctl_t *ctl,
	    atm_t *atm,
	    obs_t *obs,
	    aero_t *aero) {
  
  static int mask[NDMAX][NRMAX];
  
  int id, ir;
  double tic, toc;
  
  queue_t *qs;

  /* Save observation mask... */
  for(id=0; id<ctl->nd; id++)
    for(ir=0; ir<obs->nr; ir++)
      mask[id][ir]=!gsl_finite(obs->rad[ir][id]); //CHANGED
  
  /* Hydrostatic equilibrium... */
  // WARNING: ctl->ctm_h2o is at the moment not exectly the same as in jurassic-GPU
  // because of different read_ctl function 
  static int ig_h2o = -999;
  if((ctl->ctm_h2o) && (-999 == ig_h2o)) ig_h2o = jur_find_emitter(ctl, "H2O");

  jur_hydrostatic1d_CPU(ctl, atm, obs->nr, ig_h2o); // in this call atm might get modified
  
  /* Particles: Calculate optical properties in retrieval */
  if(ctl->retnn || ctl->retrr || ctl->retss) {
    jur_sca_get_opt_prop(ctl, aero);
  }

  if (abs(ctl->leaf_nr) > 0) { /* switch usage of queue on by setting |MAX_QUEUE| > 0 */
    
    int leaf_rays_per_ray = ctl->leaf_nr / obs->nr + 1;
    if(ctl->leaf_nr < 0) leaf_rays_per_ray = (-1 * ctl->leaf_nr) / obs->nr + 1;
    ALLOC(qs, queue_t, obs->nr);
    for(int i = 0; i < obs->nr; i++)
      jur_sca_init_queue(&qs[i], leaf_rays_per_ray);
    
    printf("DEBUG #%d %s init %d queues with %d elements\n", ctl->MPIglobrank,  __func__, obs->nr, leaf_rays_per_ray);
    /*
     *  Work Queue Architecture with three stages:
     *    Pp Prepare constructs the multiple scattering tree 
     *               and pushes work items into a queue
     *    x  Execute performs the work in the queue (lowest 
     *               tree level) without further scattering
     *    Cc Collect traverses the tree again and collects
     *               the results
     */
    ctl->queue_state = Queue_Prepare; /* activate the work queue */ 
    printf("DEBUG #%d %s start Queue_Prepare\n", ctl->MPIglobrank, __func__);


    tic = omp_get_wtime();

    /* Do first ray path sequential (to initialize model)... */
    jur_sca_formod_pencil(ctl, atm, obs, aero, ctl->sca_mult, 0, &qs[0]);
    
    toc = omp_get_wtime();
    printf("TIMER #%d Prepare 1st part time: %lf\n", ctl->MPIglobrank, toc - tic);

    tic = omp_get_wtime();
    /* Do remaining ray paths in parallel... */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for(int i=1; i<obs->nr; i++){
      jur_sca_formod_pencil(ctl, atm, obs, aero, ctl->sca_mult, i, &qs[i]);
    }
    printf("\n");

    toc = omp_get_wtime();
    printf("TIMER #%d Prepare 2nd part time: %lf\n", ctl->MPIglobrank, toc - tic);


  } else {
    // I had to add this line to avoid warnings
    qs = NULL;

    ctl->queue_state = Queue_Inactive; /* deactivate the work queue */
    printf("DEBUG #%d %s Queue_Inactive\n", ctl->MPIglobrank,  __func__);

    tic = omp_get_wtime();
    /* Do first ray path sequential (to initialize model)... */
    jur_sca_formod_pencil(ctl, atm, obs, aero, ctl->sca_mult, 0, NULL);
    
    /* Do remaining ray paths in parallel... */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for(ir=1; ir<obs->nr; ir++){
      jur_sca_formod_pencil(ctl, atm, obs, aero, ctl->sca_mult, ir, NULL);
    }

    toc = omp_get_wtime();
    printf("TIMER #%d Queue inactive execution time: %lf\n", ctl->MPIglobrank, toc - tic);
  }
  
  if (Queue_Prepare == ctl->queue_state) {
      tic = omp_get_wtime();
      if (ctl->leaf_nr < 0) { /* execute on CPU */
        ctl->queue_state = Queue_Execute;
        printf("DEBUG #%d %s start Queue_Execute [%d, %d) on CPU\n", ctl->MPIglobrank, __func__, 0, obs->nr);
        printf("DEBUG #%d only scatter CPU-execute version\n", ctl->MPIglobrank);
        /* Do first ray path sequential (to initialize model)... */
        
        for(int i = 0; i < 1; i++) {
          for(ir = qs[i].begin; ir < qs[i].begin + 1; ++ir) {
            jur_sca_formod_pencil(ctl, atm, NULL, aero, 0, ir, &qs[0]);
          }
          for(ir = qs[i].begin + 1; ir < qs[i].end; ++ir) {
            jur_sca_formod_pencil(ctl, atm, NULL, aero, 0, ir, &qs[0]);
          }
        }/* ir-loop */
        /* Do remaining ray paths in parallel... */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
        for(int i = 1; i < obs->nr; i++) {
          for(int j = qs[i].begin; j < qs[i].end; j++) {
            jur_sca_formod_pencil(ctl, atm, NULL, aero, 0, j, &qs[i]);
          } /* ir-loop */
        }
      } else {
        printf("DEBUG #%d Call advanced execute!\n", ctl->MPIglobrank);
        jur_sca_advanced_execute(ctl, atm, aero, qs, obs->nr);
        //ERRMSG("No GPU version of jur_sca_formod_pencil implemented!");
      }
      toc = omp_get_wtime();
      printf("TIMER #%d Execute time: %lf\n", ctl->MPIglobrank, toc - tic);
      
      tic = omp_get_wtime();
      ctl->queue_state = Queue_Collect;
      printf("DEBUG #%d %s start Queue_Collect\n", ctl->MPIglobrank, __func__);


      for(int i=0; i<1; i++) {
        jur_sca_formod_pencil(ctl, atm, obs, aero, ctl->sca_mult, i, &qs[i]);
      }
      toc = omp_get_wtime();
      printf("TIMER #%d Collect-1st part time: %lf\n", ctl->MPIglobrank, toc - tic);

      tic = omp_get_wtime();

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
      for(int i=1; i<obs->nr; i++) {
        jur_sca_formod_pencil(ctl, atm, obs, aero, ctl->sca_mult, i, &qs[i]);
      }
 
      toc = omp_get_wtime();
      printf("TIMER #%d Collect-2nd part time: %lf\n", ctl->MPIglobrank, toc - tic);

      tic = omp_get_wtime();
      for(int i = 0; i < obs->nr; i++)
        jur_sca_init_queue(&qs[i], -1); 
      free(qs);
      ctl->queue_state = Queue_Inactive; /* done */
      toc = omp_get_wtime();
      printf("TIMER #%d Memory free time: %lf\n", ctl->MPIglobrank, toc - tic);
  }

  /* Apply field-of-view convolution... */
  jur_sca_formod_fov(ctl, obs);
  
  /* Convert radiance to brightness temperature... */
  if(ctl->write_bbt)
    for(ir=0; ir<obs->nr; ir++)
      for(id=0; id<ctl->nd; id++)
	      obs->rad[ir][id]=brightness_core(obs->rad[ir][id], ctl->nu[id]); //CHANGED

  /* Apply observation mask... */
  for(id=0; id<ctl->nd; id++)
    for(ir=0; ir<obs->nr; ir++)
      if(mask[id][ir])
	obs->rad[ir][id]=GSL_NAN; //CHANGED
}

/*****************************************************************************/

void jur_sca_formod_fov(ctl_t *ctl,
		obs_t *obs) {
  
  static obs_t obs2;
  
  static double dz[NSHAPE], rad[NDMAX][NRMAX], tau[NDMAX][NRMAX],
    w[NSHAPE], wsum, z[NRMAX], zfov;
  
  static int init=0, i, id, idx, ir, ir2, n, nz;
  
  /* Do not take into account FOV... */
  if(ctl->fov[0]=='-')
    return;
  
  /* Initialize FOV data... */
  if(!init) {
    init=1;
    jur_sca_read_shape(ctl->fov, dz, w, &n);
  }
  
  /* Copy observation data... */
  jur_copy_obs(ctl, &obs2, obs, 0);
  
  /* Loop over ray paths... */
  for(ir=0; ir<obs->nr; ir++) {
    
    /* Get radiance and transmittance profiles... */
    nz=0;
    for(ir2=GSL_MAX(ir-NFOV, 0); ir2<GSL_MIN(ir+1+NFOV, obs->nr); ir2++)
      if(obs->time[ir2]==obs->time[ir]) {
	z[nz]=obs2.tpz[ir2];
	for(id=0; id<ctl->nd; id++) {
	  rad[id][nz]=obs2.rad[ir2][id]; //CAHNGED
	  tau[id][nz]=obs2.tau[ir2][id]; //CHANGED
	}
	nz++;
      }
    if(nz<2)
      ERRMSG("Cannot apply FOV convolution!");
    
    /* Convolute profiles with FOV... */
    wsum=0;
    for(id=0; id<ctl->nd; id++) {
      obs->rad[ir][id]=0; //CHANGED
      obs->tau[ir][id]=0; //CHANGED
    }
    for(i=0; i<n; i++) {
      zfov=obs->tpz[ir]+dz[i];
      idx=jur_locate(z, nz, zfov);
      for(id=0; id<ctl->nd; id++) {
        obs->rad[ir][id]+=w[i] //CHANGED
          *LIN(z[idx], rad[id][idx], 
	       z[idx+1], rad[id][idx+1], zfov);
        obs->tau[ir][id]+=w[i] //CHANGED
          *LIN(z[idx], tau[id][idx],
	       z[idx+1], tau[id][idx+1], zfov);
      }
      wsum+=w[i];
    }
    for(id=0; id<ctl->nd; id++) {
      obs->rad[ir][id]/=wsum; //CHANGED
      obs->tau[ir][id]/=wsum; //CHANGED
    }
  }
}

/*****************************************************************************/

void jur_sca_formod_pencil(ctl_t *ctl,
                   atm_t *atm,
                   obs_t *obs,
                   aero_t *aero,
                   int scattering,
                   int ir,
                   queue_t *q) {

  static int init=0;

  // removing tbl_t
  // static tbl_t *tbl;
  static trans_table_t *trans_tbl;

  pos_t *los;
  int np = 0;
  double tsurf;
  los = NULL; // because it has to be initialized 
  obs_t *obs2;
  
  // tau_path[NGMAX][NDMAX] --> tau_path[NDMAX][NGMAX], because of jurassic-gpu
  double beta_ctm[NDMAX], beta_ext_tot, dx[3], eps, src_all, src_planck[NDMAX],
    src_sca[NDMAX], tau_path[NDMAX][NGMAX], tau_gas[NDMAX], x[3], x0[3], x1[3];

  int i, id, ip, ip0, ip1, ig;

  int const Queue_Prepare_Leaf = Queue_Prepare << 1;
  int const Queue_Execute_Leaf = Queue_Execute << 1;
  int const Queue_Collect_Leaf = Queue_Collect << 1;
#ifdef WORK_QUEUE
  int const queue_mode = ctl->queue_state << (0 == scattering);
#else
  int const queue_mode = Queue_Inactive; /* Queue_Inactive == -1 */
#endif

#ifdef  FORMOD_DEBUG
  printf("# %s(..., %p, aero, scattering=%d, ir=%d) queue_mode = %d;\n", __func__, (void*)obs, scattering, ir, queue_mode);
#endif 
  
if ((Queue_Collect|Queue_Prepare) & queue_mode) { /* CPp */
  /* Allocate... */
  los = (pos_t*) malloc((NLOS) * sizeof(pos_t));

  /* Raytracing... */
  np = traceray(ctl, atm, obs, ir, los, &tsurf, aero, 1); // with scattering
} /* CPp */

if (Queue_Prepare_Leaf == queue_mode) { /* ==p */
  i = jur_sca_push_queue(q, (void*)obs, ir); /* push input and pointer to output */
  if (i < 0) ERRMSG("Too many queue items!"); /* failed */
  return;
} /* ==p */

if (Queue_Collect_Leaf == queue_mode) { /* ==c */
  jur_sca_pop_queue(q, (void*)&obs2, &ir); /* pop result */
  /* Copy results... */
  for(id=0; id<ctl->nd; id++) {
    obs->rad[ir][id] = obs2->rad[ir][id]; //CHANGED
    obs->tau[ir][id] = obs2->tau[ir][id]; //CAHNGED
  } /* id */

  // free obs2 only in simulations with scattering
  // in clear-air case obs2 equals obs!
  if(ctl->sca_mult > 0 && obs2->nr - 1 == ir) free(obs2);
  return;
} /* ==c */

if (Queue_Execute_Leaf == queue_mode) { /* ==x */
  jur_sca_get_queue_item(q, (void*)&obs, &ir, ir); /* get input */
  los = (pos_t*) malloc((NLOS) * sizeof(pos_t));
  np = traceray(ctl, atm, obs, ir, los, &tsurf, aero, 1); // with scattering
} /* ==x */

if ((Queue_Collect|Queue_Execute_Leaf) & queue_mode) { /* Cx */
  /* Read tables... */
  if(!init) {
    init=1;
    // removing tbl_t
    // tbl = scatter_get_tbl(ctl);
    // printf("%d\n", tbl->np[0][0]); // Have to do it, because of unused warning... 

    // bug that I had:
    // https://stackoverflow.com/questions/8552684/pointer-return-value-changes-after-function-call
    trans_tbl = get_tbl(ctl); 
  }

  /* Initialize... */
  for(id=0; id<ctl->nd; id++) {
    obs->rad[ir][id]=0.0; //CHANGED
    obs->tau[ir][id]=1.0; //CHANGED
    // added for jurassic-gpu tau_gas
    for(ig = 0; ig < NG; ig++) { // loop over gases
      tau_path[id][ig] = 1.0;
    } // ig
  }
} /* Cx */


  /* Loop over LOS points... */
  for(ip=0; ip<np; ip++) {

if ((Queue_Collect|Queue_Execute_Leaf) & queue_mode) { /* Cx */
    /* Get trace gas transmittance... */
    for(id = 0; id < ctl->nd; id++)
      tau_gas[id] = apply_ega_core(trans_tbl, &los[ip], tau_path[id], ctl->ng, id);

    /* Get continuum absorption... */
    for(id = 0; id < ctl->nd; id++)
      beta_ctm[id] = jur_continua_core_CPU(ctl, &los[ip], id) / los[ip].ds;

    /* Compute Planck function... */
    for(id = 0; id < ctl->nd; id++) {
      src_planck[id] = src_planck_core(trans_tbl, los[ip].t, id);
    }

} /* Cx */

    /* Compute radiative transfer with scattering source... */
    if(scattering>0 && los[ip].aerofac>0) {

if ((Queue_Collect|Queue_Prepare) & queue_mode) { /* CP */
      /* Compute scattering source term... */
      jur_geo2cart(los[ip].z, los[ip].lon, los[ip].lat, x);
      ip0=(ip>0 ? ip-1 : ip);
      ip1=(ip<np ? ip+1 : ip);
      jur_geo2cart(los[ip0].z, los[ip0].lon, los[ip0].lat, x0);
      jur_geo2cart(los[ip1].z, los[ip1].lon, los[ip1].lat, x1);
      for(i=0; i<3; i++)
        dx[i]=x1[i]-x0[i];

      jur_sca_srcfunc_sca(ctl,atm,aero,obs->time[ir],x,dx,los[ip].aeroi,src_sca,scattering, q);
} /* CP */

if ((Queue_Collect) & queue_mode) { /* C */
      /* Loop over channels... */
      for(id=0; id<ctl->nd; id++)
        if(tau_gas[id]>0) {

          /* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
          /* Get gas and aerosol/cloud extinctions... */
          beta_ext_tot = (-1.)*log(tau_gas[id])/los[ip].ds + beta_ctm[id] +
                         los[ip].aerofac*aero->beta_e[los[ip].aeroi][id];

          /* enthält tau_gas bereits k????????? */

          /* Get segment emissivity */
          eps = 1-exp(-1*beta_ext_tot*los[ip].ds);

          /* Compute weighted segment source */
          src_all=((beta_ext_tot - los[ip].aerofac*aero->beta_s[los[ip].aeroi][id]) *
                   src_planck[id] +
                   los[ip].aerofac*aero->beta_s[los[ip].aeroi][id]*src_sca[id]) /
                  beta_ext_tot;

          /* Compute radiance: path extinction * segment emissivity * segment source */
          obs->rad[ir][id] += obs->tau[ir][id]*eps*src_all; //CHANGED

          /* Compute path transmittance... */
          obs->tau[ir][id] *= exp(-1.*beta_ext_tot*los[ip].ds); //CAHNGED
        }
} /* C */
    }

    /* Compute radiative transfer without scattering source... */
    else {

if ((Queue_Collect|Queue_Execute_Leaf) & queue_mode) { /* Cx */
      /* Loop over channels... */
      for(id=0; id<ctl->nd; id++)
        if(tau_gas[id]>0) {

          /* Get segment emissivity... */
          if (ctl->sca_n==0 || los[ip].aerofac==0 ) {
            eps=1-tau_gas[id]*exp(-1. * beta_ctm[id] * los[ip].ds);
          }
          else if (strcmp(ctl->sca_ext, "beta_a") == 0) {

            eps=1-tau_gas[id]*exp(-1. * (beta_ctm[id] + los[ip].aerofac*
                                         aero->beta_a[los[ip].aeroi][id])*los[ip].ds);
          } else {
            eps=1-tau_gas[id]*exp(-1. * (beta_ctm[id] + los[ip].aerofac*
                                         aero->beta_e[los[ip].aeroi][id])*los[ip].ds);
          }

          /* Compute radiance... */
          obs->rad[ir][id]+=src_planck[id]*eps*obs->tau[ir][id]; //CHANGED

          /* Compute path transmittance... */
          obs->tau[ir][id]*=(1-eps); //CAHNGED
        }
    }
} /* Cx */
  }

if ((Queue_Collect|Queue_Execute_Leaf) & queue_mode) { /* Cx */
  /* Add surface... */
  if(tsurf>0) {
    /* Compute Planck function */
    for(id = 0; id < ctl->nd; id++) {
      src_planck[id] = src_planck_core(trans_tbl, tsurf, id);
    }

    for(id=0; id<ctl->nd; id++)
      obs->rad[ir][id]+=src_planck[id]*obs->tau[ir][id]; //CHANGED
  }

} /* Cx */


  if ((Queue_Collect|Queue_Execute_Leaf|Queue_Prepare) & queue_mode) { 
    /* Free... */
    free(los);
  }
}

/*****************************************************************************/

void jur_sca_read_shape(const char *filename,
		double *x,
		double *y,
		int *n) {
  
  FILE *in;
  
  char line[LEN];
  
  /* Write info... */
  /*printf("Read shape function: %s\n", filename); */
  
  /* Open file... */
  if(!(in=fopen(filename, "r")))
    ERRMSG("Cannot open file!");
  
  /* Read data... */
  *n=0;
  while(fgets(line, LEN, in))
    if(sscanf(line,"%lg %lg", &x[*n], &y[*n])==2)
      if((++(*n))>NSHAPE)
        ERRMSG("Too many data points!");
  
  /* Check number of points... */
  if(*n<1)
    ERRMSG("Could not read any data!");
  
  /* Close file... */
  fclose(in);
}

