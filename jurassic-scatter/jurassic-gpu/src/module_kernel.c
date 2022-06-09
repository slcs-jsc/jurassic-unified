#include "jurassic.h"
#include "sca_retrievalmodel.h"
#include "sca_scatter.h"
#include "sca_forwardmodel.h"

int main(int argc, char *argv[]) {
  
  static atm_t atm;

  static ctl_t ctl;
  
  static obs_t obs;

  static aero_t aero;
  
  gsl_matrix *k;
  
  size_t m, n;
  
  /* Check arguments... */
  if(argc<6)
    ERRMSG("Give parameters: <ctl> <obs> <atm> <aero> <kernel>");
  
  /* Read forward model control parameters... */
  jur_read_ctl(argc, argv, &ctl);

  // Initialization of the tables
  jur_table_initialization(&ctl); 
  
  /* Read observation geometry... */
  jur_read_obs(NULL, argv[2], &ctl, &obs);
  
  /* Read atmospheric data... */
  jur_read_atm(NULL, argv[3], &ctl, &atm);
  
  /* ============================================================= */
  /* Read aerosol and cloud data */
  if(strcmp(argv[4],"-")!=0 && ctl.sca_n>0) {
    jur_sca_read_aero(NULL, argv[4], &ctl, &aero);
    /* Get aerosol/cloud optical properties */
    jur_sca_get_opt_prop(&ctl, &aero);
  } 
  else if (strcmp(argv[4],"-")==0 && ctl.sca_n>0) { 
    ERRMSG("Please give aerosol file name or set SCA_N=0 for clear air simulation!");
  } 
  /* ============================================================= */

  /* Get sizes... */
  n=atm2x(&ctl, &atm, &aero, NULL, NULL, NULL);
  m=obs2y(&ctl, &obs, NULL, NULL, NULL);
  
  /* Check sizes... */
  if(n<=0)
    ERRMSG("No state vector elements!");
  if(m<=0)
    ERRMSG("No measurement vector elements!");
  
  /* Allocate... */
  k=gsl_matrix_alloc(m, n);
  
  /* Compute kernel matrix... */
  kernel(&ctl, &atm, &obs, &aero, k);
  
  /* Write matrix to file... */
  write_matrix(NULL, argv[5], &ctl, k, &atm, &aero, &obs, "y", "x", "r");
  
  /* Free... */
  gsl_matrix_free(k);
  
  return EXIT_SUCCESS;
}
