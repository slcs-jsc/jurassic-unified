#include "control.h"

#define __host__
#include "interface.h"

/* Module to read and handle all input parameters. */

/*****************************************************************************/

void read_ctl(int argc,
	      char *argv[],
	      ctl_t *ctl) {
  
  int id, ig, iw;
  
  /* Write info... */
  printf("\nJURASSIC - Version with Scattering\n"
	 "Executable: %s\nControl parameters: %s\n", argv[0], argv[1]);
  
  /* Emitters... */
  ctl->ng=(int)jur_scan_ctl(argc, argv, "NG", -1, "0", NULL);
  if(ctl->ng<0 || ctl->ng>NGMAX)
    ERRMSG("Set 0 <= NG <= MAX!");
  for(ig=0; ig<ctl->ng; ig++)
    jur_scan_ctl(argc, argv, "EMITTER", ig, "", ctl->emitter[ig]);
  
  /* Radiance channels... */
  ctl->nd=(int)jur_scan_ctl(argc, argv, "ND", -1, "0", NULL);
  if(ctl->nd<0 || ctl->nd>NDMAX)
    ERRMSG("Set 0 <= ND <= MAX!");
  for(id=0; id<ctl->nd; id++)
    ctl->nu[id]=jur_scan_ctl(argc, argv, "NU", id, "", NULL);

  /* Spectral windows... */
  ctl->nw=(int)jur_scan_ctl(argc, argv, "NW", -1, "1", NULL);
  if(ctl->nw<0 || ctl->nw>NWMAX)
    ERRMSG("Set 0 <= NW <= NWMAX!");
  for(id=0; id<ctl->nd; id++)
    ctl->window[id]=(int)jur_scan_ctl(argc, argv, "WINDOW", id, "0", NULL);
  
  /* Emissivity look-up tables... */
  jur_scan_ctl(argc, argv, "TBLBASE", -1, "-", ctl->tblbase);
  
  /* Hydrostatic equilibrium... */
  ctl->hydz=jur_scan_ctl(argc, argv, "HYDZ", -1, "-999", NULL);
  
  /* Continua... */
  ctl->ctm_co2=(int)jur_scan_ctl(argc, argv, "CTM_CO2", -1, "1", NULL);
  ctl->ctm_h2o=(int)jur_scan_ctl(argc, argv, "CTM_H2O", -1, "1", NULL);
  ctl->ctm_n2=(int)jur_scan_ctl(argc, argv, "CTM_N2", -1, "1", NULL);
  ctl->ctm_o2=(int)jur_scan_ctl(argc, argv, "CTM_O2", -1, "1", NULL);
  
  /* Scattering on Aerosol/Clouds ... */
  ctl->sca_n=(int)jur_scan_ctl(argc, argv, "SCA_N", -1, "0", NULL);
  if(ctl->sca_n<0 || ctl->sca_n>SCAMOD)
    ERRMSG("Set 0 <= SCA_NMOD <= MAX!");
  ctl->sca_mult=(int)jur_scan_ctl(argc, argv, "SCA_MULT", -1, "1", NULL);
  jur_scan_ctl(argc, argv, "SCA_EXT", -1, "beta_a", ctl->sca_ext);
  if(ctl->sca_n>0 && ctl->sca_mult==0 && 
     !strcasecmp(ctl->sca_ext, "beta_e") && 
     !strcasecmp(ctl->sca_ext, "beta_a"))
    ERRMSG("Please set extinction to beta_a or beta_e.");

  /* Interpolation of atmospheric data... */
  ctl->ip=(int)jur_scan_ctl(argc, argv, "IP", -1, "1", NULL);
  ctl->cz=jur_scan_ctl(argc, argv, "CZ", -1, "0", NULL);
  ctl->cx=jur_scan_ctl(argc, argv, "CX", -1, "0", NULL);
  
  /* Ray-tracing... */
  ctl->refrac=(int)jur_scan_ctl(argc, argv, "REFRAC", -1, "1", NULL);
  ctl->rayds=jur_scan_ctl(argc, argv, "RAYDS", -1, "10", NULL);
  ctl->raydz=jur_scan_ctl(argc, argv, "RAYDZ", -1, "1", NULL);
  ctl->transs=jur_scan_ctl(argc, argv, "TRANSS", -1, "0.02", NULL);
  
  /* Field of view... */
  jur_scan_ctl(argc, argv, "FOV", -1, "-", ctl->fov);
  
  /* Retrieval interface... */
  ctl->retp_zmin=jur_scan_ctl(argc, argv, "RETP_ZMIN", -1, "-999", NULL);
  ctl->retp_zmax=jur_scan_ctl(argc, argv, "RETP_ZMAX", -1, "-999", NULL);
  ctl->rett_zmin=jur_scan_ctl(argc, argv, "RETT_ZMIN", -1, "-999", NULL);
  ctl->rett_zmax=jur_scan_ctl(argc, argv, "RETT_ZMAX", -1, "-999", NULL);
  for(ig=0; ig<ctl->ng; ig++) {
    ctl->retq_zmin[ig]=jur_scan_ctl(argc, argv, "RETQ_ZMIN", ig, "-999", NULL);
    ctl->retq_zmax[ig]=jur_scan_ctl(argc, argv, "RETQ_ZMAX", ig, "-999", NULL);
  }
  for(iw=0; iw<ctl->nw; iw++) {
    ctl->retk_zmin[iw]=jur_scan_ctl(argc, argv, "RETK_ZMIN", iw, "-999", NULL);
    ctl->retk_zmax[iw]=jur_scan_ctl(argc, argv, "RETK_ZMAX", iw, "-999", NULL);
  }
  ctl->retnn=(int)jur_scan_ctl(argc, argv, "RETNN", -1, "0", NULL);
  ctl->retrr=(int)jur_scan_ctl(argc, argv, "RETRR", -1, "0", NULL);
  ctl->retss=(int)jur_scan_ctl(argc, argv, "RETSS", -1, "0", NULL);

  ctl->retnn_zmin=(int)jur_scan_ctl(argc, argv, "RETNN_ZMIN", -1, "-999", NULL);
  ctl->retnn_zmax=(int)jur_scan_ctl(argc, argv, "RETNN_ZMAX", -1, "-999", NULL);

  /* Output flags... */
  ctl->write_bbt=(int)jur_scan_ctl(argc, argv, "WRITE_BBT", -1, "0", NULL);
  ctl->write_matrix=(int)jur_scan_ctl(argc, argv, "WRITE_MATRIX", -1, "0", NULL);
  
  //Added:
  /* Number of leaf rays ... */
  ctl->leaf_nr=(int)jur_scan_ctl(argc, argv, "MAX_QUEUE", -1, "0", NULL);
  
  //Added: useGPU
  ctl->useGPU=(int)jur_scan_ctl(argc, argv, "USEGPU", -1, "0", NULL);
}

/*****************************************************************************/
