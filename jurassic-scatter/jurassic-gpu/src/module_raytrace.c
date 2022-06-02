#include "jurassic.h"
#include "sca_forwardmodel.h"

int main(int argc, char *argv[]) {

  static atm_t atm;

  static ctl_t ctl;

  static pos_t los[NLOSMAX];

  static obs_t obs;

  static aero_t aero;

  FILE *out;

  char filename[LENMAX], aerofile[LENMAX];

  int id, ig, ip, ir, iw;

  /* Check arguments... */
  if(argc<4)
    ERRMSG("Give parameters: <ctl> <obs> <atm>");

  /* Read control parameters... */
  jur_read_ctl(argc, argv, &ctl);

  /* Get aero... */
  jur_scan_ctl(argc, argv, "AEROFILE", -1, "-", aerofile);

  /* Read observation geometry... */
  jur_read_obs(NULL, argv[2], &ctl, &obs);

  /* Read atmospheric data... */
  jur_read_atm(NULL, argv[3], &ctl, &atm);

  /* Read aerosol and cloud data */
  if(aerofile[0]!='-' && ctl.sca_n>0) {
    jur_sca_read_aero(NULL, aerofile, &ctl, &aero);
    /* Get aerosol/cloud optical properties */
    jur_sca_get_opt_prop(&ctl, &aero);
  } else if (aerofile[0]=='-' && ctl.sca_n>0) {
    ERRMSG("Please give aerosol file name or set SCA_N=0 for clear air simulation!");
  }

  /* Loop over rays... */
  for(ir=0; ir<obs.nr; ir++) {

    /* Raytracing... */
    double tsurf;
    int np = jur_traceray(&ctl, &atm, &obs, ir, los, &tsurf, &aero, ctl.sca_n>0);

    /* Create file... */
    sprintf(filename, "los.%d", ir);
    if(!(out=fopen(filename, "w")))
      ERRMSG("Cannot create los.tab!");

    /* Write header... */
    fprintf(out,
        "# $1 = time (seconds since 2000-01-01T00:00Z)\n"
        "# $2 = LOS point altitude [km]\n"
        "# $3 = LOS point longitude [deg]\n"
        "# $4 = LOS point latitude [deg]\n"
        "# $5 = LOS point pressure [hPa]\n"
        "# $6 = LOS point temperature [K]\n");
    for(ig=0; ig<ctl.ng; ig++)
      fprintf(out, "# $%d = LOS point %s volume mixing ratio \n",
          7+ig, ctl.emitter[ig]);
    for(iw=0; iw<ctl.nw; iw++)
      fprintf(out, "# $%d = LOS point window %d extinction [1/km]\n",
          7+ctl.ng+iw, iw);
    for(ig=0; ig<ctl.ng; ig++)
      fprintf(out, "# $%d = LOS point %s column density [molec/cm^2] \n",
          7+ctl.ng+ctl.nw+ig, ctl.emitter[ig]);
    for(id=0; id<ctl.nd; id++)
      fprintf(out, "# $%d = LOS point beta_e(%.4f cm-1) [km-1] \n", 7+ctl.ng+ctl.nw+ctl.ng+id, ctl.nu[id]);
    for(id=0; id<ctl.nd; id++)
      fprintf(out, "# $%d = LOS point beta_s(%.4f cm-1) [km-1] \n", 7+ctl.ng+ctl.nw+ctl.ng+ctl.nd+id, ctl.nu[id]);
    for(id=0; id<ctl.nd; id++)
      fprintf(out, "# $%d = LOS point beta_a(%.4f cm-1) [km-1] \n", 7+ctl.ng+ctl.nw+ctl.ng+2*ctl.nd+id, ctl.nu[id]);
    fprintf(out, "# $%d = LOS segement length [km] \n", 7+ctl.ng+ctl.nw+ctl.ng+3*ctl.nd);

    fprintf(out, "\n");

    /* Loop over LOS points... */
    for(ip=0; ip<np; ip++) {
      fprintf(out, "%.2f %g %g %g %g %g", obs.time[ir],
          los[ip].z, los[ip].lon, los[ip].lat,
          los[ip].p, los[ip].t);
      for(ig=0; ig<ctl.ng; ig++)
        fprintf(out, " %g", los[ip].q[ig]);
      for(iw=0; iw<ctl.nw; iw++)
        fprintf(out, " %g", los[ip].k[iw]);
      for(ig=0; ig<ctl.ng; ig++)
        fprintf(out, " %g", los[ip].u[ig]);
      for(id=0; id<ctl.nd; id++)
        fprintf(out, " %g", aero.beta_e[los[ip].aeroi][id]*los[ip].aerofac);
      for(id=0; id<ctl.nd; id++)
        fprintf(out, " %g", aero.beta_s[los[ip].aeroi][id]*los[ip].aerofac);
      for(id=0; id<ctl.nd; id++)
        fprintf(out, " %g", aero.beta_a[los[ip].aeroi][id]*los[ip].aerofac);
      fprintf(out, " %g", los[ip].ds);

      fprintf(out, "\n");
    }

    /* Close file... */
    fclose(out);

    printf("Wrote output to %s \n",filename);
  }

  return EXIT_SUCCESS;
}
