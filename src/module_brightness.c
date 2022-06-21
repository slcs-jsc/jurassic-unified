#include "sca_forwardmodel.h"

int main(int argc, char *argv[]) {

  double nu, rad;
  
  /* Check arguments... */
  if(argc<3)
    ERRMSG("Give parameters: <rad> <nu>");

  /* Read arguments... */
  rad=atof(argv[1]);
  nu=atof(argv[2]);

  /* Compute brightness temperature... */
  printf("%g\n", jur_brightness_core(rad, nu));

  return EXIT_SUCCESS;
}
