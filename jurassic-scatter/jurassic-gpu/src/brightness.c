/*
  This file is part of JURASSIC.
  
  JURASSIC is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  JURASSIC is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with JURASSIC. If not, see <http://www.gnu.org/licenses/>.
  
  Copright (C) 2003-2015 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Convert radiance to brightness temperature.
*/

#include "jurassic.h"

int main(
  int argc,
  char *argv[]) {

  double nu, rad;

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Give parameters: <rad> <nu>");

  /* Read arguments... */
  rad = atof(argv[1]);
  nu = atof(argv[2]);

  /* Compute brightness temperature... */
  printf("%.10g\n", jur_brightness(rad, nu));

  return EXIT_SUCCESS;
}
