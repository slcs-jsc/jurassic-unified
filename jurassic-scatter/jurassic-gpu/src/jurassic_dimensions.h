#ifndef JURASSIC_DIMENSIONS_H
#define JURASSIC_DIMENSIONS_H

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/* Maximum length of ASCII data lines. */
#define LENMAX 10000

/* Maximum size of measurement vector. */
#define MMAX (NRMAX*NDMAX)

/* Maximum size of state vector. */
#define NMAX (NQMAX*NPMAX)

/* Maximum number of quantities. */
#define NQMAX (2+NGMAX+NWMAX+20)

/* Maximum number of spectral windows. */
#define NWMAX 5

/* Maximum number of radiance channels. */
#define NDMAX 70

/* Maximum number of emitters. */
#define NGMAX 23

/* Maximum number of LOS points. */
#define NLOSMAX 10000

/* Maximum number of atmospheric data points. */
#define NPMAX 1000

/* Maximum number of ray paths. */
#define NRMAX 1000

/* Maximum number of shape function grid points. */
#define NSHAPEMAX 10000

/* Number of ray paths used for FOV calculations. */
#define NFOVMAX 50

/* Maximum number of pressure levels in emissivity tables. */
#define TBLNPMAX 45

/* Maximum number of source function temperature levels. */
#define TBLNSMAX 1201

/* Maximum number of temperatures in emissivity tables. */
#define TBLNTMAX 30

/* Maximum number of column densities in emissivity tables. */
#define TBLNUMAX 430

/* Maximum number of scattering models. */
#define SCAMODMAX 30

/* Maximum number of aerosol/cloud layers */
#define NLMAX 10

/* Number of scattering angles (from 0 to 180 deg). */
#define NTHETAMAX 181

/* Number of points for Gauss-Hermite integration. */
#define NRADMAX 170

/* Maximum number of refractive indices. */
#define REFMAX 5000

/* Maximum number of RFM spectral grid points. */
#define RFMNPTSMAX 10000000

/* Maximum length of RFM data lines. */
#define RFMLINEMAX 100000

#endif
