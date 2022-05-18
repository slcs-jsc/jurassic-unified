#ifndef JURASSIC_DIMENSIONS_H
#define JURASSIC_DIMENSIONS_H

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/* Maximum length of ASCII data lines. */
#define LEN 5000 //Ok!

/* Maximum size of measurement vector. */
#define MMAX (NRMAX*NDMAX) //TODO: M (NR*ND)
#define M MMAX

/* Maximum size of state vector. */
#define NMAX (NQMAX*NPMAX) //TODO: N (NQ*NP)
#define N NMAX

/* Maximum number of quantities. */
#define NQMAX (2+NGMAX+NWMAX) //TODO: NQ (2+NG+NW)
#define NQ NQMAX

/* Maximum number of spectral windows. */
#define NWMAX 5 //TODO: NW 1 
#define NW NWMAX

/* Maximum number of radiance channels. */
#define NDMAX 70 //TODO: ifndef ND 100
#ifndef ND
  #define ND NDMAX 
#endif

/* Maximum number of emitters. */
#define NGMAX 23 //TODO: ifndef NG 30
#ifndef NG
  #define NG NGMAX
#endif  

/* Maximum number of LOS points. */
#define NLOS 10000 //TODO: NLOS 400

/* Maximum number of atmospheric data points. */
#define NPMAX 1000 //TODO: NP 9600
#define NP NPMAX

/* Maximum number of ray paths. */
#define NRMAX 1000 //TODO: NR 1088
#define NR NRMAX

/* Maximum number of shape function grid points. */
#define NSHAPE 10000 //TODO: NSHAPE 2048 | TODO: ?

/* Number of ray paths used for FOV calculations. */
#define NFOV 50 //TODO: NFOV 5

/* Maximum number of pressure levels in emissivity tables. */
#define TBLNPMAX 45 //TODO: TBLNP 40
#define TBLNP TBLNPMAX

/* Maximum number of source function temperature levels. */
#define TBLNSMAX 1201 //TODO: TBLNS 1201
#define TBLNS TBLNSMAX

/* Maximum number of temperatures in emissivity tables. */
#define TBLNTMAX 30 //TODO: TBLNT 30
#define TBLNT TBLNTMAX

/* Maximum number of column densities in emissivity tables. */
#define TBLNUMAX 430 //TODO: TBLNU 304
#define TBLNU TBLNUMAX

/* Maximum number of scattering models. */
#define SCAMOD 30

/* Maximum number of aerosol/cloud layers */
#define NLMAX 10

/* Number of scattering angles (from 0 to 180 deg). */
#define NTHETA 181

/* Number of points for Gauss-Hermite integration. */
#define NRAD 170

/* Maximum number of refractive indices. */
#define REFMAX 5000

//Added:

/* Maximum number of RFM spectral grid points. */
#define RFMNPTS 10000000

/* Maximum length of RFM data lines. */
#define RFMLINE 100000

#endif
