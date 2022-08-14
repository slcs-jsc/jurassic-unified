#ifndef JURASSIC_STRUCTS_H
#define JURASSIC_STRUCTS_H

/* ------------------------------------------------------------
   Structs...
   ------------------------------------------------------------ */
typedef struct {
  void* result;
  int ir;
} queue_item_t;

typedef struct {
  queue_item_t* items;
  int capacity;
  int begin;
  int end;
} queue_t;

typedef struct { /// Atmospheric data. /////////////////////////////////////////
  double time[NPMAX];      /// Time (seconds since 2000-01-01T00:00Z).
  double z[NPMAX];           /// Altitude [km].
  double lon[NPMAX];         /// Longitude [deg].
  double lat[NPMAX];         /// Latitude [deg].
  double p[NPMAX];           /// Pressure [hPa].
  double t[NPMAX];           /// Temperature [K].
  double q[NGMAX][NPMAX];      /// Volume mixing ratio.
  double k[NWMAX][NPMAX];      /// Extinction [1/km].
  int np;                /// Number of data points.
  int init;              /// Init flag for interpolation (0=no, 1=yes).
} atm_t; ///////////////////////////////////////////////////////////////////////

typedef struct { /// Forward model control parameters. /////////////////////////
  int ng;                 /// Number of emitters.
  char emitter[NGMAX][LENMAX];  /// Name of each emitter.
  int nd;                 /// Number of radiance channels.
  int nw;                 /// Number of spectral windows.
  double nu[NDMAX];          /// Centroid wavenumber of each channel [cm^-1].
  int window[NDMAX];         /// Window index of each channel.
  char tblbase[LENMAX];      /// Basename for table files and filter function files.
  double hydz;            /// Reference height for hydrostatic pressure profile (-999 to skip) [km].
  int ctm_co2;            /// Compute CO2 continuum (0=no, 1=yes).
  int ctm_h2o;            /// Compute H2O continuum (0=no, 1=yes).
  int ctm_n2;             /// Compute N2 continuum (0=no, 1=yes).
  int ctm_o2;             /// Compute O2 continuum (0=no, 1=yes).
  int ip;                 /// Interpolation method (1=profile, 2=satellite track, 3=Lagrangian grid).
  double cz;              /// Influence length for vertical interpolation [km].
  double cx;              /// Influence length for horizontal interpolation [km].
  int refrac;             /// Take into account refractivity (0=no, 1=yes).
  double rayds;           /// Maximum step length for raytracing [km].
  double raydz;           /// Vertical step length for raytracing [km].
  char fov[LENMAX];          /// Field-of-view data file.
  double retp_zmin;       /// Minimum altitude for pressure retrieval [km].
  double retp_zmax;       /// Maximum altitude for pressure retrieval [km].
  double rett_zmin;       /// Minimum altitude for temperature retrieval [km].
  double rett_zmax;       /// Maximum altitude for temperature retrieval [km].
  double retq_zmin[NGMAX];   /// Minimum altitude for volume mixing ratio retrieval [km].
  double retq_zmax[NGMAX];   /// Maximum altitude for volume mixing ratio retrieval [km].
  double retk_zmin[NWMAX];   /// Minimum altitude for extinction retrieval [km].
  double retk_zmax[NWMAX];   /// Maximum altitude for extinction retrieval [km].
  int write_bbt;          /// Use brightness temperature instead of radiance (0=no, 1=yes).
  int write_matrix;       /// Write matrix file (0=no, 1=yes).
  int formod;             /// Forward model (1=CGA, 2=EGA, 3=RFM).
  char rfmbin[LENMAX];       /// Path to RFM binary.
  char rfmhit[LENMAX];       /// HITRAN file for RFM.
  char rfmxsc[NGMAX][LENMAX];   /// Emitter cross-section files for RFM.
  int useGPU;             /// Use GPU-accelerated formod implementation (0=no, 1=yes)
  int checkmode;          /// do not perform input, computation, nor output, just make sure files are there
  int MPIglobrank;        /// MPI global rank
  int MPIlocalrank;       /// MPI node-local Rank
  int gpu_nbytes_shared_memory; /// Shared memory controler for GPU kernels

  /// ---------------- for scattering ------------------
  int sca_n;              /// Number of scattering models.
  int sca_mult;           /// Number of recursions for multiple scattering.
  /// (0=no scattering, 1=single scattering, 2<=multiple scattering)
  char sca_ext[LENMAX];      /// Extinction coefficient type if sca_mult=0
  double transs;          /// Sampling step for transition layers [km].
  double retnn_zmin;      /// Minimum altitude for particle [km].
  double retnn_zmax;      /// Maximum altitude for particle retrieval [km].
  int retnn;              /// Retrieval of particle concentration (0=no, 1=yes)
  int retrr;              /// Retrieval of particle size (0=no, 1=yes)
  int retss;              /// Retrieval of particle size distribution width (0=no, 1=yes)
  int leaf_nr;            /// Number of leaf rays, for example number of secondary rays if sca_mult=1
  int queue_state;        /// We have multiple queues, but all of them are in same state
} ctl_t; ///////////////////////////////////////////////////////////////////////

typedef struct {    /// Point on the Line-of-sight data without storing //////////
  double z;         /// Altitude [km].
  double lon;       /// Longitude [deg].
  double lat;       /// Latitude [deg].
  double p;         /// Pressure [hPa].
  double t;         /// Temperature [K].
  double q[NGMAX];      /// Volume mixing ratio.
  double k[NWMAX];      /// Extinction [1/km].
  int aeroi;        /// Aerosol/cloud layer index
  double aerofac;   /// Aerosol/cloud layer scaling factor for transition layer
  double ds;        /// Segment length [km].
  double u[NGMAX];      /// Column density [molecules/cm^2].
#ifdef CURTIS_GODSON
  double cgp[NGMAX];    /// Curtis-Godson pressure [hPa].
  double cgt[NGMAX];    /// Curtis-Godson temperature [K].
  double cgu[NGMAX];    /// Curtis-Godson column density [molecules/cm^2].
#endif
#ifdef GPUDEBUG
  int ip, ir;       /// debug helpers
#endif
} pos_t; //////////////////////////////////////////////////////////////////////

typedef struct { /// Observation geometry and radiance data. //////////////////
  double time[NRMAX];   /// Time (seconds since 2000-01-01T00:00Z).
  double obsz[NRMAX];   /// Observer altitude [km].
  double obslon[NRMAX]; /// Observer longitude [deg].
  double obslat[NRMAX]; /// Observer latitude [deg].
  double vpz[NRMAX];      /// View point altitude [km].
  double vplon[NRMAX];    /// View point longitude [deg].
  double vplat[NRMAX];    /// View point latitude [deg].
  double tpz[NRMAX];      /// Tangent point altitude [km].
  double tplon[NRMAX];    /// Tangent point longitude [deg].
  double tplat[NRMAX];    /// Tangent point latitude [deg].
  double tau[NRMAX][NDMAX]; /// Transmittance of ray path.    // transposed
  double rad[NRMAX][NDMAX]; /// Radiance [W/(m^2 sr cm^-1)].  // transposed
  int nr;             /// Number of ray paths.
} obs_t; ///////////////////////////////////////////////////////////////////////

typedef float real_tblND_t;

typedef struct {  /// Transposed emissivity look-up tables. - GPU version  /////
  int32_t np[NGMAX][NDMAX];                             /// Number of pressure levels.
  int32_t nt[NGMAX][TBLNPMAX][NDMAX];                      /// Number of temperatures.
  int32_t nu[NGMAX][TBLNPMAX][TBLNTMAX][NDMAX];               /// Number of column densities.
  double p[NGMAX][TBLNPMAX][NDMAX];                        /// Pressure [hPa].
  double t[NGMAX][TBLNPMAX][TBLNTMAX][NDMAX];                 /// Temperature [K].
  real_tblND_t u[NGMAX][TBLNPMAX][TBLNTMAX][TBLNUMAX][NDMAX];    /// Column density [molecules/cm^2].
  real_tblND_t eps[NGMAX][TBLNPMAX][TBLNTMAX][TBLNUMAX][NDMAX];  /// Emissivity.
  double sr[TBLNSMAX][NDMAX];                           /// Source function radiance [W/(m^2 sr cm^-1)].
  double st[TBLNSMAX];                               /// Source function temperature [K].
#ifdef  FAST_INVERSE_OF_U
  /// u0inv[g][p][t][d] * u[g][p][t][0][d] == 1 must hold!  /// FAST_INVERSE_OF_U
  double u0inv[NGMAX][TBLNPMAX][TBLNTMAX][NDMAX];                       /// FAST_INVERSE_OF_U
  /// We assume a logarithmic increment by 2^(1/6)          /// FAST_INVERSE_OF_U
#endif
} trans_table_t; ///////////////////////////////////////////////////////////////////////

/// --------------------------------- for scattering -----------------------------------
typedef struct { /// Aerosol and Cloud properties. /////////////////////////////////////
  /// Aerosol and cloud input parameters
  int nm;                         /// Number of aerosol/cloud models
  double top_mod[SCAMODMAX];         /// Model top altitude [km]
  double bottom_mod[SCAMODMAX];      /// Model bottom altitude [km]
  double trans_mod[SCAMODMAX];       /// Model transition layer thickness [km]
  char type[SCAMODMAX][LENMAX];         /// Optical properties source
  char filepath[SCAMODMAX][LENMAX];     /// Refractive index file or optical properties file
  double nn[SCAMODMAX];              /// Number concentration [cm-3] or extinction coefficient [km-1]
  double rr[SCAMODMAX];              /// Median radius of log-normal size distribution [mum]
  double ss[SCAMODMAX];              /// Width of log-normal size distribution
  /// Aerosol and cloud optical properties for radiative transfer
  int nl;                         /// Number of aerosol/cloud layers
  int nmod[NLMAX];                /// Number of modes per layer
  double top[NLMAX];              /// Layer top altitude [km]
  double bottom[NLMAX];           /// Layer bottom altitude [km]
  double trans[NLMAX];            /// Transition layer thickness [km]
  double beta_e[NLMAX][NDMAX];    /// Extinction coefficient [1/km]
  double beta_s[NLMAX][NDMAX];    /// Scattering coefficient [1/km]
  double beta_a[NLMAX][NDMAX];    /// Absorption coefficient [1/km]
  double p[NLMAX][NDMAX][NTHETAMAX]; /// Phase function for each layer, angle and wave number
} aero_t; /////////////////////////////////////////////////////////////////////////////

typedef struct { /// Retrieval control parameters. /////////////////////////////////////
  char dir[LENMAX];            /// Working directory.
  int kernel_recomp;        /// Recomputation of kernel matrix (number of iterations).
  int conv_itmax;           /// Maximum number of iterations.
  double conv_dmin;         /// Minimum normalized step size in state space.
  double resmax;            /// Threshold for radiance residuals [%] (-999 to skip filtering).
  int err_ana;              /// Carry out error analysis (0=no, 1=yes).
  double err_formod[NDMAX]; /// Forward model error [%].
  double err_noise[NDMAX];  /// Noise error [W/(m^2 sr cm^-1)].
  double err_press;         /// Pressure error [%].
  double err_press_cz;      /// Vertical correlation length for pressure error [km].
  double err_press_ch;      /// Horizontal correlation length for pressure error [km].
  double err_temp;          /// Temperature error [K].
  double err_temp_cz;       /// Vertical correlation length for temperature error [km].
  double err_temp_ch;       /// Horizontal correlation length for temperature error [km].
  double err_nn;            /// Particle concentration error [cm-3].
  double err_rr;            /// Particle radius error [m-6].
  double err_ss;            /// Particle size distribution width error.
  double err_q[NGMAX];      /// Volume mixing ratio error [%].
  double err_q_cz[NGMAX];   /// Vertical correlation length for volume mixing ratio error [km].
  double err_q_ch[NGMAX];   /// Horizontal correlation length for volume mixing ratio error [km].
  double err_k[NWMAX];      /// Extinction error [1/km].
  double err_k_cz[NWMAX];   /// Vertical correlation length for extinction error [km].
  double err_k_ch[NWMAX];   /// Horizontal correlation length for extinction error [km].
} ret_t;

#endif
