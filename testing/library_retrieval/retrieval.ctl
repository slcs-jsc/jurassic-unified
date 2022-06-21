# Emitters...
NG = 13
EMITTER[0] = CCl4
EMITTER[1] = ClO
EMITTER[2] = ClONO2
EMITTER[3] = CO2
EMITTER[4] = F11
EMITTER[5] = F113
EMITTER[6] = F22
EMITTER[7] = H2O
EMITTER[8] = HNO3
EMITTER[9] = HNO4
EMITTER[10] = NO2
EMITTER[11] = O3
EMITTER[12] = PAN

# Channels...
ND = 5
NU[0] = 785.00
NU[1] = 785.42
NU[2] = 785.84
NU[3] = 786.26
NU[4] = 786.68
# NU[5] = 787.10
# NU[6] = 787.52
# NU[7] = 787.94
# NU[8] = 788.36
# NU[9] = 788.78
# NU[10] = 789.20
# NU[11] = 789.62
# NU[12] = 790.04
# NU[13] = 790.46
# NU[14] = 790.88
# NU[15] = 791.30
# NU[16] = 791.72
# NU[17] = 792.14
# NU[18] = 792.56
# NU[19] = 792.98
# NU[20] = 793.40
# NU[21] = 793.82
# NU[22] = 794.24
# NU[23] = 794.66
# NU[24] = 795.08
# NU[25] = 795.50
# NU[26] = 795.92
# NU[27] = 796.34
# NU[28] = 796.76
# NU[29] = 797.18
# NU[30] = 797.60
# NU[31] = 798.02
# NU[32] = 798.44

TBLBASE = /p/fastdata/slmet/slmet111/model_data/jurassic/tab/crista-nf/ascii_785_964/crista_nf 

# Continua...
CTM_CO2 = 1
CTM_H2O = 1
CTM_N2 = 0
CTM_O2 = 0

# Aerosol/Clouds...
SCA_N = 1
SCA_MULT = 1
SCA_EXT = beta_e

MAX_QUEUE = 1e5

# Raytracing...
RAYDS = 10.
RAYDZ = 0.1

# Atmosphere/Climatology...
ZMIN = 0
ZMAX = 80
DZ = 1
CLIMPATH = /p/project/chwu36/hwu361/crista-nf/setup/atmo/clim_pscs.tab

# use the GPU: 0:never, 1:always, -1:if possible
USEGPU = -1

# ======================================================================
# Preprocessor...
# ======================================================================

# Retrieval grid...
NZ = 27
Z[0] = 0
Z[1] = 3
Z[2] = 6
Z[3] = 9
Z[4] = 12
Z[5] = 15
Z[6] = 18
Z[7] = 21
Z[8] = 24
Z[9] = 27
Z[10] = 30
Z[11] = 33
Z[12] = 36
Z[13] = 39
Z[14] = 42
Z[15] = 45
Z[16] = 48
Z[17] = 51
Z[18] = 54
Z[19] = 57
Z[20] = 60
Z[21] = 65
Z[22] = 70
Z[23] = 75
Z[24] = 80
Z[25] = 85
Z[26] = 90

# ======================================================================
# Retrieval...
# ======================================================================

# Iteration control...
CONV_ITMAX = 30
KERNEL_RECOMP = 3

# Temperature a priori...
ERR_TEMP = 20
ERR_TEMP_CZ = 50
ERR_TEMP_CH = 1e10

# Error analysis...
ERR_NOISE[0] = 0.0004998
ERR_NOISE[1] = 0.0004777
ERR_NOISE[2] = 0.0005071
ERR_NOISE[3] = 0.0005084
ERR_NOISE[4] = 0.0004924
ERR_NOISE[5] = 0.0004985
ERR_NOISE[6] = 0.0005169
ERR_NOISE[7] = 0.0005647
ERR_NOISE[8] = 0.0004961
ERR_NOISE[9] = 0.0005059
ERR_NOISE[10] = 0.0005046
ERR_NOISE[11] = 1.498e-06
ERR_NOISE[12] = 1.438e-06
ERR_NOISE[13] = 1.82e-06
ERR_NOISE[14] = 1.508e-06
ERR_NOISE[15] = 1.522e-06
ERR_NOISE[16] = 1.398e-06
ERR_NOISE[17] = 1.524e-06
ERR_NOISE[18] = 1.935e-06
ERR_NOISE[19] = 2.242e-06
ERR_NOISE[20] = 1.447e-06
ERR_NOISE[21] = 1.861e-06
ERR_NOISE[22] = 1.53e-06
ERR_NOISE[23] = 1.461e-06
ERR_NOISE[24] = 1.701e-06
ERR_NOISE[25] = 1.423e-06
ERR_NOISE[26] = 1.891e-06
ERR_NOISE[27] = 1.525e-06
ERR_NOISE[28] = 1.535e-06
ERR_NOISE[29] = 1.78e-06
ERR_NOISE[30] = 1.508e-06
ERR_NOISE[31] = 1.94e-06
ERR_NOISE[32] = 1.539e-06
ERR_NOISE[33] = 1.491e-06

