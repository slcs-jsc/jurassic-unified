
# Emitters...
NG = 2
EMITTER[0] = CCl4
EMITTER[1] = ClO

# Channels...
ND = 2
NU[0] = 785.00
NU[1] = 785.42

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
DZ = 10
CLIMPATH = /p/project/chwu36/hwu361/crista-nf/setup/atmo/clim_pscs.tab

# use the GPU: 0:never, 1:always, -1:if possible
USEGPU = -1
