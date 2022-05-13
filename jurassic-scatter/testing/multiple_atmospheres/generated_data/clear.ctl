# Emitters...
NG = 3
EMITTER[0] = C2H2
EMITTER[1] = C2H6
EMITTER[2] = CCl4

# Channels...
ND = 40
NU[0] = 785.0000
NU[1] = 785.4200
NU[2] = 785.8400
NU[3] = 786.2600
NU[4] = 786.6800
NU[5] = 787.1000
NU[6] = 787.5200
NU[7] = 787.9400
NU[8] = 788.3600
NU[9] = 788.7800
NU[10] = 789.2000
NU[11] = 789.6200
NU[12] = 790.0400
NU[13] = 790.4600
NU[14] = 790.8800
NU[15] = 791.3000
NU[16] = 791.7200
NU[17] = 792.1400
NU[18] = 792.5600
NU[19] = 792.9800
NU[20] = 793.4000
NU[21] = 793.8200
NU[22] = 794.2400
NU[23] = 794.6600
NU[24] = 795.0800
NU[25] = 795.5000
NU[26] = 795.9200
NU[27] = 796.3400
NU[28] = 796.7600
NU[29] = 797.1800
NU[30] = 797.6000
NU[31] = 798.0200
NU[32] = 798.4400
NU[33] = 798.8600
NU[34] = 799.2800
NU[35] = 799.7000
NU[36] = 800.1200
NU[37] = 800.5400
NU[38] = 800.9600
NU[39] = 801.3800

TBLBASE = /p/fastdata/slmet/slmet111/model_data/jurassic/tab/crista-nf/ascii_785_964/crista_nf 
# TBLBASE = /p/fastdata/slmet/slmet111/model_data/jurassic/tab/crista-nf/binary_785_964/crista_nf

# Continua...
CTM_CO2 = 1
CTM_H2O = 1
CTM_N2 = 0
CTM_O2 = 0

# Aerosol/Clouds...
SCA_N = 0
SCA_MULT = 0
SCA_EXT = beta_e

# Raytracing...
RAYDS = 10.
RAYDZ = 0.1

# Atmosphere/Climatology...
ZMIN = 0
ZMAX = 80
DZ = 1

READ_BINARY = 0
WRITE_BINARY = 0

MAX_QUEUE = 1e5
USEGPU = 0
