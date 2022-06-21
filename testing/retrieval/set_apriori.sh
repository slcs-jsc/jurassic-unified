#!/bin/bash

echo "# $1 = aerosol layer top altitude [km]"> aero_apr.tab
echo "# $2 = aerosol layer bottom altitude [km]" >> aero_apr.tab
echo "# $3 = transition layer thickness [km]" >> aero_apr.tab
echo "# $4 = source for optical properties" >> aero_apr.tab
echo "# $5 = refractive index file" >> aero_apr.tab
echo "# $6 = particle concentration of log-normal mode [cm-3]" >> aero_apr.tab
echo "# $7 = median radius of log-normal mode [mum]" >> aero_apr.tab
echo "# $8 = width of log-normal mode" >> aero_apr.tab
echo "18 17 0.01 MIE /private/sgrie/ESA-study/MIPAS/h2so4-215K-shettle.dat 1.617 1.0 1.6" >> aero_apr.tab

for path in h2so4/midln/18km/*
do
    #cp aero_apr.tab $path/
    cp $path/aero.tab $path/aero_apr.tab
    cp $path/atm.tab $path/atm_apr.tab
    cp $path/rad-retrieval.tab $path/obs_meas.tab 
done
rm aero_apr.tab
