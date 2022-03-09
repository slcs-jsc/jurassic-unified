#! /bin/bash

# ----------------------------------------------------------------------

function compare {
    
    # Write info...
    echo -n "Compare $1: "
    
    # Compare...
    diff -q $1 org/$1 && echo "OK" || exit
    
    # Remove...
    rm $1
}

# ----------------------------------------------------------------------

function info {
    
    # Write info...
    echo
    echo "=================================================="
    echo $1
    echo "=================================================="
    echo
}

# ----------------------------------------------------------------------

# Setup...
src=../src

#info "Create atmosphere..."
$src/climatology aerosol0.ctl - atm.tab CLIMZONE pwin || exit

#info "Create observation geometry..."
$src/nadir aerosol0.ctl 800 0 10 1 obs_nadir.tab || exit

info "Call forward model..."
$src/formod aerosol0.ctl obs_nadir.tab atm.tab rad_aero0.tab AEROFILE aero0.tab|| exit

info "Compare files..."
compare obs_nadir.tab
compare rad_aero0.tab
