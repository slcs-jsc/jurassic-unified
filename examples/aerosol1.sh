#!/bin/bash -x
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=12
#SBATCH --output=out
#SBATCH --error=err
#SBATCH --time=00:05:00
#SBATCH --partition=booster
#SBATCH --gres=gpu:4
#SBATCH --account=slmet

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

info "Create atmosphere..."
$src/module_climatology aerosol1.ctl - atm.tab || exit

info "Create geometry..."
$src/module_limb aerosol1.ctl 800 5 15 1 obs.tab || exit

info "Call forward model..."
$src/sca_formod aerosol1.ctl obs.tab atm.tab rad_aero1.tab AEROFILE aero1.tab|| exit

# info "Compare files..."
compare rad_aero1.tab

