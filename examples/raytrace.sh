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
$src/module_climatology clear-air.ctl - atm.tab CLIMZONE pwin|| exit

info "Create observation geometry..."
$src/module_limb clear-air.ctl 800 8 8 1 obs_raytrace.tab || exit

info "Call raytrace module..."
$src/module_raytrace clear-air.ctl obs_raytrace.tab atm.tab || exit

info "Compare files..."
compare atm.tab
compare obs_raytrace.tab
compare los.0
