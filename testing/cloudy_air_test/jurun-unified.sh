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

jurassic_formod=../../src/sca_formod

$jurassic_formod cloudy-785-812-ascii.ctl obs.tab atm.tab rad.tab AEROFILE aero.tab
