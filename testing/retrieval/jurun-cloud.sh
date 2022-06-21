#!/bin/bash -x
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --output=out
#SBATCH --error=err
#SBATCH --time=00:20:00
#SBATCH --partition=gpus
#SBATCH --gres=gpu:1
#SBATCH --account=slmet

module load Python/3.8.5
module load GCC/9.3.0
module load OpenMPI/4.1.0rc1
module load CUDA/11.3
module load Nsight-Systems/2021.1.1
module load Nsight-Compute/2020.3.0
module load ParaStationMPI/5.4.7-1
module load mpi-settings/CUDA
module load GCCcore/.10.3.0
module load GCCcore/.9.3.0
module load GSL/2.6 

# Setup...
name=$1


############################################################

# echo "Start forward model..."

# time $src/formod cloud-retrieval.ctl obs.tab atm.tab rad-retrieval.tab AEROFILE aero.tab DIRLIST dirlist-ret.asc

echo "Start retrieval with scattering..."
time $name cloud-retrieval.ctl dirlist-ret.asc

newdir=r6
for path in $(cat dirlist-ret.asc)
do
    mkdir $path/$newdir
    cp $path/*_*.tab $path/$newdir
    cp $path/costs.tab $path/$newdir
done
