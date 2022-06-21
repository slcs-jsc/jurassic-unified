#!/bin/bash -x
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --output=out
#SBATCH --error=err
#SBATCH --time=00:05:00
#SBATCH --partition=gpus
#SBATCH --gres=gpu:1
#SBATCH --account=slen

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

export OMP_NUM_THREADS=24

### set up
src=$1

### run forward model
time srun $src/formod cloud.ctl obs.tab atm.tab rad.tab
