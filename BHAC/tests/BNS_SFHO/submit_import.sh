#! /bin/bash

#SBATCH -J BHACp
#SBATCH --nodes=8 #3
#SBATCH --ntasks-per-node=30 #10
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --exclusive
#SBATCH --mem=0 #250G

#SBATCH --time=300:00:00
#SBATCH --error="slurm_run.err"
#SBATCH --output="slurm_run.out"
#SBATCH --partition=calea
#SBATCH --constraint=infiniband

starttime=$(date +"%T")
echo "Current time : $starttime"
start=`date +%s`

source ~/.bhac_env.sh
#module restore FUKA_GNU
#module restore fuka_gnu_16_02_2025 

#module purge
#module restore BHAC_MODULES

# Increase limits
#ulimit -n 65535
ulimit -n 8192
ulimit -a  | grep files 

#srun --mpi=pmi2 -n 20 ./bhac -i amrvac.par
#srun --mpi=pmix -n 20 ./bhac -i amrvac.par
mpirun -np 240 ./bhac -i import.par
# amrvac_foo4.par #amrvac_m1.par #amrvac_m1_1var.par #amrvac_m1.par



#export SLURM_MPI_TYPE=pmi
#export PMI_NO_FORK=1
#srun --mpi=pmi -n 20 ./bhac -i amrvac.par

finaltime=$(date +"%T")
echo "Finishing, current time : $finaltime"


end=`date +%s`
runtime=$((end-start))


echo "Total runtime was: $runtime"

