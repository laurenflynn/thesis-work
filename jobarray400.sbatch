#!/bin/bash
#SBATCH -J 10000families1to400array # A single job name for the array
#SBATCH -p serial_requeue # Partition
#SBATCH -c 1 # one core
#SBATCH -t 0-2:30 # Running time of 2.5 hours
#SBATCH --mem 20000 # Memory request of 20 gb
#SBATCH -o slurm_log_files/jobarray1_400/10000families_%A_%a.out # Standard output
#SBATCH -e slurm_log_files/jobarray1_400/10000families_%A_%a.err # Standard error
#SBATCH --array=138,229
singularity exec --cleanenv --env R_LIBS_USER=$HOME/R/ifxrstudio/RELEASE_3_15 /n/singularity_images/informatics/ifxrstudio/ifxrstudio:RELEASE_3_15.sif Rscript src/2_slurm_with_loops.R ${SLURM_ARRAY_TASK_ID}
