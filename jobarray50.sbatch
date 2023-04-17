#!/bin/bash
#SBATCH -J 10000families1to50arrayRetryFailedJobs # A single job name for the array
#SBATCH -p serial_requeue # Partition
#SBATCH -c 1 # one core
#SBATCH -t 0-2:30 # Running time of 4.5 hours
#SBATCH --mem 20000 # Memory request of 20 gb
#SBATCH -o slurm_log_files/jobarray1_200/10000families_%A_%a.out # Standard output
#SBATCH -e slurm_log_files/jobarray1_200/10000families_%A_%a.err # Standard error
#SBATCH --array=101,129,130,140,148,165,167,211,213,219,229,231,233,251,262,266,279,286,294,295,297,300,302,318,344,361,365,367,372,386,387,389,397,399,221
singularity exec --cleanenv --env R_LIBS_USER=$HOME/R/ifxrstudio/RELEASE_3_15 /n/singularity_images/informatics/ifxrstudio/ifxrstudio:RELEASE_3_15.sif Rscript 2_slurm_with_loops.R ${SLURM_ARRAY_TASK_ID}
