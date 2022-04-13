#!/bin/sh
#
#
# Replace <ACCOUNT> with your account name before submitting.
#
#SBATCH --account=astro      # The account name for the job.
#SBATCH --job-name=I     # The job name.
#SBATCH -o "/rigel/astro/users/cx2204/stars/nuc_burning/Execute-%j.out" 
#SBATCH --mail-type  ALL
#SBATCH --mail-user  cx2204@columbia.edu
#SBATCH -c 10                     # The number of cpu cores to use.
#SBATCH --time=200:00              # The time the job will take to run (here, 60 min)
#SBATCH --mem-per-cpu=1gb        # The memory the job will use per cpu core.

./rn
