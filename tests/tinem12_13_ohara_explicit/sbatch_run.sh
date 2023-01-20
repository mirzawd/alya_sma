#!/bin/bash
#
#  Submit jobs in MN-IV
#     sbatch < job.sh
#
#SBATCH --job-name=atrium
#SBATCH -D .
#SBATCH --error=%j.err
#SBATCH --output=%j.out
#SBATCH --ntasks=96
#SBATCH --tasks-per-node=48
#SBATCH --time=1:00:00
#SBATCH --qos=debug
#
#
# Launches ALYA
#
echo '--| ALYA: JOB SUBMITTED WITH SRUN AT: ' `date`
time srun /home/bsc21/bsc21371/Software/alya/build/src/alya/alya tine_slab
echo '--| ALYA: JOB FINISHED AT: ' `date`

