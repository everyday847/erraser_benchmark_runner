
sbatch_header = """#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=erraser_{the_pdb}_{dirnum}
#################
#a file for job output, you can check job progress
#SBATCH --output=erraser_{the_pdb}_{dirnum}.log
#################
# a file for errors from the job
#SBATCH --error=erraser_{the_pdb}_{dirnum}.err
#################
#time you think you need; default is 2 hours
#format could be dd-hh:mm:ss, hh:mm:ss, mm:ss, or mm
#SBATCH --time=48:00:00
#################
#quality of service; think of it as job priority
#SBATCH --qos=normal
#SBATCH -p biochem,normal,owners
#################
#number of nodes you are requesting
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#################
#memory per node; default is 4000 MB per CPU
# Consider more for big maps??
#SBATCH --mem-per-cpu=8000
#################
 
cd {pwd}\n\n"""
