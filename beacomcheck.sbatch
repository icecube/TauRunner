#!/bin/bash
# beacomcheck.sbatch
#
#SBATCH -J beacomcheck # A single job name for the array
#SBATCH -p arguelles_delgado # Partition
#SBATCH -n 3 # three cores
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # Running time of 10 hours
#SBATCH --mem 4000 # Memory request of 4 GB
#SBATCH -o trash/beacomcheck_%A_%a.out # Standard output
#SBATCH -e trash/beacomcheck_%A_%a.err # Standard error
cd TauRunner/python
python main.py -n 1000 -t $SLURM_ARRAY_TASK_ID -s 1 -e 100000 -secondaries True -save ~/TauRunner/python/beacomcheck
