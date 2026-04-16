#! /bin/bash -x

source 

module load apptainer
apptainer run "./src/Experiments/source/optimal_cut.sif" "./src/Experiments/source/setup_experiments.py"
