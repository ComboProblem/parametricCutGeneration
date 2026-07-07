#! /bin/bash -x

source src/Experiments/source/cluster_enviroment.sh

module load apptainer
apptainer run $OPTIMAL_CUT_CONTAINER "./src/Experiments/source/setup_experiments.py"
