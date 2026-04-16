#!/bin/bash -x

#SLURM JOB INFORMATION GLOBALS
export CLUSTER_ACCOUNT=your_account
export PARTITION=your_partition
# SLURM JOB SETUP PARAMETERS
export SETUP_MEM=1GB
export SETUP_TIME=10 # in minutes

# PATHS
export APPTAINER_DEF_PATH="src/Experiments/source/Apptainer.def"
export OPTIMAL_CUT_CONTAINER="src/Experiments/source/optimal_cut.sif" 
export PARAMETRIC_EXPS_BASE="src/Experiments/ConductedExperimenets"
export PARAM_FILE_PATH="$EXPRIMENTS_PATH_BASE/paramFiles"
export DATA_OUT_FILE_PATH="$EXPRIMENTS_PATH_BASE/{experiment}/Data"
export TEMP="src/TEMP"
export MODEL_FILES="src/Models"
# Paths ExperimentalParameters
export EXP_PARAM_SCIP_PATH="$PARAM_FILE_PATH/scip_experimental_settings.toml"
export EXP_PARAM_SCIPY_PATH="$PARAM_FILE_PATH/SciPy_experimental_settings.toml"
export EXP_PARAM_PATH="$PARAMETRIC_EXPS_BASE/{experiment}/{experiment}.toml"
#Trial SLURM parameters
export TRIAL_MEM=16GB
export TRIAL_TIME=60 # in minutes
