#!/bin/bash -x

#SLURM JOB INFORMATION GLOBALS
export CLUSTER_ACCOUNT=your_account
export PARTITION=your_partition
# SLURM JOB SETUP PARAMETERS
export SETUP_MEM=1GB
export SETUP_TIME=10 # in minutes
# PATHS
# PATHS.CONTAINER
export APPTAINER_DEF_PATH="/"
export EXPRIMENTS_PATH_BASE="/"
export PARAM_FILE_PATH="$EXPERIMENTS_PATH/paramFiles"
export DATA_OUT_FILE_PATH="$EXPERIMENTS_PATH/Data"
export TEMP=""
export MODEL_FILES="$TEMP/Models"
# Paths ExperimentalParameters
export EXP_PARAM_SCIP_PATH="$PARAM_FILE_PATH/scip_experimental_settings.toml"
export EXP_PARAM_SCIPY_PATH="$PARAM_FILE_PATH/SciPy_experimental_settings.toml"
export EXP_PARAM_PATH="$PARAM_FILE_PATH/experimental_parameters.toml"
#Trial SLURM parameters
export TRIAL_MEM=1GB
export TRIAL_TIME=10 # in minutes
