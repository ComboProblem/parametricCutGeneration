#!/bin/bash -x

#SLURM JOB INFORMATION
export CLUSTER_ACCOUNT=math-grp
export PARTITION=high
export MEM=4gb
export SETUP_TIME=10 # in minutes
export MAX_TIME_PER_EXP=60  # in minutes
export OVERHEAD_TIME=5 # in minutes
#PATHS
export APPTAINER_DEF_PATH="~/MinimalFunctionCache/src/minimalFunctionCache/minimalFunctionCacheGen/Apptainer.def"
export EXPRIMENTS_PATH="~/MinimalFunctionCache/src/parametricCutGen/Experiments/"
export PARAM_FILE_PATH="$EXPERIMENTS_PATH/paramFiles"
export DATA_OUT_FILE_PATH="$EXPERIMENTS_PATH/Data"
export EXP_TEMP="~/MinimalFunctionCache/TEMP"
export MODEL_FILES="$EXP_TEMP/Models"
export TRIAL_FILES_TO_RUN="$EXP_TEMP/Run"
#Experimental Parameters
export EXP_PARAM_SCIP_PATH="$PARAM_FILE_PATH/scip_experimental_settings.toml"
export EXP_PARAM_SCIPY_PATH="$PARAM_FILE_PATH/SciPy_experimental_settings.toml"
export EXP_PARAM_PATH="$PARAM_FILE_PATH/experimental_parameters.toml"

