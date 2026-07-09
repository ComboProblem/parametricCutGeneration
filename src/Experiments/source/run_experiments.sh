#!/bin/bash

export MODEL_NAME="$1"
echo "Running on $(hostname)"
echo "Starting Trial: Model: $MODEL_NAME ExperimentID: $SLURM_ARRAY_TASK_ID"

module load apptainer
apptainer run container/parametricCutGen.sif python3 source/run_trial.py

echo "Finished ExperimentID: $SLURM_ARRAY_TASK_ID for $MODEL_NAME"

