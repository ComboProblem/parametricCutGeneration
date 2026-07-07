#!/bin/bash

export PROBLEM = "$1"

echo "Running on $(hostname)"
echo "Task: $SLURM_ARRAY_TASK_ID"
apptainer run container/optimal_cut.sif python3 run_trial.py
