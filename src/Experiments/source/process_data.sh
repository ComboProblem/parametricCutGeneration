#!/bin/bash

export MODEL_NAME="$1"
export EST_RUN="n"

echo "Running on $(hostname)"
echo "Processing Data: Model: $MODEL_NAME"

mkdir TEMP/$MODEL_NAME
module load apptainer
apptainer run container/parametricCutGen.sif python3 source/process_data.py

echo "Finished Processing data for $MODEL_NAME"
echo "Cleaning temp files."
rm -rf TEMP/$MODEL_NAME

