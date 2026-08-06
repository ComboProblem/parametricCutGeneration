#!/bin/bash

PARTITION="partition_name"
CLUSTER_ACCOUNT="account_name"


run_experiments(){
for model in model_files/*
do
    model_name=$(basename $model .mps)
    echo "Queueing experiments for $model_name."
    sbatch --array=0-64 --account=$CLUSTER_ACCOUNT --partition=$PARTITION --time=240:00 --output="TEMP/${model_name}_ID_%a.out" source/run_experiments.sh $model_name 
done
}

main() {
chmod +x source/run_experiments.sh
chmod +x setup.sh
setup
mkdir TEMP
run_experiments
}

main

