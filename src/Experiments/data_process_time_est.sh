#!/bin/bash

PARTITION="partition_name"
CLUSTER_ACCOUNT="account_name"


process_data_timing(){
JOB_COUNT=0
for model in model_files/*
do
    model_name=$(basename $model .mps)
    echo "Queing data processing pipeline for $model_name."
    chmod +x "./TEMP/{$model_name}_data_run.sh"
    ./TEMP/{$model_name}_data_run.sh $CLUSTER_ACCOUNT $PARTITON
    ((JOB_COUNT=JOB_COUNT+1))
    if (( JOB_COUNT == 80 )); then
        echo "Giving some time for other jobs to finish."
        sleep 5m
        JOB_COUNT=0
    fi
done
}

main() {
chmod +x source/estimate_process_time.sh
#mkdir TEMP
process_data_timing

}

main

