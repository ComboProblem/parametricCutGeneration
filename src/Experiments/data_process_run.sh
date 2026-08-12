#!/bin/bash

PARTITION="partition_name"
CLUSTER_ACCOUNT="account_name"


process_data(){
JOB_COUNT=0
for model in model_files/*
do
    model_name=$(basename $model .mps)
    echo "Queing data processing pipeline for $model_name."
    srun --account=$CLUSTER_ACCOUNT --partition=$PARTITION --mem=5G --time=480:00 --output="TEMP/${model_name}_result_ID_%a.out" source/process_data.sh $model_name
    ((JOB_COUNT=JOB_COUNT+1))
    if (( JOB_COUNT == 80 )); then
        echo "Giving some time for other jobs to finish."
        sleep 4h
        JOB_COUNT=0
    fi
done
}

main() {
chmod +x source/process_data.sh
#mkdir TEMP
process_data
}

main

