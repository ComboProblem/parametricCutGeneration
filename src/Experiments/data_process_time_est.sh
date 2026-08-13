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
    sbatch --account=$CLUSTER_ACCOUNT --partition=$PARTITION --mem=5G --time=$SLURM_TIME:00 --output="TEMP/${model_name}_result_ID_%a.out" source/process_data.sh $model_name $EST_RUN $SLUM_TIME
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

