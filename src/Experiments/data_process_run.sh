#!/bin/bash

PARTITION="partition_name"
CLUSTER_ACCOUNT="account_name"


process_data(){
JOB_COUNT=0
for model in model_files/*
do
    model_name=$(basename $model .mps)
    echo "Queing data processing pipeline for $model_name."
    data_run_file_name="${model_name}_data_run.sh"
    if [ -f "TEMP/${data_run_file}" ]; then
        echo "Time estimated file exists."
        chmod +x ./TEMP/$data_run_file_name
         ./TEMP/$data_run_file_name $CLUSTER_ACCOUNT $PARTITON
    else
        sbatch --account=$CLUSTER_ACCOUNT --partition=$PARTITION --mem=5G --time=120:00 --output="TEMP/${model_name}_result_ID_%a.out" source/process_data.sh $model_name i
    fi
    
    ((JOB_COUNT=JOB_COUNT+1))
    if (( JOB_COUNT == 80 )); then
        echo "Giving some time for other jobs to finish."
        sleep 2h
        JOB_COUNT=0
    fi
done
}


main() {
chmod +x source/estimate_process_time.sh
chmod +x source/process_data.sh
#mkdir TEMP
process_data

}

main

