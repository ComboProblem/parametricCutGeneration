#!/bin/bash

PARTITION="partition_name"
CLUSTER_ACCOUNT="account_name"

get_models() {
if ! [ -f "downloads/benchmark.zip" ] 
then
    echo "Downloading benchmark instances"
    wget http://miplib.zib.de/downloads/benchmark.zip -P downloads/
    unzip -u downloads/benchmark.zip -d model_files
    for compressed_model in model_files/*
    do 
        gzip -d $compressed_model
    done
else
    echo "Benchmark instance archive already exists."
fi
}
get_solutions() {
if ! [ -f "downloads/solutions.zip" ] 
then
    echo "Downloading all solutions."
    wget https://miplib.zib.de/downloads/solutions.zip -P downloads/
    unzip -u solution_files/solutions.zip -d solution_files
    for model in solution_files/solutions/*
    do
    	model_name=($basename $model)
	    gzip -d solution_files/solutions/$model_name/1/$model_name.sol.gz
    done
else
    echo "Solution archive already exists."
fi
}


check_apptainer() {
module load apptainer
if [ ! -f "container/parametricCutGen.sif" ]; then
  echo "Building Apptainer"
  apptainer build container/parametricCutGen.sif source/Apptainer.def
else
   echo "Optimal cut image already exists."
fi
}

run_experiments(){
for model in model_files/*
do
    model_name=$(basename $model .mps)
    echo "Queueing experiments for $model_name."
    sbatch --array=0--256 --account=$CLUSTER_ACCOUNT --partition=$PARTITION --time=20	0:00 --output="TEMP/$model_name_expID_%a.out"source/run_experiments.sh $model_name 
done
}

main() {
chmod +x source/run_experiments.sh
get_models
get_solutions
check_apptainer
mkdir TEMP
run_experiments
}

main
