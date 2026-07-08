#!/bin/bash -x

export PARTITION=""
export CLUSTER_ACCOUNT=""

get_models() {
if ! [ -f "model_files/benchmark.zip" ] 
then
    echo "Downloading benchmark instances"
    wget http://miplib.zib.de/downloads/benchmark.zip -P instances/
else
    echo "Benchmark instance archive already exists."
fi

unzip -u model_files/benchmark.zip -d model_files

fi
}
get_solutions() {
if ! [ -f "solution_files/benchmark.zip" ] 
then
    echo "Downloading all solutions."
    wget http://miplib.zib.de/downloads/benchmark.zip -P instances/
else
    echo "Solution archive already exists."
fi

unzip -u model_files/benchmark.zip -d model_files

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
# The number of * is based on the number of non global parameters for experiments
# Currently the directory paths represent algorithm/cut_score/max_number_of_bkpts/max_number_of_cuts/
# see fun:setup_experiment_paths in setup_experiments.py 
for experiment in filename="${fullfile##*/}"
do 
    chmod +x $experiment
    echo "Dispatching trial runs for $experiment."
    ./$experiment
done
}

main(){
echo "Loading cluster environment." 
source "./src/Experiments/source/cluster_environment.sh"
file_system_setup()
get_models()
check_apptainer()
run_experiments()
}

main()
