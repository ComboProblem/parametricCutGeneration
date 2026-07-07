#!/bin/bash -x

get_models() {
if ! [ -f "model_files/benchmark.zip" ] 
then
    echo "Downloading benchmark instances"
    wget http://miplib.zib.de/downloads/benchmark.zip -P instances/
else
    echo "Benchmark instance archive already exists"
fi

unzip -u model_files/benchmark.zip -d model_files

fi
}

check_apptainer() {
module load apptainer
if [ ! -f $OPTIMAL_CUT_CONTAINER ]; then
  echo "Building Apptainer"
  apptainer build $OPTIMAL_CUT_CONTAINER source/Apptainer.def
else
   echo "Optimal cut image already exists."
fi
}

setup_experiments() {
srun --partition=$PARTITION --account=$CLUSTER_ACCOUNT --ntasks=1 --cpus-per-task=1 --time="$SETUP_TIME:00" --mem=$SETUP_MEM --wait ./src/Experiments/source/setup_experiments.sh
}

run_experiments(){
# The number of * is based on the number of non global parameters for experiments
# Currently the directory paths represent algorithm/cut_score/max_number_of_bkpts/max_number_of_cuts/
# see fun:setup_experiment_paths in setup_experiments.py 
for experiment in $PARAMETRIC_EXPS_BASE/*/*/*/*/TrialPrograms/run_trails.sh;
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
