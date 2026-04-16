#!/bin/bash -x

file_system_setup() {
mkdir -p $MODEL_FILES
mkdir -p $PARAMETRIC_EXPS
}

get_models() {
if ! [ -f "$MODEL_FILES/benchmark.zip" ]
then
    echo "Downloading benchmark instances"
    wget http://miplib.zib.de/downloads/benchmark.zip -P instances/
else
    echo "Benchmark instance archive already exists"
fi

unzip -u $MODEL_FILES/benchmark.zip -d $MODEL_FILES

fi
}

check_apptainer() {
module load apptainer
if [ ! -f "$./src/Experiments/source/optimal_cut.sif" ]; then
  echo "Building Apptainer"
  apptainer build "./src/Experiments/source/optimal_cut.sif" "./src/Experiments/source/Apptainer.def"
else
   echo "Optimal cut image already exists."
fi
}

setup_experiments() {
srun --partition=$PARTITION --account=$CLUSTER_ACCOUNT --ntasks=1 --cpus-per-task=1 --time="$INITAL_TIME:00" --mem=$BKPT_MEM --wait ./src/Experiments/source/setup_experiments.sh
}

run_experiments(){

}

main(){
echo "Loading cluster environment." 
source "./src/Experiments/source/cluster_environment.sh"
file_system_setup()
get_models()
check_apptainer()
}

main()
