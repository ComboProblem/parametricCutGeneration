#!/bin/bash -x

file_system_setup() {
mkdir -p $MODEL_FILES
mkdir -p $PARAMETRIC_EXPS_BASE
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
if [ ! -f $OPTIMAL_CUT_CONTAINER ]; then
  echo "Building Apptainer"
  apptainer build $OPTIMAL_CUT_CONTAINER $APPTAINER_DEF_PATH
else
   echo "Optimal cut image already exists."
fi
}

setup_experiments() {
srun --partition=$PARTITION --account=$CLUSTER_ACCOUNT --ntasks=1 --cpus-per-task=1 --time="$SETUP_TIME:00" --mem=$SETUP_MEM --wait ./src/Experiments/source/setup_experiments.sh
}

run_experiments(){
# for loop over somethign like a star command see example  Text/*/simple.sh see about globbing

}

main(){
echo "Loading cluster environment." 
source "./src/Experiments/source/cluster_environment.sh"
file_system_setup()
get_models()
check_apptainer()
}

main()
