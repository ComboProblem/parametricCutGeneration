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

get_solutions() {
if ! [ -f "solutions_files/solutions.zip" ] 
then
    echo "Downloading benchmark instances"
    wget http://miplib.zib.de/downloads/solutions.zip -P instances/
else
    echo "Solution instance archive already exists"
fi

unzip -u solutions_files/solutions.zip -d solutions_files

fi
}

check_apptainer() {
module load apptainer
if [ ! -f container/optimal_cut.sif ]; then
  echo "Building Apptainer"
  apptainer build container/optimal_cut.sif source/Apptainer.def
else
   echo "Optimal cut image already exists."
fi
}

run_experiments(){
for problem in problem_names:
    sbatch --array=-1--287 --partition=$PARTITION --account=$CLUSTER_ACCOUNT --ntasks=1 --time="390:00" --mem=10GB run_trial.sh $problem
}

main(){
echo "Loading cluster environment." 
source "./src/Experiments/source/cluster_environment.sh"
get_models()
check_apptainer()
run_experiments()
}

main()
