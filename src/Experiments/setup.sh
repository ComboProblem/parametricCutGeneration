#!/bin/bash

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
    unzip -u downloads/solutions.zip -d solution_files
    for model in solution_files/solutions/*
    do
    	model_name=$(basename $model)
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


main(){
get_models
get_solutions
check_apptainer
}

main
