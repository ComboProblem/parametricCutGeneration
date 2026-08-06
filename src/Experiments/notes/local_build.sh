#!/usr/bin/bash

python_setup(){
python3 -m venv /cgp-env/venv
source /cgp-env/venv/bin/activate
git clone --branch MinFunStable https://github.com/ComboProblem/cutgeneratingfunctionology.git
cd cutgeneratingfunctionology
pip install ".[passagemath]"
pip install cvxpy
pip install scipy
pip install pyscipopt
pip install pplitepy
}

container_setup(){
apptainer build container/cgp.sif src/Experiments/source/Apptainer.def
}

setup(){
while getopts `cp` flag; do
case "${flag}" in 
   c) container_setup()
   p) python_setup()

}

