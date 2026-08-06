===========
Experiments
===========

Cut generation is performed based on parameters specified by the experiment idetification nonnegative integer. 
Generally this will be some version of ``exp_ID`` 
``exp_ID`` is converted to bytes then mapped to an experiment. 
The greatest two bits corrospond to different cuts scores: ``'parallelism'``, ``'cut_off_distance'``, ``'violation'``, ``'realitive_violation'``.
The remaining bits corrospond to a denominator of the chart espislon, that is : epsilon = 1/(2*[1,17)), sequentially. 
64 is reserved for a baseline with no optimally generated cuts. 

Model are pulled in from the MIP 2017 Benchmark Libary. (https://miplib.zib.de/)

Model files are assumed to have the ``path model_files/$model_name.mps``.
Model are generally referenced by name either "model_name" or "model" when name is implied.

Solution files are used for testing are assumed to have the path ``solution_files/solutions/$model_name/1/$model_name.sol``

cgp problem output is recored as ``data/model_name.bkpt_as_param.cut_score.eps_denom.number_of_cuts.txt``

scip output is recored as ``data/model_name.bkpt_as_param.cut_score.eps_denom.number_of_cuts.stats.json``

==============
On the Cluster
==============

To start, clone parametricCutGeneration from source.  

```git clone https://github.com/ComboProblem/parametricCutGeneration.git```


The cluster must be using ``SLURM`` cluster (https://slurm.schedmd.com/overview.html).
and supports ``apptainer`` (https://apptainer.org/).

Use your favorite editor such as ``vim`` to add ``src/Experiments/run_single_trial.sh`` and ``src/Experiments/run_experiments.sh`` to the ``.gitignore``.


```vim .gitignore```

Go to the experiments subdirectory.

```cd parametricCutGen/src/Experiements```

Use your favorite editor such as ``vim`` to edit your SLURM account infromatiom  in ``run_single_trial.sh`` and ``run_experiments.sh``. 

```vim run_single_trial.sh```.
```vim run_experiments.sh```.

Runing a single trial of an experiment is simple. For example

```chmod +x ./run_single_trial.sh
   ./run_single_trial.sh 64 gen-ip002
``` 

Runs trial 64 on the model gen-ip002 from the 2017 MIP LIB benchmark set. 
This trial corrosponds to the no cut generation baseline.
The current configuration assumes 200 minutes of run time per trial using a single core and 10G of memory.
Modifying this could have adverse effects.
In future versions, alternative run time configurations will be supported. 

In general ```./run_single_trial.sh exp_ID model_name``. 

If you wish to reproduce all experiments (about 3 weeks in real time with at ~ 50 jobs running concurrently),

```chmod +x ./run_experiments.sh
   ./run_experiments.sh
``` 

========
WARNINGS
========

Current imementation runs slowly.

Each trial can take greater than the specified 3 hours SCIP time limit due to ``scip`` signals not being passed to ``parametricCutGen``.
When running locally, jobs may not stop in time.
``SLURM`` will kill jobs that are over time.
Future implemenations may fix this.
Be warned that ``trial command may exceed ``scip`` set time (mesured in seconds). 

``data`` is a tracked folder. 

====================
Optimal Cut Data Set
====================

``data`` is a tracked folder. 


