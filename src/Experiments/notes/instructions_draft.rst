=======
Locally
=======

The instructions assume you are using bash.

To start, clone parametricCutGeneration from source.  

```git clone https://github.com/ComboProblem/parametricCutGeneration.git```

Build ``parametricCutGeneration``.

``chmod +x parametricCutGeneration/local_build.sh``
``./parametricCutGeneration/local_build.sh -p``

Add ``cgp-env`` to the ``PATH``.

``PATH=$PATH:~/parametricCutGeneration/cgp-env``

If you wish to make this shortcut permate, add the last line to your ``~/.bashrc``.

Now start the ``cgp-env`` experimntal enviroment to begin experimenting.
Output is directed to the data folder. 


``cgp-env``
``(cgp-env)``
``(cgp-env) run_trial 64 gen-ip002 --scip_log ``
``(cgp-env) ls data``
``data/gen-ip002.no_cuts.cgfl data/gen-ip002.no_cuts.stats.json``
``(cgp-env) run_trial 0 gen-ip002``
``(cgp-env) ls data``
``gen-ip002.no_cuts.txt gen-ip002.no_cuts.stats.json gen-ip002.no_cuts.cgfl gen-ip002.bkpt_as_param.parallelism.2.1.1.cgfl``
``(cgp-env) trial_cgp --algorithm=bkpt_as_param --cut_score=cut_off_distance --eps_denom=4 --number_of_cuts=1 --model_name=gen-ip002 --scip_time=60``
``(cgp-env) ls data``
``gen-ip002.no_cuts.txt gen-ip002.no_cuts.stats.json gen-ip002.no_cuts.cgfl gen-ip002.bkpt_as_param.parallelism.2.1.1.cgfl gen-ip002.bkpt_as_param.parallelism.4.1.1.cgfl``
The last one runs a test of a ``cgp`` for given parameters, model, and scip time and ouputs the data into the data directory. 

Passing the ``--container`` or ``--c`` flag to ``cgp--env`` runs all commands using ``container/cgp.sif``.
This container will build on first run.
