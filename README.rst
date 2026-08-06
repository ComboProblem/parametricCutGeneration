Parametric Cut Generation
~~~~~~~~~~~~~~~~~~~~~~~~~

Parametric Cut Generation contains the source code for the python package ``parametricCutGen`` based on the dissertation of Acadia Larsen.

``parametricCutGen`` implements a single row optimal cut selection for Mixed Integer Programs over the domain (and restricted domains) of continuous minimal functions with at most k breakpoints. 

This repository is in an alpha state. Version 0.1 coming soon!

Current Intents
~~~~~~~~~~~~~~~

* Illustrate concept of parametric cut generation and optimal cut generation as a proof of concept for MIP solvers.
* Reproducibility of experimental data using a HPC.
* Demonstrate use of ``passsagemath`` and ``cutgeneratingfunctionology`` in application; in particular illustrate application of cutting edge mathematics to application of MIPs.
* Documentation is intended support to my dissertation.

Installation
~~~~~~~~~~~~

This repository is currerntly only available from source.


An installation of ``cutgeneratingfunctionlogy``, ``passagemath``, ``pplitepy``, ``pyscipopt``, ``scipy``,  and ``cvxpy``  are required.
To install::

    git clone https://github.com/ComboProblem/parametricCutGeneration.git
    cd parametricCutGeneration
    python3 -m venv /cgp-env/venv
    source /cgp-env/venv/bin/activate
    git clone --branch MinFunStable https://github.com/ComboProblem/cutgeneratingfunctionology.git
    cd cutgeneratingfunctionology
    pip install '.[passagemath]'
    pip install cvxpy
    pip install scipy
    pip install pyscipopt
    pip install pplitepy    
    pip install .

If you wish to run a container, an ``apptainer`` ``.def`` file is provided. The recommened build is given.::

    git clone https://github.com/ComboProblem/parametricCutGeneration.git
    cd parametricCutGeneration
    apptainer build src/Experiments/source/Apptainer.def src/Experiments/container/cgp.sif
    apptainer run src/Experiments/container/cgp.sif bash

See ``src/Experiments/readme.rst`` for details about use with a cluster.


Examples
~~~~~~~~

Cut generation problems can be used in ``pyscipopt`` via optimal cut generation.
In python, ``OptimalCut`` can be added in the following way.::

   from parametricCutGen.optimal_cut_generation import OptimalCut
   from pyscipopt import Model
   model = Model()
   sepa = OptimalCut(cgp_kwds={'algorithm':'bkpt_as_param', 'backend':'pplite', 'cut_score':'parallelism',  'epsilon': 1/4, 'M':1e6})
   model.includeSepa(sepa, 'optima_cut', 'Optimally generated cuts using breakpoints as parameters algorithm', priority=1000, freq=1)


See `src/Experiments/README.rst` for a quick primier on optimal cut parameters.

License 
~~~~~~~
The code is released under the GNU General Public License, version 2, or any later version as published by the Free Software Foundation.
