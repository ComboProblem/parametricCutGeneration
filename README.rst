Parametric Cut Generation
~~~~~~~~~~~~~~~~~~~~~~~~~

Parametric Cut Generation contains the source code for the python package ``parametricCutGen`` based on the dissertation of Acadia Larsen.

``parametricCutGen`` implements a single row optimal cut selection for Mixed Integer Programs over the domain (and restricted domains) of continuous minimal functions with at most k breakpoints. This repository is in an alpha state. 

Current Intents
~~~~~~~~~~~~~~~
* Illustrate concept of explicit optimal cut selection as a proof of concept for MIP solvers.
* Reproducibility of experimental data using a HPC.
* Demonstrate use of ``passsagemath`` and ``cutgeneratingfunctionology`` in application; in particular illustrate application of cutting edge mathematics to application of MIPs.
* Documentation and testing is minimal is intended to serve as a supplement to my dissertation. 

Installation
~~~~~~~~~~~~

This repository is only available from source with an installation of ``cutgeneratingfunctionlogy`` and ``passagemath`` required.::

    git clone --branch MinFunStable https://github.com/ComboProblem/cutgeneratingfunctionology.git
    cd cutgeneratingfunctionology
    python3 -m venv /cgf-venv/venv
    source /cgf-venv/venv/bin/activate
    pip install ".[passagemath]"
    pip install pplitepy
    git clone https://github.com/ComboProblem/MinimalFunctionCache.git
    cd MinimalFunctionCache
    pip install .
    cd ..
    git clone https://github.com/ComboProblem/parametricCutGeneration.git
    cd parametricCutGeneration
    pip install .

Examples
~~~~~~~~

From the python interpreter or sage interpreter the following code gives an example of solving cut generation problem::

    from parametricCutGen.cut_generation_problem import *
    from cutgeneratingfunctionology.igp import *
    import logging
    logging.disable()
    cgp_full_2 = cutGenerationProblem(algorithm="full", cut_score="steepest_direction", max_num_of_bkpts=2) # equiv to gmic
    cgp_value_polyhedron_2bkpt = cutGenerationProblem(algorithm="value_poly_lp", cut_score="steepest_direction", max_num_of_bkpts=2) # equiv to gmic
    cgp_bkpt_as_param_steepest_direction = cutGenerationProblem(algorithm="bkpt_as_param", cut_score="steepest_direction")
    cgp_bkpt_as_param_parallelism = cutGenerationProblem(algorithm="bkpt_as_param", cut_score="parallelism")
    cgp_value_polyhedron = cutGenerationProblem(algorithm="value_poly_lp", backend="pplite", cut_score="steepest_direction")
    cpg_gmic = [cgp_full_2, cgp_value_polyhedron_2bkpt]
    gpc_non_triv = [cgp_bkpt_as_param_steepest_direction , cgp_bkpt_as_param_parallelism, cgp_value_polyhedron]
    binvarow = [3.27, 4.66, 5.53, .56]
    binvc = [-1.2, -4.4, -5.6, -.1]
    f = 1.8 # aka b of the row

    print("gmic is recovered when expected.")
    print([inf_norm_of_cont_pwl(gmic(f=4/5), g) for g in [cpg.solve(binvarow, binvc, f) for cpg in cpg_gmic]])

    print("None trival solutions close to true solution.") # LP solvers and NL solvers used are not exact solvers, this is expected.
    sol_to_non_triv_problems = piecewise_function_from_breakpoints_and_values([0, .27, .53, .8, 1], [0, .9, .1, 1, 0])
    print([float(inf_norm_of_cont_pwl(sol_to_non_triv_problems, g)) for g in [cpg.solve(binvarow, binvc, f) for cpg in gpc_non_triv]])

Cut generation problems can be used in ``pyscipopt`` via optimal cut generation::

   from parametricCutGen.optimal_cut_generation import OptimalCut
   from piscipopt import Model
   model = Model()
   sepa = OptimalCut(cgp_kwds={"algorithm":"bkpt_as_param"})
   model.includeSepa(sepa, "optima_cut", "breakpoints a parameters problem", priority=1000, freq=1)


License 
~~~~~~~
The code is released under the GNU General Public License, version 2, or any later version as published by the Free Software Foundation.
