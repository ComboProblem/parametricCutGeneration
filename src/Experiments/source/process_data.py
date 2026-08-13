from cutgeneratingfunctionology.igp import *
from parametricCutGen.cut_generation_problem import *
from parametricCutGen.cgf_specializations import *
from parametricCutGen.experimental_utils.pyscipopt_data_collection_events import *
from pyscipopt import Model, quicksum, SCIP_PARAMSETTING, exp, log, sqrt, sin, SCIP_HEURTIMING
import json
import os
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

def scip():
    scip = Model()
    scip.hideOutput()
    return scip

pcg_exp_logger = logging.getLogger("parametricCutGen.experimental_utils.pyscipopt_data_collection_events")
pcg_exp_logger.setLevel(logging.ERROR)
functions_logger = logging.getLogger("cutgeneratingfunctionology.igp.functions")
functions_logger.setLevel(logging.ERROR)
r"""
Transform raw data (generated cgf log from an experimental .txt output) to a useable format as a .cgfl (cut generating function libary)

The .cgfl contains the following information

%model_name %metadata

%model_name is the name of the model
%metadata is resvered for future use.

    To include eventually::
    %version = { package : version for package in parametricatCutGeneration }
    %source = url_to_model.mps
    %solution = url_to_solution.sol
    %max_constraint_violation = real >= 0

For
 
%cut_gen_info %generation_params %stats

row is row of the MIP the cut was generated with. 

%cut_gen_info = ((b,v), row) the point in RR^{2n} corrosponding to a cut generating function (or approximation there of) and the row it generates a cut for.

%generation_params = {"cut_score": list_of_str, "chart_epsilon" : list_of_floats, "problem_dim": list_of_int }

%stats = {"max_constraint_volation" : float, "is_gmic": bool, "realitive_change_in_gap": float, "tree_depth_approx": int}
"""

cut_score_names = ['parallelism', 'cut_off_distance', 'violation', 'realitive_violation']
trial_eps_denom = [2*i for i in range(1,17)]
trial_cuts = [1]

bits_for_id = 6


def find_possible_f_index(fun):
    r"""
    Assume a model function, fun, with a finite number of breakpoints is given.
    Finds the index i such that minimizes |fun(lambda_i) - 1| 

    INPUT:
    - igp piecewise linear function

    OUTPUT:
    - integer or None
    """
    values = [fun(x) for x in  fun.end_points()]
    try:
        f_index = values.index(1)
    except Exception as e:
        f_index = np.argmin([abs(v-1) for v in values])
    return f_index
        

def as_sage_rational(dct):
    """
    TESTS::
    >>> sage_rational_zero = as_sage_rational({'__sage.rings.rational.Rational__': True, 'numerator': 0, 'denominator': 1})
    >>> sage_rational_zero in QQ
    True
    """
    if "__sage.rings.rational.Rational__" in dct:
        return QQ(dct["numerator"]/dct["denominator"])
    return dct

def load_and_process_data(model_name, sol_path, numerical_epsilon=1e-9, scip_time=90):
    oracle_sol, no_cuts_dual_value, tree_depth_without_cuts = get_oracle_sol_and_no_cuts_dual_value(model_name, sol_path, scip_time=scip_time, count=None)
    logger.debug(f"dual_bound_no_cuts: {no_cuts_dual_value}, tree depth: {tree_depth_without_cuts}")
    exp_data = {"bkpt_val_cgf":[None, None], "generation_params":[{"cut_score": None, "chart_epsilon":None, "problem_dim": None }], "stats":[{"max_constraint_violation" : None, "is_gmic": None , "dual_bound": no_cuts_dual_value, "tree_depth_approx": tree_depth_without_cuts}]}
    good_fun_count = 0
    for exp_id in range(64):
        bit_string =  '0'*(bits_for_id-exp_id.bit_length()) + bin(exp_id)[2:]
        cut_score_index = int(bit_string[0:2], 2)
        esp_index =  int(bit_string[2:], 2)
        number_of_cuts_index = 0 #int(bit_string[5:], 2)
        logger.debug(f"loading exp_id: {exp_id}")
        try:
            logger.debug(f"attempting to load: data/{model_name}.bkpt_as_param.{cut_score_names[cut_score_index]}.{trial_eps_denom[esp_index]}.{trial_cuts[number_of_cuts_index]}.txt")
            with open(f"data/{model_name}.bkpt_as_param.{cut_score_names[cut_score_index]}.{trial_eps_denom[esp_index]}.{trial_cuts[number_of_cuts_index]}.txt", 'r') as data_file:
                data = json.load(data_file)
            for datum in data["cgf_log"]:
                b = [ as_sage_rational(dct) for dct in datum[0][0] ]
                v = [ as_sage_rational(dct) for dct in datum[0][1] ]
                row = datum[4]
                fun = piecewise_function_from_breakpoints_and_values(b, v)
                max_constraint_violation = minimality_constraint_violation(fun)
                logger.debug(f"parm data: {(b,v)}, {row}")
                try:
                    if max_constraint_violation <= numerical_epsilon:
                        if ((b,v), row) not in exp_data["bkpt_val_cgf"]: 
                            tree_depth_approx, dual_bound = realitive_gap_and_tree_approx(fun, row, model_name, sol_path, oracle_sol, no_cuts_dual_value, scip_time=scip_time, count=good_fun_count)
                            good_fun_count += 1
                            logger.debug(f"function is approx minimal")
                            exp_data["bkpt_val_cgf"].append(((b,v), row))
                            exp_data["generation_params"].append({"cut_score": [cut_score_names[cut_score_index]], "chart_epsilon": [QQ(1/trial_eps_denom[esp_index])], "problem_dim": [datum[2]] })
                            exp_data["stats"].append({"max_constraint_violation" : max_constraint_violation, "is_gmic": is_gmic(fun, numerical_epsilon), "dual_bound": dual_bound, "tree_depth_approx": tree_depth_approx})
                            logger.debug(f"compute result: {tree_depth_approx}, {dual_bound-no_cuts_dual_value}")
                        else:
                            logger.debug(f"functions stats already computed, adding generation data.")
                            ind = exp_data["bkpt_val_cgf"].index(((b,v), row))
                            exp_data["generation_params"][ind]["cut_score"].append(cut_score_names[cut_score_index])
                            exp_data["generation_params"][ind]["chart_epsilon"].append(QQ(1/trial_eps_denom[esp_index]))
                            exp_data["generation_params"][ind]["problem_dim"].append(datum[2])
                    else:
                        if ((b,v), row) not in exp_data["bkpt_val_cgf"]: 
                            dual_bound = "DidNotCompute"
                            tree_depth_approx = "DidNotCompute"
                            logger.debug(f"function is not approx minimal, no computations")
                            exp_data["bkpt_val_cgf"].append(((b,v), row))
                            exp_data["generation_params"].append({"cut_score": [cut_score_names[cut_score_index]], "chart_epsilon": [QQ(1/trial_eps_denom[esp_index])], "problem_dim": [datum[2]] })
                            exp_data["stats"].append({"max_constraint_violation" : max_constraint_violation, "is_gmic": is_gmic(fun, numerical_epsilon), "dual_bound": dual_bound, "tree_depth_approx": tree_depth_approx})
                        else:
                            ind = exp_data["bkpt_val_cgf"].index(((b,v), row))
                            exp_data["generation_params"][ind]["cut_score"].append(cut_score_names[cut_score_index])
                            exp_data["generation_params"][ind]["chart_epsilon"].append(QQ(1/trial_eps_denom[esp_index]))
                            exp_data["generation_params"][ind]["problem_dim"].append(datum[2])
                    
                except Exception as e:
                    logger.error("error in running calcuations")
                    logger.error(e)
        except Exception as e:
            logger.info(f"No record found for experiment id {exp_id}")
    return exp_data

def est_process_time(model_name, sol_path, numerical_epsilon=1e-9, scip_time=90, seed=None, sample_size_per_exp=3, max_est_time=600):
    import random
    import time
    if seed is None:
        seed = random.seed()
    else:
        seed = random.seed(seed)
    logger.debug(f"Python random seed: {seed}")
    number_of_samples_recored=0
    total_number_of_trials_to_process=0
    total_time = 0
    start_no_cuts = time.time()
    oracle_sol, no_cuts_dual_value, tree_depth_without_cuts = get_oracle_sol_and_no_cuts_dual_value(model_name, sol_path, scip_time=scip_time, count=None)
    logger.debug(f"dual_bound_no_cuts: {no_cuts_dual_value}, tree depth: {tree_depth_without_cuts}")
    exp_data = {"bkpt_val_cgf":[None, None], "generation_params":[{"cut_score": None, "chart_epsilon":None, "problem_dim": None }], "stats":[{"max_constraint_violation" : None, "is_gmic": None , "dual_bound": no_cuts_dual_value, "tree_depth_approx": tree_depth_without_cuts}]}
    end_no_cuts = time.time()
    total_time += end_no_cuts - start_no_cuts 
    number_of_samples_recored+=1
    good_fun_count = 0
    for exp_id in range(64):
        bit_string =  '0'*(bits_for_id-exp_id.bit_length()) + bin(exp_id)[2:]
        cut_score_index = int(bit_string[0:2], 2)
        esp_index =  int(bit_string[2:], 2)
        number_of_cuts_index = 0 #int(bit_string[5:], 2)
        logger.debug(f"loading exp_id: {exp_id}")
        try:
            logger.debug(f"attempting to load: data/{model_name}.bkpt_as_param.{cut_score_names[cut_score_index]}.{trial_eps_denom[esp_index]}.{trial_cuts[number_of_cuts_index]}.txt")
            with open(f"data/{model_name}.bkpt_as_param.{cut_score_names[cut_score_index]}.{trial_eps_denom[esp_index]}.{trial_cuts[number_of_cuts_index]}.txt", 'r') as data_file:
                data = json.load(data_file)
            num_fun = len(data["cgf_log"])
            total_number_of_trials_to_process += num_fun
            logger.info(f"exp_id {exp_id} recorded {num_fun} functions.")
            if total_time < max_est_time:            
                indices = random.sample(range(num_fun), sample_size_per_exp)
                for datum in [data["cgf_log"][i] for i in indices]:
                    sample_start_time = time.time()
                    b = [ as_sage_rational(dct) for dct in datum[0][0] ]
                    v = [ as_sage_rational(dct) for dct in datum[0][1] ]
                    row = datum[4]
                    fun = piecewise_function_from_breakpoints_and_values(b, v)
                    max_constraint_violation = minimality_constraint_violation(fun)
                    logger.debug(f"parm data: {(b,v)}, {row}")
                    try:
                        if max_constraint_violation <= numerical_epsilon:
                            tree_depth_approx, dual_bound = realitive_gap_and_tree_approx(fun, row, model_name, sol_path, oracle_sol, no_cuts_dual_value, scip_time=scip_time, count=good_fun_count)
                            good_fun_count += 1
                            logger.debug(f"function is approx minimal")
                            exp_data["bkpt_val_cgf"].append(((b,v), row))
                            exp_data["generation_params"].append({"cut_score": [cut_score_names[cut_score_index]], "chart_epsilon": [QQ(1/trial_eps_denom[esp_index])], "problem_dim": [datum[2]] })
                            exp_data["stats"].append({"max_constraint_violation" : max_constraint_violation, "is_gmic": is_gmic(fun, numerical_epsilon), "dual_bound": dual_bound, "tree_depth_approx": tree_depth_approx})
                        else:
                            dual_bound = "DidNotCompute"
                            tree_depth_approx = "DidNotCompute"
                            logger.debug(f"function is not approx minimal, no computations")
                            exp_data["bkpt_val_cgf"].append(((b,v), row))
                            exp_data["generation_params"].append({"cut_score": [cut_score_names[cut_score_index]], "chart_epsilon": [QQ(1/trial_eps_denom[esp_index])], "problem_dim": [datum[2]] })
                            exp_data["stats"].append({"max_constraint_violation" : max_constraint_violation, "is_gmic": is_gmic(fun, numerical_epsilon), "dual_bound": dual_bound, "tree_depth_approx": tree_depth_approx})
                        
                    except Exception as e:
                        logger.error("error in running calcuations")
                        logger.error(e)
                    sample_end_time = time.time()
                    total_time += sample_end_time - sample_start_time
                    number_of_samples_recored+=1
                else:
                    continue

        except Exception as e:
            logger.info(f"No record found for experiment id {exp_id}")
    mean_time = total_time/number_of_samples_recored
    logger.debug(f"mean_time: {mean_time}, total_number_of_trials_to_process:{total_number_of_trials_to_process}")
    return (mean_time)*total_number_of_trials_to_process

def minimality_constraint_violation(fun):
    """
    Assume there exists an f in (0,1) s.t. fun(f) == 1.
    
    max_constraint_violation >= 0
    
    delta fun(x,y) := fun(x) + fun(y) - fun(x+y) 
    
    Evaluates the truth of the statement:
    evaluates min(delta fun(x,y)) for all x,y in R and
    For all x+y equiv f pmod 1, evaluate max(abs(delta fun(x, y) - 1)). 
    Report the least of these two values.

    If the value returend is 0, then fun is minimal.
    Otherwise, the value is maximum miniality constraint violation.

    Parameters::
    `fun` - 1DPWL 
    `max_constraint_violation` - element of [0,infty). 
    TESTS::
    >>> minimality_constraint_violation(gmic(4/5))
    0
    >>> f = piecewise_function_from_breakpoints_and_values([0, 1/3, 2/3, 1], [0, 1, 1/4, 0])
    >>> minimality_constraint_violation(f)
    .5
    """
    f_index = find_possible_f_index(fun)
    f_val = fun.end_points()[f_index]
    sym = [abs(float(delta_pi_eval) -1) for delta_pi_eval in generate_symmetric_vertices_continuous_expr(fun, f_val)]
    type_1 = [float(lhs) for lhs in generate_type_1_vertices_continuous_expr(fun)]
    type_2 = [float(lhs) for lhs in generate_type_2_vertices_continuous_expr(fun)]
    m_1 = min(type_1+type_2) # m_1 = 0  if and only if fun is subadditive
    if m_1 < 0:
        m_1 = abs(m_1)
    m_2 = max(sym) # m_2 = 0 if and only if fun is symmetric
    return max(m_1,m_2)

def is_gmic(fun, fun_equality=1e-9):
    """
    Suppose a function f in PWL{<=n} is given. 
    assume fun(f) == 1 and f is a breakpoint.
    Evaluate (||fun - gmic(f)||_inf <= fun_equality).
    
    Parameters::
    `fun` - 1DPWL 
    `fun_equality` - element of [0,infty). 
    """
    try:
        f_index = find_possible_f_index(fun)
    except Exception as e:
        f_index = None
    if f_index is not None:
        f = fun.end_points()[f_index]
        if inf_norm_of_cont_pwl(fun, gmic(f)) <= fun_equality: # from topology of PWL
            return True
    return False

def realitive_gap_and_tree_approx(fun, row, model_name, sol_path, oracle_sol, no_cuts_dual_value, scip_time=60, count=None):
    model = scip()
    model.readProblem(filename=f"model_files/{model_name}.mps")
    model.setParam("limits/time", scip_time)
    model.setSeparating(SCIP_PARAMSETTING.OFF)
    sepa = backTest(count, row, fun, model_name)
    model.includeSepa(sepa, "optimal_cut", "optimal cut over space of paramaterized cut generating functions", priority=10000, freq=0)
    model.setParam("separating/maxcutsroot", 1)
    model.setParam("separating/maxroundsroot", 1)
    model.setHeuristics(SCIP_PARAMSETTING.OFF)
    heuristic = OracleHeurisitc(sol_path)
    model.setPresolve(SCIP_PARAMSETTING.OFF)
    model.includeHeur(heuristic, "OracleHeurisitc", "for observing changes in dual bound from cuts", "Y", timingmask=SCIP_HEURTIMING.DURINGLPLOOP)
    model.setParam("limits/nodes", 1)
    try:
        model.optimize()
    except Exception as e:
        raise e
    try:
        model.writeMIP(filename=f"TEMP/{model_name}/{count}.lp")
    except Exception as e:
        logger.debug(e)
    model_lp_with_cut = scip()
    model_lp_with_cut.readProblem(filename=f"TEMP/{model_name}/{count}.lp")
    model_lp_with_cut.relax()
    model_lp_with_cut.optimize()
    sol = model_lp_with_cut.getVarDict()    
    sol_processed = {key : sol["t_"+key] for key in oracle_sol.keys() }
    dual_value = model_lp_with_cut.getDualbound()
    tree_depth_result = tree_depth_approx(sol_processed, oracle_sol)
    return tree_depth_result, dual_value

def get_oracle_sol_and_no_cuts_dual_value(model_name, sol_path, scip_time=60, count=None):
    model = scip()
    model.readProblem(filename=f"model_files/{model_name}.mps")
    model.setParam("limits/time", scip_time)
    model.setSeparating(SCIP_PARAMSETTING.OFF)
    sepa = backTest(count, None, None, model_name)
    model.includeSepa(sepa, "optimal_cut", "optimal cut over space of paramaterized cut generating functions", priority=10000, freq=0)
    model.setParam("separating/maxcutsroot", 1)
    model.setParam("separating/maxroundsroot", 1)
    model.setHeuristics(SCIP_PARAMSETTING.OFF)
    heuristic = OracleHeurisitc(sol_path)
    model.setPresolve(SCIP_PARAMSETTING.OFF)
    model.includeHeur(heuristic, "OracleHeurisitc", "for observing changes in dual bound from cuts", "Y", timingmask=SCIP_HEURTIMING.DURINGLPLOOP)
    model.setParam("limits/nodes", 1)
    model.optimize()
    oracle_sol = model.getVarDict()
    model.writeMIP(filename=f"TEMP/{model_name}/{count}.lp")
    model_root_lp_with_bounds = scip()
    model_root_lp_with_bounds.readProblem(filename=f"TEMP/{model_name}/{count}.lp")
    model_root_lp_with_bounds.relax()
    model_root_lp_with_bounds.optimize()
    sol = model_root_lp_with_bounds.getVarDict()
    sol_processed = {key : sol["t_"+key] for key in oracle_sol.keys() }
    no_cuts_dual_value = model_root_lp_with_bounds.getDualbound()
    tree_depth_without_cuts = tree_depth_approx(oracle_sol, sol_processed)
    return oracle_sol, no_cuts_dual_value, tree_depth_without_cuts

def tree_depth_approx(oracle_sol, sol_processed):
    return sum(int(abs(oracle_sol[key]-sol_processed[key]))+1 for  key in oracle_sol.keys())

def sage_rational_to_json(obj):
    if isinstance(obj, sage.rings.rational.Rational):
        return {'__sage.rings.rational.Rational__': True, 'numerator': int(obj.numerator()), 'denominator': int(obj.denominator())}
    raise TypeError(f'Cannot serialize object of {type(obj)}')

def __main__():
    model_name = "sp97ar" #os.getenv("MODEL_NAME")
    est_time = 2 #int(os.getenv("RUN_TIME")) #in minutes
    estimate_run = "y" #os.getenv("EST_RUN")
    sol_path = f"solution_files/solutions/{model_name}/1/{model_name}.sol"
    if estimate_run == "y":
        time_est = est_process_time(model_name, sol_path, max_est_time=(est_time-1.5)*60)
        alloc_slurm_time = int(time_est/60)+ 60
        with open(f"TEMP/{model_name}_data_run.sh", "w") as full_run:
            full_run.write(f"#!/bin/bash\nCLUSTER_ACCOUNT=$1\nPARTITION=$1\nsbatch --account=$CLUSTER_ACCOUNT --partition=$PARTITION --mem=5G --time={alloc_slurm_time}:00 --output=\"TEMP/{model_name}_result_ID_%a.out\" source/process_data.sh \"{model_name}\"")
            
    else:
        data = load_and_process_data(model_name, sol_path, scip_time=scip_time)
        data["metadata"] = {"cgp_version":"0.0.1alpha", "numerical_epsilon":1e-9, "model_name":model_name}
        with open(f"result/{model_name}.cgfl", "w") as cgfl:
            json.dump(data, cgfl, default=sage_rational_to_json)

__main__()
