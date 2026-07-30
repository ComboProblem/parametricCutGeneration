from cutgeneratingfunctionology.igp import *
from parametricCutGen.cut_generation_problem import *
from parametricCutGen.cgf_specializations import *
from pyscipopt import Model, quicksum, SCIP_PARAMSETTING, exp, log, sqrt, sin, SCIP_HEURTIMING
import os


model_name = os.getenv("MODEL")

cut_score_names = ['parallelism', 'cut_off_distance', 'violation', 'realitive_violation']
trial_eps_denom = [2*i for i in range(1,17)]
trial_cuts = [1]

bits_for_id = 6

def as_sage_rational(dct):
    if "__sage.rings.rational.Rational__" in dct:
        return QQ(dct["numerator"]/dct["denominator"])
	return dct

def load_exp(exp_id):
    if exp_id < 64:
        bit_string =  '0'*(bits_for_id-exp_id.bit_length()) + bin(exp_id)[2:]
        cut_score_index = int(bit_string[0:2], 2)
        esp_index =  int(bit_string[2:], 2)
        number_of_cuts_index = 0 #int(bit_string[5:], 2)
	    dual_log_with_cuts = {"row" : [], "ncuts": [], "dual_value" : [] }
	    fun_data = {"row" : [], "ncuts" : [], "fun": [], "score":[], "problem_dim":[]}
	    preprocessed_data = {"exp_id": exp_id,  "dual_log": dual_log_with_cuts, "fun_data" : fun_data} 
        with open(f"data/{model_name}.bkpt_as_param.{cut_score_names[cut_score_index]}.{trial_eps_denom[esp_index]}.{trial_cuts[number_of_cuts_index]}.txt", 'r') as data_file:
            data = json.loads(model.data, data_file, object_hook=sage_rational_decoder)
        for log in data["dual_log"]:
            if len(log[2]) == 1:
                cut_name = log[3].spit(":")[0]
                ncut, row = cut_name.split("optimal_cut")[1].split("_x")
                dual_log_with_cuts["ncuts"].append(int(ncut))
                dual_log_with_cuts["row"].append(int(row))
                dual_log_with_cuts["dual_value"].append(log[1])
        for paramaterized_data in data["cgf_log"]:
            fun_data["row"].append(data["cgf_log"][4])
            fun_data["ncuts"].append(data["cgf_log"][3])
            b = data["cgf_log"][0][0]
            v = data["cgf_log"][0][1]
            fun_data["fun"].append(piecewise_function_from_breakpoints_and_values(b, v))
            fun_data["problem_dim"].append(data["cgf_log"][2])
            fun_data["score"].append(data["cgf_log"][1])
    elif exp_id == 64:
	    dual_log_with_cuts = {"row" : None, "ncuts": None, "dual_value" : [] }
	    fun_data = {"row" : None, "ncuts" : None, "fun": None, "score": None, "problem_dim":None}
	    preprocessed_data = {"exp_id": exp_id,  "dual_log": dual_log_with_cuts, "fun_data" : fun_data}
	return preprocessed_trial
	
def preprocess_data(model_name):
	# load all trials
    preprocessed_data = {"exp_id": []}
	for exp_id in range(65):
		try:
            preprocessed_trial = load_exp(exp_id)
			preprocessed_data["exp_id"].append(exp_id)
            preprocessed_data[exp_id] = preprocessed_trial
		except Exception as e: # if the trial has failed to run, we will get a loading error, just contiue here.
			continue
    # compare accross cut scores, fix eps and a row. 
    return preprocessed_data

def process_minimality_constraint_violation(preprocessed_data, max_minimality_constraint_violation=1e-9):
	"""Compute max constraint violation of minimality test. """
    def minimality_constraint_violation(fun):
	    for lhs in generate_type_1_vertices_continuous_expr(fun):
            if lhs > -1*max_minimality_constraint_violation:
                return False
	    for lhs in generate_type_2_vertices_continuous_expr(fun):
            if lhs > -1*max_minimality_constraint_violation:
                return False
	    f_index = find_f_index(fun)
	    f = fun.end_points()[f_index]    
        for lhs in generate_symmetric_vertices_continuous_expr(fun, f):
            if abs(lhs - 1) > max_minimality_constraint_violation:
                return False
        return True
    functions_are_approx_minimal = lambda exp_id : [minimality_constraint_violation(fun) for fun in preprocessed_data[exp_id]["fun_data"]["fun"]]
    # report success rate of this. (should be 100%)
    result = {exp_id : functions_are_approx_minimal(exp_id) for exp_id in preprocessed_data["exp_id"]}
    return result

def process_compare_functions(preprocess_data, is_all_gmic, fun_equality=1e-9):
    result = {} # row result will be 
    for eps_index in range(2, 17):
        exp_ids = [int(bin(i)[2:]+bin(eps_index)[2:], 2) for i in range(4)]
        if all((exp_id in preprocess_data) for exp_id in exp_ids):
            if not all( [all(is_all_gmic[exp_id]) for exp_id in preprocess_data["exp_id"]]):
                min_no_cgp_solved = min([len(preprocessed_data[exp_id]["fun_data"]) for exp_id in exp_ids])
                result[trial_eps_denom[eps_index]] = { "min_no_cgp_solved": min_no_cgp_solved, "row_result" : []}                 
                for i in range(min_no_cgp_solved):
                    # use zip or iterator combinatorics.
                    res = []
                    for pair in permuations(range(4), 2):
                        f = preprocessed_data[exp_ids[pair[0]]]["fun_data"]["fun"][i]
                        f_tilde = preprocessed_data[exp_ids[pair[1]]]["fun_data"]["fun"][i]
                        if function_norm(f, f_tilde) < fun_equality:
                            res.append(True)
                        else:
                            res.append(False)
                    result[trial_eps_denom[eps_index]]["row_result"].append(res)
        else:
            continue
    return result




def process_if_functions_are_gmic(preprocessed_data, fun_equality=1e-9):
    def is_gmic(fun):
	    f_index = find_f_index(fun)
	    f = fun.end_points()[f_index]
	    if inf_norm_of_cont_pwl(fun, gmic(f)) < fun_equality:
		    return True
	    return False
    functions_are_gmic = lambda exp_id : [is_gmic(fun) for fun in preprocessed_data[exp_id]["fun_data"]["fun"]]
	result = {exp_id : functions_are_gmic(exp_id) for exp_id in preprocessed_data["exp_id"]}
    return result	

def 


def approx_tree_data(preprocessed_data, model_name, sol_path):
    result = {}
    model = Model()
    model.readProblem(filename=f"model_files/{model_name}.mps)
    model.setParam("limits/time", scip_time)
    model.setSeparating(SCIP_PARAMSETTING.OFF)
    sepa = approxTreeDepth(None, None, model_name)
    model.includeSepa(sepa, "optimal_cut", "optimal cut over space of paramaterized cut generating functions", priority=10000, freq=0)
    model.setParam("separating/maxcutsroot", 1)
    model.setParam("separating/maxroundsroot", 1)
    model.setHeuristics(SCIP_PARAMSETTING.OFF)
    heuristic = OracleHeurisitc(sol_path)
    model.setPresolve(SCIP_PARAMSETTING.OFF)
    model.includeHeur(heuristic, "OracleHeurisitc", "for observing changes in dual bound from cuts", "Y", timingmask=SCIP_HEURTIMING.DURINGLPLOOP)
    model.setParam("limits/nodes", 1)
    model.optimize()
    model_lp_with_cut = Model()
    model_lp_with_cut.readProblem(filename=f"TEMP/{model_name}.{row}.lp")
    model_lp_with_cut.optimize()
    no_cuts_sol = model_lp_with_cut.getSols()
    no_cuts_dual_value = model_lp_with_cut.getDualbound()
    result[exp_id]
    for exp_id in preprocessed_data["exp_id"] if exp_id != 64:
        for row,fun in zip(preprocessed_data[exp_id]["fun_data"]["row"], preprocessed_data[exp_id]["fun_data"]["fun"])
            model = Model()
            model.readProblem(filename=f"model_files/{model_name}.mps)

            #base tree depth - fun_base_tree > 0 is good meaning smaller estimated tree size,  <= no poential savings,
            # factors into tree depth, problem size, cut gen problem size
            model.setParam("limits/time", scip_time)
            model.setSeparating(SCIP_PARAMSETTING.OFF)
            sepa = approxTreeDepth(row, fun, model_name)
            model.includeSepa(sepa, "optimal_cut", "optimal cut over space of paramaterized cut generating functions", priority=10000, freq=0)
            model.setParam("separating/maxcutsroot", 1)
            model.setParam("separating/maxroundsroot", 1)
            model.setHeuristics(SCIP_PARAMSETTING.OFF)
            heuristic = OracleHeurisitc(sol_path)
            model.setPresolve(SCIP_PARAMSETTING.OFF)
            model.includeHeur(heuristic, "OracleHeurisitc", "for observing changes in dual bound from cuts", "Y", timingmask=SCIP_HEURTIMING.DURINGLPLOOP)
            model.setParam("limits/nodes", 1)
            model.optimize()
            model_lp_with_cut = Model()
            model_lp_with_cut.readProblem(filename=f"TEMP/{model_name}.{exp_id}.{row}.lp")
            model_lp_with_cut.optimize()
            sol = model_lp_with_cut.getSols()
            dual_value = model_lp_with_cut.getDualbound()
            result[exp_id]
    return result
    
