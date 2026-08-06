import logging
import os
import json
import tomllib
from pyscipopt import Model
from parametricCutGen.optimal_cut_generation import OptimalCut
from parametricCutGen.utils import validate_paths, parse_logger
from parametricCutGen.scip_data_collection_events import record_data

trial_logger = parse_logger(__name__)

shell_paths = { "experiment_trial_programs_path":os.getenv("EXPERIMENT_TRIAL_PROGRAMS_PATH"), "model":os.getenv("MODEL"),"data_target_path":os.getenv("DATA_TARGET_PATH"), "container":os.getenv("OPTIMAL_CUT_CONTAINER"), "experiment_params":os.getenv("EXP_PARAM_PATH"), "model_file":os.getenv("MODEL_FILE"), "conduct_experiment_base":os.getenv("PARAMETRIC_EXPS_BASE")  }

paths = validate_paths(shell_paths, trial_logger)

def run_trial_node_evolution(paths):
    model = Model()
    cgp_experiment_kwrds = json.loads(os.path.join(paths["experiment_trial_programs_path"], "experiment_parameters.json")
    write_path = paths["metadata_write_path"]
    seapa = OptimalCut(cgp_kwds=cgp_experiment_kwrds["cpg_kwds"])
    model.setSeparating(SCIP_PARAMSETTING.OFF)
    model.setHeuristics(SCIP_PARAMSETTING.OFF)
    model.setPresolve(SCIP_PARAMSETTING.OFF)
    model.includeSepa(sepa, "optimal_cut", "experiment cuts", priority=10000, freq=0)
    # Add exactly k cuts at the root.
    model.setParam("separating/maxcutsroot", 1)
    model.setParam("separating/maxroundsroot", cgp_experiment_kwrds["max_number_of_cuts"])
    model.setParam("limits/nodes", 10000000)
    model = record_data(model, write_path)
    record_data = record_data(model) 
    model.readProblem(paths["model_file"])
    try:
        model.optimize()
    except Exception as e:
        
    model.data

    
#elif exp_id > 255:
#    exp_id = exp_id - 255
#    bit_string =  '0'*(6-exp_id.bit_length()) + bin(exp_id)[2:]
#    number_breakpoints_index =  int(bit_string[0:3], 2)
#    number_of_cuts_index =  int(bit_string[3:], 2)
#    cpg_kwds = {'algorithm':'value_poly_lp', "max_bkpt":trial_eps_denom[esp_index]}
#    numb_cuts =  trial_cuts[number_of_cuts_index]
#    numb_cuts = trial_cuts[number_of_cuts_index]
#    model = Model()
#    model.readProblem(filename=model_path)
#    model.setParam("limits/time", scip_time)
#    model.setSeparating(SCIP_PARAMSETTING.OFF)
#    sepa = OptimalCut(write_cgf_data=True, cgp_kwds=cpg_kwds)
#    model.includeSepa(sepa, "optimal_cut", "optimal cut over space of paramaterized cut generating functions", priority=10000, freq=0)
#    model.setParam("separating/maxcutsroot", 1)
#    model.setParam("separating/maxroundsroot", numb_cuts)
#    model.setHeuristics(SCIP_PARAMSETTING.OFF)
#    heuristic = OracleHeurisitc(sol_path)
#    model.includeHeur(heuristic, "OracleHeurisitc", "for observing changes in dual bound from cuts", "Y", timingmask=SCIP_HEURTIMING.DURINGLPLOOP)
#    model.setPresolve(SCIP_PARAMSETTING.OFF)
#    model.setParam("limits/nodes", 1)
#    model=record_data(model)
#    model.hideOutput()
#    model.optimize()
#    model.writeStatistics(filename=f"data/{model_name}.bkpt_as_param.value_poly_lp.bkpt_as_param.{cut_score_names[cut_score_index]}.{trial_eps_denom[esp_index]}.{numb_cuts}.out")
#    with open(f"data/{model_name}.value_poly_lp.bkpt_as_param.{cut_score_names[cut_score_index]}.{trial_eps_denom[esp_index]}.{numb_cuts}.txt", 'w') as data_file:
#        json.dump(model.data, data_file, default=sage_rational_to_json)
#
#except Exception as e:
#    logging.debug(e)
#    print()
#    try:
#        data = {"Exception": str(e)} | model.data
#    except Exception as e2:
#        data = {"Exception": str(e), "failure": str(e2)}
