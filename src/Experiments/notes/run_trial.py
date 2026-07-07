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
