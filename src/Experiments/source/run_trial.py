import logging
import os
import json
import tomllib
from pyscipopt import Model
from parametricCutGen.optimal_cut_generation import OptimalCut
from parametricCutGen.utils import validate_paths, parse_logger
from parametricCutGen.scip_data_collection_events import CutGapDataRecording

trial_logger = parse_logger(__name__)

shell_paths = { "experiment_trial_programs_path":os.getenv("EXPERIMENT_TRIAL_PROGRAMS_PATH"), "model":os.getenv("MODEL"),"data_target_path":os.getenv("DATA_TARGET_PATH"), "container":os.getenv("OPTIMAL_CUT_CONTAINER"), "experiment_params":os.getenv("EXP_PARAM_PATH"), "model_files":os.getenv("MODEL_FILES"), "conduct_experiment_base":os.getenv("PARAMETRIC_EXPS_BASE")  }

paths = validate_paths(shell_paths, trial_logger)

def scip_parameter_parser_and_model_loader(paths):
    """configures the model"""
    model = Model()
    # model.readParams(os.path.join(paths["experiment_trial_programs_path"], "scip_experiment_parmaters.set"))
    model.readParams(os.path.join(paths["experiment_params"], "scip_experimental_settings.set"))
    model.readProblem(filename=os.path.join(paths["model_files"], paths["model"])])
    return model

def run_trial(paths):
    model = scip_parameter_parser_and_model_loader(model, paths)
    # non_linear_solver = SciPy_paramater_parser(paths)
    experiment_parameters = json.loads(os.path.join(paths["experiment_trial_programs_path"], "experiment_parameters.json")
    seapa = OptimalCut(algorithm=experiment_parameters["algorithm"], backend=experiment_parameters["backend"], cut_score=experiment_parameters["cut_score"],  epsilon=experiment_parameters["epsilon"], M=experiment_parameters["M"], max_cgp_solver_time=experiment_parameters["max_cgp_solver_time"], max_num_of_bkpts=experiment_parameters["max_num_of_bkpts"], multithread=experiment_parameters["multithread"],
       prove_seperator=experiment_parameters["prove_seperator"] rel_tol=experiment_parameters["rel_tol"], show_proof=experiment_parameters["show_proof"])
    model.includeSepa(sepa, "optimal_cut_exp", "exp_params:{experiment}", priority=1000, freq=1)
    data_record = CutGapDataRecording(model, "optimal_cut_exp", experiment_parameters["max_number_of_cuts"])
    model.includeEventhdlr(data_record, "record_gap_data", "Records dual gap data when optimal_cut_exp is called" )
    model.optimize()
    data_record.write_data(paths["data_target_path"])
