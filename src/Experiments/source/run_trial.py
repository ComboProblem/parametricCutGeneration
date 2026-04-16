import logging
import os
import json
import tomllib
from pyscipopt import Model
from parametricCutGen.optimal_cut_generation import OptimalCut, scipyCutGenProbelmSolverInterface
from parametricCutGen.utils import validate_paths, parse_logger
from parametricCutGen.scip_data_collection_events import CutGapDataRecording

trial_logger = parse_logger(__name__)

shell_paths = { "experiment_trial_programs_path":os.getenv("EXPERIMENT_TRIAL_PROGRAMS_PATH"), "model":os.getenv("MODEL"),"data_target_path":os.getenv("DATA_TARGET_PATH"), "container":os.getenv("OPTIMAL_CUT_CONTAINER"), "experiment_params":os.getenv("EXP_PARAM_PATH"), "model_files":os.getenv("MODEL_FILES"), "conduct_experiment_base":os.getenv("PARAMETRIC_EXPS_BASE")  }

paths = validate_paths(shell_paths, trial_logger)

def experiment_parameters_parser(paths)
    return json.loads(os.path.join(paths["experiment_trial_programs_path"], "experiment_parameters.json")

def scip_parameter_parser_and_model_loader(paths):
    """configures the model"""
    with open(os.path.join(paths["experiment_params"], "scip_experimental_settings.toml") as f:
        parsed_toml_file = tomllib.loads(f.read())
    
    return model

def run_trial(paths):
    model = scip_parameter_parser_and_model_loader(model, paths)
    non_linear_solver = SciPy_paramater_parser(paths)
    experiment_parameters = experiment_parameters_parser(paths)
    seapa = OptimalCut(algorithm_name=experiment_parameters["algorithm"], cut_score=experiment_parameters["cutScore"], num_bkpt=experiment_parameters["num_bkpt"], multithread=False, prove_seperator=experiment_parameters["prove_seperator"], show_proof = False, epsilon=experiment_parameters["model_epsilon"], M = experiment_parameters["liptitz_constant"], nonlinear_backend=non_linear_solver)
     
    # add seperator to model
    # load data
    # solve model
    # record data
    # write data
    # First line is experiment parameters; it should contain spepific pamaters, name of model file, ect. 
    # rest is data; i should collect the .out file

