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
    cgp_experiment_kwrds = json.loads(os.path.join(paths["experiment_trial_programs_path"], "experiment_parameters.json")
    paths["metadata_write_path"]=path
    seapa = OptimalCut(write_mip_and_cut=True, cgp_kwds=cgp_experiment_kwrds, paths=paths)
    model.includeSepa(sepa, "optimal_cut_exp", "exp_params:{experiment}", priority=1000, freq=1)
    data_record = CutGapDataRecording(model, "optimal_cut_exp", experiment_parameters["max_number_of_cuts"])
    model.includeEventhdlr(data_record, "record_gap_data", "Records dual gap data when optimal_cut_exp is called" )
    model.optimize()
    data_record.write_data(paths["data_target_path"])
