import logging
import os
from pyscipopt import Model
from parametricCutGen.optimal_cut_generation import OptimalCut, scipyCutGenProbelmSolverInterface
from parametricCutGen.utils import validate_paths, parse_logger
from parametricCutGen.scip_data_collection_events import CutGapDataRecording

trial_logger = parse_logger(__name__)

paths = validate_paths({}}, set_experiments_logger)

def experiment_parameters_parser(paths)
    return experiment_parameters

def scip_parameter_parser_and_model_loader(paths):
    """configures the model"""
    
    return model
    
def SciPy_paramater_parser(paths):
    return non_linear_solver

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

