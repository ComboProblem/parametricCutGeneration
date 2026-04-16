import tomllib
import logging
import os
from pyscipopt import Model
from parametricCutGen.optimal_cut_generation import OptimalCut, scipyCutGenProbelmSolverInterface

experiments_logger = logging.getLogger(__name__)

def validate_paths():
    paths = {}
    paths_are_valid = True
    bad_paths = {}
    for path in paths:
        if os.path.exists(paths[path]):
            rep_elm_worker_logger.debug(f"Path {path} with dir {paths[path]} exists.")
        else:
            rep_elm_worker_logger.debug(f"Path {path} with dir {paths[path]} does not exists.")
            paths_are_valid = False
            bad_paths[path] = paths[path]
    if paths_are_valid is False:
        experiments_logger.error(f"All paths should exists. We have the following bad paths: {bad_paths}")
    return paths

def scip_parameter_parser(model, path_to_scip_config_toml):
    """configures the model"""
    lo
    pass
    
def SciPy_paramater_parser(path_to_SciPy_config_toml):
    pass

def run_trial(path_to_scip_config_toml, path_to_SciPy_config_toml, path_to_problem_file, **experiment_parameters):
    experimental_model = Model()
    scip_parameter_parser(experimental_model, scip_parameter_parser)
    non_linear_solver = SciPy_paramater_parser(path_to_SciPy_config_toml)
    seapa = OptimalCut(algorithm_name=experiment_parameters["algorithm"], cut_score=experiment_parameters["cutScore"], num_bkpt=experiment_parameters["num_bkpt"], multithread=False, prove_seperator=experiment_parameters["prove_seperator"], show_proof = False, epsilon=experiment_parameters["model_epsilon"], M = experiment_parameters["liptitz_constant"], nonlinear_backend=non_linear_solver)
    # add seperator to model
    # load data
    # solve model
    # record data
    # write data
    # First line is experiment parameters; it should contain spepific pamaters, name of model file, ect. 
    # rest is data; i should collect the .out file

def unpack_parameters(paths):
    
    return trial_list

# this gets ran as the first request. 
def write_trials_for_experiments(paths):
    """
    Writes all experimental trials to trial directory and defines a scrip to run all trials.
    """
    trial_number = 0
    # while loop for writing trial_$NUMBER.sh = f"trail_{trial_number}.sh"
    trial_list = unpack_parameters(paths) # path_to_experiment_tom_file, this should unpack experimental parameters and define all trials
    while len(trial_list) > 0:
        trial_info = trial_list.pop()
        trial_file_name = f"trial_{trial_number}.sh"
        with open(os.path.join(paths[], trial_file_name), w) as trial_file:
            trial_file.write("#!/bin/bash\n")
            trial_file.write("source cluster_enviroment.sh\n")
            trial_file.write("echo RUNNING TRIAL NUMBER: $SLURM_ARRAY_TASK_ID \n")
            trial_file.write("module load apptainer\n")
            apptainer_confg_line = f"apptainer exec optimal_cut.sif run_trial({path_to_scip_config_toml}, {path_to_SciPy_config_toml},... )"
            trial_file.write(apptainer_confg_line)
        trial_number += 1
    with open(os.path.join(paths[], "run_trials.sh"), w) as run_trials:
        run_trials.write("#!/bin/bash\n")
        run_trials.write("source cluster_enviroment.sh\n")
        run_trials.write("echo \"Running on $(hostname)\"\n echo \"TASK: $SLURM_ARRAY_TASK_ID\" \n")
        run_trials.write("chmod +x \"$TRIAL_FILES_TO_RUN\\trial_$SLURM_ARRAY_TASK_ID.sh\"\n")
        run_trials.write(".\\$TRIAL_FILES_TO_RUN\\trial_$SLURM_ARRAY_TASK_ID.sh\"\n")
    with open(os.path.join(paths[], "trials_env.sh"), w) as trial_env:
        trial_env.write("#!/bin/bash\n")
        trial_env.write(f"export NUM_TRIAL={data_queue_length}")
    
