import tomllib
import logging
import os
from pyscipopt import Model
from parametricCutGen.optimal_cut_generation import OptimalCut, scipyCutGenProbelmSolverInterface

experiments_logger = logging.getLogger(__name__)

def scip_parameter_parser(model, path_to_scip_config_toml):
    """configures the model"""
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

# this gets ran as the first request. 
def write_trials_for_experiments(path_to_experiment_tolm_file, path_to_model_files, experimal_run_TEMP):
    trial_number = 0
    # while loop for writing trial_$NUMBER.sh = f"trail_{trial_number}.sh"
    data_queue = None # path_to_experiment_tom_file, this should unpack experimental parameters and define all trials
    data_queue_length = data_queue.length()
    os.chdir(os.getenv("TRIAL_FILES_TO_RUN"))
    while not data_queue.is_empty():
        trial_info = data_queue.get_next_trial()
        trial_file_name = f"trial_{trial_number}.sh"
        with open(trial_file_name, w) as trial_file:
            trial_file.write("#!/bin/bash\n")
            trial_file.write("source cluster_enviroment.sh\n")
            trial_file.write("echo RUNNING TRIAL NUMBER: $SLURM_ARRAY_TASK_ID \n")
            trial_file.write("module load apptainer\n")
            apptainer_confg_line = f"apptainer exec optimal_cut.sif run_trial({path_to_scip_config_toml}, {path_to_SciPy_config_toml},... )"
            trial_file.write(apptainer_confg_line)
        trial_number += 1
    if trial_number+1 != data_queue_length: # pretty sure this is the right number 
        experiments_logger.ERROR(f"Wrote: {trial_number+1} trials\n Had data queue start length : {data_queue_length}")
    os.chdir(os.genenv("EXP_TEMP"))
    with open("run_trial.sh", w) as run_trial:
        run_trial.write("#!/bin/bash\n")
        run_trial.write("source cluster_enviroment.sh\n")
        run_trial.write("echo \"Running on $(hostname)\"\n echo \"TASK: $SLURM_ARRAY_TASK_ID\" \n")
        run_trial.write("chmod +x \"$TRIAL_FILES_TO_RUN\\trial_$SLURM_ARRAY_TASK_ID.sh\"\n")
        run_trial.write(".\\$TRIAL_FILES_TO_RUN\\trial_$SLURM_ARRAY_TASK_ID.sh\"\n")
    with open("trials_env.sh", w) as trial_env:
        trial_env.write("#!/bin/bash\n")
        trial_env.write(f"export NUM_TRIAL={data_queue_length}")
