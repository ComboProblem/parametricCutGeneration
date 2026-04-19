import logging
import os
import json
import tomllib
from parametricCutGen.utils import validate_paths, parse_logger

setup_experiments_logger = parse_logger(__name__)

_shell_paths = {"container":os.getenv("OPTIMAL_CUT_CONTAINER"), "experiment_params":os.getenv("EXP_PARAM_PATH"), "model_files":os.getenv("MODEL_FILES"), "conduct_experiment_base":os.getenv("PARAMETRIC_EXPS_BASE") }


def setup_parametric_experiments(paths):
    with open(paths["experiment_params"]) as f:
        parsed_toml_file = tomllib.loads(f.read())
    global_settings = {setting:parsed_toml_file['experiment']['settings'][setting]  for setting in parsed_toml_file['experiment']['settings']}
    experiment_list = [ {param:parsed_toml_file['experiment']['parameters'][algo][param] for param in parsed_toml_file['experiment']['parameters'][algo]}  for algo in parsed_toml_file['experiment']['parameters']['algorithm'] ]
    return experiment_list, global_settings # list of dictionaries with particular parameters for each experiments and a dictionary of global settings

def setup_experiment_paths(paths, experiment):
    """
    Experiment paths are named in the following manner $PARAMETRIC_EXPS_BASE/algorithm/cut_score/max_number_of_bkpts/max_number_of_cuts
    where the last bit is an if applicable type of deal. Each experimental path has a folder named TrialPrograms and Data. 
    """
    # we order directory names by ordering the names of the parameters
    experimental_paramaters = ["algorithm", "cut_score", "max_num_of_bkpts", "number_of_cuts"]
    working_dir =  paths["conduct_experiment_base"]
    for param in experimental_paramaters:
        path = os.path.join(working_dir, experiment[param])
        if not os.path.exists(path):
            os.mkdir(path)
        else:
            setup_experiments_logger.debug(f"Path {path} already exists.")
        working_dir = path
    trial_program_dir = os.path.join(working_dir, "TrialPrograms")
    if not os.path.exists(trial_program_dir):
        os.mkdir(trial_program_dir)
    else:
        setup_experiments_logger.warning(f"Trial Program paths exists. Proceed with caution path: {trial_program_dir}")
    exp_data_dir =  os.path.join(working_dir, "Data")
    if not os.path.exists(exp_data_dir):
        os.mkdir(exp_data_dir)
    else:
        setup_experiments_logger.warning(f"Trial Program paths exists. Proceed with caution path: {exp_data_dir}")
    paths["trial_programs"] = trial_program_dir
    paths["exp_data"] = exp_data_dir
    return paths

# potential to extend, for now we'll use a universal number_of_cuts or calls to the paramaterized_problem_solver
def define_trial_scip_settings(paths, experiment):
    raise NotImplementedError
    with open(os.join(paths["trial_programs"], "scip_experiment_parmaters.set")) as scip_exp_param_file:
        # maximal number of separation rounds per node (-1: unlimited)
        # [type: int, advanced: FALSE, range: [-1,2147483647], default: -1]
        scip_exp_param_file.write(f"separating/maxrounds = {experiment[""]}\n")
        # maximal number of separation rounds in the root node (-1: unlimited)
        # [type: int, advanced: FALSE, range: [-1,2147483647], default: -1]
        scip_exp_param_file.write(f"separating/maxroundsroot = {experiment[""]}\n")
        # maximal number of separation rounds in the root node of a subsequent run (-1: unlimited)
        # [type: int, advanced: TRUE, range: [-1,2147483647], default: -1]
        scip_exp_param_file.write(f"separating/maxroundsrootsubrun = {experiment[""]}\n")
        # maximal additional number of separation rounds in subsequent price-and-cut loops (-1: no additional restriction)
        # [type: int, advanced: TRUE, range: [-1,2147483647], default: 1]
        scip_exp_param_file.write(f"separating/maxaddrounds = {experiment[""]}\n")
        # maximal number of consecutive separation rounds without objective or integrality improvement in local nodes (-1: no additional restriction)
        # [type: int, advanced: FALSE, range: [-1,2147483647], default: 1]
        scip_exp_param_file.write(f"separating/maxstallrounds = {experiment[""]}\n")
        # maximal number of consecutive separation rounds without objective or integrality improvement in the root node (-1: no additional restriction)
        # [type: int, advanced: FALSE, range: [-1,2147483647], default: 10]
        scip_exp_param_file.write(f"separating/maxstallroundsroot = {experiment[""]}\n")

def write_trials_for_experiments(paths, experiment):
    """
    Writes all experimental trials programs to trial directory and defines a program to run all trial programs.
    """
    # setup paths
    paths = setup_experiment_paths(paths, experiment, global_settings)
    setup_experiments_logger.info(f"Writing config file for {experiment}")
    with open(os.path.join(paths["trial_programs"], "experiment_parameters.json")) as exp_config:
        json.dump(experiment | global_settings, exp_config)
    setup_experiments_logger.info(f"Writing trials for experiment with parameters {experiment}")
    # main loop
    trial_number = 0
    for model_file in os.listdir(paths["model_files"]):
        trial_file_name = f"trial_{trial_number}.sh"
        with open(os.path.join(paths["trial_programs"], trial_file_name), w) as trial_file:
            trial_file.write("#!/bin/bash\n")
            trial_file.write("source cluster_environment.sh\n")
            trial_file.write(f"export EXPERIMENT_TRIAL_PROGRAMS_PATH={paths["trial_programs"]}\n") # contains experiment_parameters.json
            trial_file.write(f"export MODEL={model_file}\n")
            trial_file.write(f"export DATA_TARGET_PATH={paths["exp_data"]}\n")
            trial_file.write(f"echo RUNNING TRIAL NUMBER $SLURM_ARRAY_TASK_ID for experiment parameters {experiment}\n")
            trial_file.write("module load apptainer\n")
            apptainer_confg_line = f"apptainer run {paths["container"]} {os.path.join(paths["experiment"],"run_trial.py"}"
            trial_file.write(apptainer_confg_line)
        trial_number += 1
    setup_experiments_logger.info(f"{trial_number+1} trials have been written for {experiment}")
    with open(os.path.join(paths["trial_programs"], "run_trials.sh"), w) as run_trials:
        run_trials.write("#!/bin/bash\n")
        run_trials.write("source cluster_enviroment.sh\n")
        run_trials.write("echo \"Running on $(hostname)\"\n echo \"TASK: $SLURM_ARRAY_TASK_ID\" \n")
        run_trials.write("chmod +x \"$TRIAL_FILES_TO_RUN\\trial_$SLURM_ARRAY_TASK_ID.sh\"\n")
        run_trials.write(".\\$TRIAL_FILES_TO_RUN\\trial_$SLURM_ARRAY_TASK_ID.sh\"\n")
        run_trials.write(f"sbatch array=0-$NUMBER_OF_TRIALS account=$CLUSTER_ACCOUNT partition=$PARTITION time=$TRIAL_TIME:00 mem=$TRIAL_MEM {os.path.joint(paths["trial_programs"], "trial_$SLURM_ARRAY_TASK_ID.sh"}") 

def __main__():
    paths = validate_paths(_shell_paths, setup_experiments_logger)
    experiments, global_settings = setup_parametric_experiments(paths)
    for experiment in experiments:
        write_trials_for_experiments(paths, experiment, global_settings)
