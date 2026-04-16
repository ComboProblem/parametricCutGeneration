import logging
import os
import json
import tomllib
from parametricCutGen.utils import validate_paths, parse_logger

setup_experiments_logger = parse_logger(__name__)

paths = ({"container":os.getenv("OPTIMAL_CUT_CONTAINER"), "experiment_params":os.getenv("EXP_PARAM_PATH"), "model_files":os.getenv("MODEL_FILES")}, "conduct_experiment_base":os.getenv("PARAMETRIC_EXPS_BASE"))
paths = validate_paths(paths, setup_experiments_logger)

def setup_parametric_experiments(paths):
    
    return experiment_list # list of dictionaries with particular parameters for each experiments

def setup_experiment_paths(paths, experiment):
    """
    Experiment paths are named in the following manner $PARAMETRIC_EXPS_BASE/algorithm/cut_score/number_of_bkpts
    where the last bit is an if applicable type of deal. Each experimental path has a folder named TrialPrograms and Data. 
    """
    experiment_dir_algo = os.path.join(paths["conduct_experiment_base"], experiment["algorithm"])
    if not os.path.exists(experiment_dir_algo):
        os.mkdir(experiment_dir_algo)
    experiment_dir_algo_cut_score = os.path.join(experiment_dir_algo, experiment["cut_score"])
    if not os.path.exists(experiment_dir_algo_cut_score):
        os.mkdir(experiment_dir_algo_cut_score)
    experiment_dir_algo_cut_score_num_bkpt = os.path.join(experiment_dir_algo_cut_score, experiment["max_number_of_bkpts"])
    if not os.path.exists(experiment_dir_algo_cut_score_num_bkpt):
        os.mkdir(experiment_dir_algo_cut_score_num_bkpt)
    trial_programs_dir = os.path.join(experiment_dir_algo_cut_score_num_bkpt, "TrialPrograms")
    if not os.path.exists(trial_programs_dir):
        os.mkdir(trial_programs_dir)
    exp_data_dir = os.path.join(experiment_dir_algo_cut_score_num_bkpt, "Data")
    if not os.path.exists(exp_data_dir):
        os.mkdir(exp_data_dir)
    paths["trial_programs"] = trial_programs_dir
    paths["exp_data"] = exp_data_dir
    return paths

def write_trials_for_experiments(paths, experiment):
    """
    Writes all experimental trials to trial directory and defines a scrip to run all trials.
    """
    # setup paths
    paths = setup_experiment_paths(paths, experiment)
    setup_experiments_logger.info(f"Writing config file for {experiment}")
    with open(os.path.join(paths["trial_programs"], "experiment.json")) as exp_config:
        json.dump(experiment, exp_config)
    setup_experiments_logger.info(f"Writing trials for experiment with parameters {experiment}")
    # main loop
    trial_number = 0
    for model_file in os.listdir(paths["model_files"]):
        trial_file_name = f"trial_{trial_number}.sh"
        with open(os.path.join(paths["trial_programs"], trial_file_name), w) as trial_file:
            trial_file.write("#!/bin/bash\n")
            trial_file.write("source cluster_environment.sh\n")
            trial_file.write(f"export MODEL={model_file}\n")
            trial_file.write(f"export DATA_TARGET={paths["exp_data"]}\n")
            trial_file.write(f"echo RUNNING TRIAL NUMBER $SLURM_ARRAY_TASK_ID for experiment {experiment}\n")
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
        run_trials.write(f"sbatch array=0-$NUMBER_OF_TRIALS account=$CLUSTER_ACCOUNT partition=$PARTITION time=$TRIAL_TIME:00 mem=$TRIAL_MEM {os.path.joint(paths["trial_programs"], "trial_$.sh"}") 

def __main__():
    for experiment in setup_parametric_experiments(paths):
        write_trials_for_experiments(paths, experiment)
