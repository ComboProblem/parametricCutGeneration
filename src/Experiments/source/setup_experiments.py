import logging
import os
from parametricCutGen.utils import validate_paths, parse_logger

set_experiments_logger = parse_logger(__name__)

paths = validate_paths({"container":os.getenv("OPTIMAL_CUT_CONTAINER"), "experiment":os.getenv("EXP_PARAM_PATH"),"model_files":os.getenv("MODEL_FILES")}, "data_out_path":os.getenv("DATA_OUT_FILE_PATH"), ""}, set_experiments_logger)



def setup_parametric_experimenets(paths):
    return experiment_list # list of dictionaries with particular parameters for each experimenet


# this gets ran as the first request. 
def write_trials_for_experiments(paths, experiment):
    """
    Writes all experimental trials to trial directory and defines a scrip to run all trials.
    """
    trial_number = 0
    # while loop for writing trial_$NUMBER.sh = f"trail_{trial_number}.sh"
    for model_file in os.listdir(path["model_files"]):
        trial_file_name = f"trial_{trial_number}.sh"
        with open(os.path.join(paths["trial_scripts"], trial_file_name), w) as trial_file:
            trial_file.write("#!/bin/bash\n")
            trial_file.write("source cluster_enviroment.sh\n")
            trial_file.write(f"echo RUNNING TRIAL NUMBER $SLURM_ARRAY_TASK_ID for experiment {experiment}\n")
            trial_file.write("module load apptainer\n")
            apptainer_confg_line = f"apptainer run {paths["container"]} {os.path.join(paths["experiment"],"run_trial.py"}"
            trial_file.write(apptainer_confg_line)
        trial_number += 1
    with open(os.path.join(paths["trail_scripts"], "run_trials.sh"), w) as run_trials:
        run_trials.write("#!/bin/bash\n")
        run_trials.write("source cluster_enviroment.sh\n")
        run_trials.write("echo \"Running on $(hostname)\"\n echo \"TASK: $SLURM_ARRAY_TASK_ID\" \n")
        run_trials.write("chmod +x \"$TRIAL_FILES_TO_RUN\\trial_$SLURM_ARRAY_TASK_ID.sh\"\n")
        run_trials.write(".\\$TRIAL_FILES_TO_RUN\\trial_$SLURM_ARRAY_TASK_ID.sh\"\n")
    with open(os.path.join(paths[], "trials_env.sh"), w) as trial_env:
        trial_env.write("#!/bin/bash\n")
        trial_env.write(f"export NUM_TRIAL={data_queue_length}")

def __main__():
    for experiment in setup_parametric_experimenets(paths):
        write_trials_for_experiments(paths, experiment)
