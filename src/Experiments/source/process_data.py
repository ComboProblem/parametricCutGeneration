import logging
import os
import json
import tomllib
from pyscipopt import Model

data_processer_logger = parse_logger(__name__)

paths = {}

# at this point all data should be written. Just need solve lps and record data points. 

def process_trial(paths):
    data_out = {"prim_sol_value_before_cut":[], "dual_sol_value_before_cut": [], "prim_sol_value_after_cut":[], "dual_sol_value_after_cut": []}
    for metadata_file in paths["trial_meta_data"]:
        trial_metadata = json.loads(metadata_file)
        before_cut = Model()
        before_cut.readProblem(trial_metadata["mip_base_path"])
        before_cut.relax()
        before_cut.optimze()
        after_cut = Model()
        after_cut.readProble(trial_metadata["cut_path"])
        after_cut.relax()
        after_cut.optimze()
        data_out["prim_sol_value_before_cut"].append()
        data_out["dual_sol_value_before_cut"].append()
        data_out["prim_sol_value_after_cut"].append()
        data_out["dual_sol_value_after_cut"].append()
    with open(os.path.join(paths["process_trial_data_target"], )) as trial_data:
        data_out = data_out | trial_metadata
        json.dumps(data_out, trial_data) 
        
