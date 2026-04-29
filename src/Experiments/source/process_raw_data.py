import logging
import os
import json
import tomllib
from pyscipopt import Model
from parametricCutGen.cluster_utils import parse_logger

#data_processer_logger = parse_logger(__name__)

test_paths = {"trial_temp": "src/Experiments/Test/DataCollectionTest", "trial_metadata":"src/Experiments/Test/DataCollectionTest/TestMetadata","processed_trial_data_target": "src/Experiments/Test/TestDataTarget"}

# at this point all data should be written. Just need solve lps and record data points. 

def process_trial(paths, trial_name):
    data_out = {"prim_sol_value_before_cut":[], "dual_sol_value_before_cut": [], "prim_sol_value_after_cut":[], "dual_sol_value_after_cut": []}
    # 
    for metadata_file in os.listdir(paths["trial_metadata"]):
        print(metadata_file)
        trial_metadata = json.load(open(os.path.join(paths["trial_metadata"], metadata_file), "r"))
        before_cut = Model()
        before_cut.readProblem(trial_metadata["mip_base_path"])
        before_cut.relax()
        before_cut.optimize()
        path_to_combined_file = join_cut_and_mip_in_lp_format(trial_metadata["mip_base_path"], trial_metadata["cut_path"], paths)
        after_cut = Model()
        after_cut.readProblem(path_to_combined_file)
        after_cut.relax()
        after_cut.optimize()
        data_out["prim_sol_value_before_cut"].append(before_cut.get)
        data_out["dual_sol_value_before_cut"].append()
        data_out["prim_sol_value_after_cut"].append()
        data_out["dual_sol_value_after_cut"].append()
        os.remove(path_to_combined_file)
    with open(os.path.join(paths["processed_trial_data_target"], f"trial_{trial_name}.json"), "w") as trial_data:
        json.dump(data_out, trial_data) 

def join_cut_and_mip_in_lp_format(path_to_mip_lp_file, path_to_cut_lp_file, paths):
    mip_file =  open(path_to_mip_lp_file, "r")
    cut_file =  open(path_to_cut_lp_file, "r")
    read_mip_file = mip_file.readlines()
    print(read_mip_file)
    insert_at_index = read_mip_file.index("Subject to\n")
    cut_line = cut_file.readlines()
    read_mip_file.insert(insert_at_index+1, cut_line[0])
    os.path.join(paths["trial_temp"], "temp.lp")
    with open(os.path.join(paths["trial_temp"], "temp.lp"), "w") as temp_mip:
        read_mip_file = "".join(read_mip_file)
        temp_mip.write(read_mip_file)

def test_process_trial():
    test_paths = {"trial_temp": "src/Experiments/Test/DataCollectionTest", "trial_metadata":"src/Experiments/Test/DataCollectionTest/TestMetadata","processed_trial_data_target": "src/Experiments/Test/TestDataTarget"}
    process_trial(test_paths, "test")

test_process_trial()
