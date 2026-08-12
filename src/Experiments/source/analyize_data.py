import os
import json
import pandas



# table format might be helpful? 
# model row cgf cut_score cgp_epsilon max_constraint_viloation is_gmic realitive_gap tree_deth_approx


# cgf_data_table.groupby("is_gmic")["realitive_gap"].mean()
def load_cgfl_to_dataframe(path):
    """
    .cgfl -> 
    """
    with os.open(path, "r") as cgfl:
        

def load_data():
    processed_data = {}
    for model in model_files:
        processed_data[f"{model}"] = json.loads()
    return processed_data

def report_generation(processed_data):    
    return sum( processed_data[model]["number_of_successful_trials"] for model in processed_data.keys() )

# mip
# Analysis. for model, trial, if trial is not gmic and model, trial cgf has small constraint violation, then relative dual gap model, trial!= retaliative dual gap model, gmic
# Analysis. For model, trial, if trial is not gmic, model, trial cgf has small constraint violation, then  tree size estimate < tree size estimate model, gmic
# cgp
# Analysis. Proportion of cgf problems solved that have a solution of gmic.
# Analysis. Proportion of cgf for model, eps, cut_score[i], row !=  model, esp, cut_score[j], row
# Analysis. Proportion of cgf which have small constraint violation. 

