import os
import json
import pandas as pd
from pathlib import Path
import logging
from cutgeneratingfunctionology.igp import *
from parametricCutGen.cut_generation_problem import inf_norm_of_cont_pwl
from pyscipopt import Model

# metaanalysis and output results. 
# logging levels
# debug level = detailed analyis
# info = meta_analysis only
# warning = poential issues with data
# error = failure

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
fh = logging.FileHandler("cgfl.report")
fh.setLevel(logging.DEBUG)
logger.addHandler(fh)
#ch = logging.StreamHandler()
#ch.setLevel(logging.INFO)
#logger.addHandler(ch)

pcg_exp_logger = logging.getLogger("parametricCutGen.experimental_utils.pyscipopt_data_collection_events")
pcg_exp_logger.setLevel(logging.ERROR)
functions_logger = logging.getLogger("cutgeneratingfunctionology.igp.functions")
functions_logger.setLevel(logging.ERROR)

cut_score_names = ['parallelism', 'cut_off_distance', 'violation', 'realitive_violation']
trial_eps_denom = [2*i for i in range(1,17)]

def as_sage_rational(dct):
    """
    TESTS::
    >>> sage_rational_zero = as_sage_rational({'__sage.rings.rational.Rational__': True, 'numerator': 0, 'denominator': 1})
    >>> sage_rational_zero in QQ
    True
    """
    if "__sage.rings.rational.Rational__" in dct:
        return QQ(dct["numerator"]/dct["denominator"])
    return dct

# The columns of the data frame will look like.
# model row cgf cut_score cgp_epsilon max_constraint_viloation is_gmic dual_value tree_deth_approx scip_efficacy


def load_cgf_data():
    data = {"model_name":[], "pwl1D":[], "row": [], "cut_score_name" : [], "chart_epsilon":[], "problem_dim":[], "max_constraint_violation":[], "is_gmic":[], "dual_bound":[],"tree_depth_approx":[], "scip_efficacy": [], "cgp_version":[], "numerical_epsilon":[]}
    index = 0
    p = Path(".") # should be ran from parametricCutGeneration/src/Experiments
    for cgf_log_file_name in list(p.glob('result/*.cgfl')):
        logger.debug(f"Loading: {cgf_log_file_name}")
        with open(cgf_log_file_name, "r") as cgf_log:
            cgfl_dict = json.load(cgf_log)
            if cgfl_dict["bkpt_val_cgf"][0] is None: # to do, once this set of data is processed, make a fix. This line is due to an old writting of not cuts as [None, None] not [(None,None)]
                cgfl_dict["bkpt_val_cgf"].pop(0)
                cgfl_dict["bkpt_val_cgf"][0] = (None, None)
            for ind in range(len(cgfl_dict["bkpt_val_cgf"])):
                if cgfl_dict["generation_params"][ind]["cut_score"] is not None:
                    for i in range(len(cgfl_dict["generation_params"][ind]["cut_score"])):
                        cgfl_row = {}
                        b = [ as_sage_rational(dct) for dct in cgfl_dict["bkpt_val_cgf"][ind][0][0] ]
                        v = [ as_sage_rational(dct) for dct in cgfl_dict["bkpt_val_cgf"][ind][0][1] ]
                        fun = piecewise_function_from_breakpoints_and_values(b, v)
                        data["model_name"].append(cgfl_dict["metadata"]["model_name"])
                        data["pwl1D"].append(fun)
                        data["row"].append(cgfl_dict["bkpt_val_cgf"][ind][1])
                        data["cut_score_name"].append(cgfl_dict["generation_params"][ind]["cut_score"][i])
                        eps = as_sage_rational(cgfl_dict["generation_params"][ind]["chart_epsilon"][i])
                        data["chart_epsilon"].append(eps)
                        data["problem_dim"].append(cgfl_dict["generation_params"][ind]["problem_dim"][i])
                        data["max_constraint_violation"].append(cgfl_dict["stats"][ind]["max_constraint_violation"])
                        data["is_gmic"].append(cgfl_dict["stats"][ind]["is_gmic"])
                        if cgfl_dict["stats"][ind]["dual_bound"] == "DidNotCompute":
                            data["dual_bound"].append(cgfl_dict["stats"][ind]["dual_bound"])
                        else: 
                            data["dual_bound"].append(float(cgfl_dict["stats"][ind]["dual_bound"]))                            
                        data["tree_depth_approx"].append(cgfl_dict["stats"][ind]["tree_depth_approx"])
                        data["scip_efficacy"].append(cgfl_dict["stats"][ind]["cut_efficacy"])
                        data["cgp_version"].append(cgfl_dict["metadata"]["cgp_version"])
                        data["numerical_epsilon"].append(cgfl_dict["metadata"]["numerical_epsilon"])
        logger.debug(f"Loaded: {cgf_log_file_name}")
    cgf_df = pd.DataFrame(data)
    return cgf_df
        

def report_different_cgf_generation_rate(successful_generation_df, generation_parameter_names=None, fun_equiv=1e-9):
    data = []
    if generation_parameter_names is None:
        generation_rate_name_grouping = ["row", "model_name"]
    else:
        generation_rate_name_grouping = ["row", "model_name"] + generation_parameter_names
    for df in successful_generation_df.groupby(generation_rate_name_grouping):
        not_gmic_df = df[1].query("not `is_gmic`")
        if len(not_gmic_df) > 0: # if 0 -> all cgfs are gmic for the row. We care about non-gmic functions.
            count = 0
            for fun in not_gmic_df["pwl1D"]:
                if all([ inf_norm_of_cont_pwl(fun, pwl1d) > fun_equiv for pwl1d in not_gmic_df["pwl1D"] if inf_norm_of_cont_pwl(fun, pwl1d) != 0]):
                    count +=1
            data.append(float(QQ(count/len(not_gmic_df))))
    return pd.Series(data, name=f"Generation_success_rate_by_{generation_parameter_names}")

def dual_bound_stats_by_model_row_cut_score(successful_generation_df, mininum_generations=5):
    result = {}
    for model_name_group in successful_generation_df.groupby("model_name"):
        model = Model()
        model_name, model_name_df = model_name_group[0], model_name_group[1]
        result[model_name] = {}
        model.readProblem(f"model_files/{model_name}.mps")
        objsen = model.getObjectiveSense()
        logger.debug(f"Checking dual bounds for {model_name}.\nObjective sense: {objsen}")
        for row in model_name_df.groupby("row"):
            if objsen == "maximize":
                op = operator.lt
            else:
                op = operator.gt
            updated_porportions = dual_value_stat(row[1], op, prev_porportions=result[model_name], mininum_generations=mininum_generations)
            if updated_porportions is not None:
                result[model_name] = updated_porportions
        for cut_score_name in result[model_name]:
            if len(result[model_name][cut_score_name]) >= 1:
                result[model_name][cut_score_name] = pd.Series(result[model_name][cut_score_name], name=f"{model_name}_{cut_score_name}_dual_bound_comparision_series")
                logger.debug(f"{result[model_name][cut_score_name].describe()}")
            else:
                result[model_name][cut_score_name] = None
    return result

def dual_value_stat(grouped_by_row_and_model_df, op, prev_porportions = {}, mininum_generations=5):
    porportions = prev_porportions.copy()
    for grouped in  grouped_by_row_and_model_df.groupby(["cut_score_name"]):
        cut_score_name, df = grouped[0][0], grouped[1]
        if cut_score_name not in porportions:
            porportions[cut_score_name] = []
        esp_2_df = df[df["chart_epsilon"]==1/2]
        if len(esp_2_df) > 1:
            logger.debug(f"dual_value_stat: {len(esp_2_df)} solutions found for chart_epsilon==1/2. ")
        for esp_2_dual_bound in esp_2_df["dual_bound"]:
            boolian_data = [True if op(dual_value, esp_2_dual_bound) else False for dual_value in df["dual_bound"]]
            if len(boolian_data) >= mininum_generations:
                proportion = float(QQ((boolian_data.count(True))/(len(boolian_data))))
                porportions[cut_score_name].append(proportion)
            else:
                return None
    return porportions

def cut_efficacy_by_model_and_cut_score(successful_generation_df, mininum_generations=5):
    result = {}
    for grouped in successful_generation_df.groupby(["model_name","row","cut_score_name"]):
        model_name = grouped[0][0]
        cut_score_name = grouped[0][2]
        if model_name not in result:
            result[model_name] = {}
        if cut_score_name not in result[model_name]:
            result[model_name][cut_score_name] = []
        df = grouped[1]
        if len(df) >= mininum_generations:
            esp_2_df = df[df["chart_epsilon"]==1/2]
            if len(esp_2_df) > 1:
                logger.debug(f"dual_value_stat: {len(esp_2_df)} solutions found for chart_epsilon==1/2.")
            for esp_2_cut_efficacy in esp_2_df["scip_efficacy"]:
                more_effective_df = df[df["scip_efficacy"] >= esp_2_cut_efficacy]
                proportion = float(QQ((len(more_effective_df) -1)/(len(df)-1)))
                result[model_name][cut_score_name].append(proportion)
    for model_name in result:
        for cut_score_name in result[model_name]:
            if len(result[model_name][cut_score_name]) >= 1:
                result[model_name][cut_score_name] = pd.Series(result[model_name][cut_score_name], name=f"{model_name}, {cut_score_name} series for cut efficacy")
                logger.debug(f"{result[model_name][cut_score_name].describe()}")
            else:
                result[model_name][cut_score_name] = None
    return result

def __main__():
    result_dir = os.path.join(os.getcwd(), "result")
    logger.info(f"Loading .cgfl files from directory: {result_dir}")
    cgf_df = load_cgf_data()
    logger.info("Cut generation data loaded.")
    successful_generation_df = cgf_df[cgf_df["max_constraint_violation"] < cgf_df["numerical_epsilon"]]
    all_models_cgf_gen_rate = float(QQ(len(successful_generation_df)/len(cgf_df)))
    logger.info(f"Successful cut generation rate across all models: {all_models_cgf_gen_rate:%}")
    gmic_df = successful_generation_df[successful_generation_df["is_gmic"]]
    all_models_gmic_gen_rate = float(QQ(len(gmic_df)/len(cgf_df)))
    logger.info(f"gmic generation rate across all models: {all_models_gmic_gen_rate:%}")
    for generation_parameter in [None, "cut_score_name", "problem_dim"]:
        if generation_parameter is None:
            diff_series = report_different_cgf_generation_rate(successful_generation_df, None)
            logger.info(f"Generating different functions success rate accorss all models and row:\n {diff_series.describe()}")
        else:
            diff_series = report_different_cgf_generation_rate(successful_generation_df, [generation_parameter])
            logger.info(f"Generating different functions success rate accorss all models, row, {generation_parameter}:\n {diff_series.describe()}")
    successful_generation_with_stats_df = successful_generation_df[successful_generation_df["dual_bound"] != "DidNotCompute"]
    dual_bound_series_dict = dual_bound_stats_by_model_row_cut_score(successful_generation_with_stats_df)
    efficacy_series_dict = cut_efficacy_by_model_and_cut_score(successful_generation_with_stats_df)
    models_to_look_at_dual = set()
    models_to_look_at_efficacy = set()
    logger.info("Dual bound calcuations complete; begining analsis.")
    dual_bound_cut_score_mean = {cut_score_name:[] for cut_score_name in cut_score_names}
    for model_name in dual_bound_series_dict:
        for cut_score_name in dual_bound_series_dict[model_name]:
            if dual_bound_series_dict[model_name][cut_score_name] is not None:
                mean = dual_bound_series_dict[model_name][cut_score_name].mean()
                dual_bound_cut_score_mean[cut_score_name].append(mean)
                if mean >= .5:
                    models_to_look_at_dual.add(model_name)
    for cut_score_name in cut_score_names:
        dual_bound_cut_score_mean[cut_score_name] = pd.Series(dual_bound_cut_score_mean[cut_score_name], name=f"Found dual bound imporvement over eps =1/2 by using objective function {cut_score_name}")
        logger.info(f"{dual_bound_cut_score_mean[cut_score_name].describe()}")
    logger.info(f"Names of models to look at for dual bound improvement: {models_to_look_at_dual}")
    efficacy_cut_score_mean ={cut_score_name:[] for cut_score_name in cut_score_names}
    for model_name in efficacy_series_dict:
        for cut_score_name in efficacy_series_dict[model_name]:
            if efficacy_series_dict[model_name][cut_score_name] is not None:
                mean = efficacy_series_dict[model_name][cut_score_name].mean()
                efficacy_cut_score_mean[cut_score_name].append(mean)
                if mean >= .5:
                    models_to_look_at_efficacy.add(model_name)
    logger.info(f"Names of models to look at for efficacy improvement: {models_to_look_at_efficacy}")
    for cut_score_name in cut_score_names:
        efficacy_cut_score_mean[cut_score_name] = pd.Series(efficacy_cut_score_mean[cut_score_name], name=f"Found efficacy of cut imporvement over eps =1/2 by using objective function {cut_score_name}")
        logger.info(f"{efficacy_cut_score_mean[cut_score_name].describe()}")
        
__main__()
