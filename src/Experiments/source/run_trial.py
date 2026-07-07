from parametricCutGen.pyscipopt_optimal_cut_generation import OptimalCut
from parametricCutGen.experimental_utils.pyscipopt_data_collection_events import *
from pyscipopt import Model
import json
import os

model_name = os.getenv("PROBLEM")
model_path = os.path.join("model_files", model_name)
sol_path = os.path.join("solution_files", model_name)
exp_id = int(os.getenv("SLURM_ARRAY_TASK_ID"))

cut_score_names = ['parallelism', 'cut_off', 'violation', 'retaliative_violation']
trial_bkpts = [2**i for i in range(1,9)]
trial_cuts = [range(1,9)]

def sage_rational_to_json(obj):
    if isinstance(obj, sage.rings.rational.Rational):
        return {'__sage.rings.rational.Rational__': True, 'numerator': int(obj.numerator()), 'denominator': int(obj.denominator())}
    raise TypeError(f'Cannot serialize object of {type(obj)}')

def as_sage_rational(dct):
    if '__sage.rings.rational.Rational__' in dct:
        return QQ(dct["numerator"]/dct["denominator"])
    return dct

if exp_id == -1:
	model = Model()
	model.readProblem(filename=model_path)
	# model.setParam("limits/time", 3600)
	model.setSeparating(SCIP_PARAMSETTING.OFF)
	model.setHeuristics(SCIP_PARAMSETTING.OFF)
    heuristic = OracleHeurisitc(sol_path)
    model.includeHeur(heuristic, "OracleHeurisitc", "for observing changes in dual bound from cuts", "Y", timingmask=SCIP_HEURTIMING.DURINGLPLOOP)
	model.setPresolve(SCIP_PARAMSETTING.OFF)
    model.setParam("limits/nodes", 1)
	model=record_data(model)
	model.hideOutput()
	model.optimize()
	model.writeStatistics(filename=os.path.join(os.getenv("DATA"),f"{model_name}.no_cuts.out"))
    with open(os.path.join(os.getenv("DATA"),f"{model_name}.no_cuts.txt"), 'w') as data_file:
        json.dump(model.data, data_file, default=sage_rational_to_json)
	
elif 0 <= exp_id <= 255:
	bit_string =  '0'*(8-exp_id.bit_length()) + bin(exp_id)[2:]
	cut_score_index = int(bit_string[0:2], 2)
	number_breakpoints_index =  int(bit_string[2:5], 2)
	number_of_cuts_index =  int(bit_string[5:], 2)
	cpg_kwds = {'algorithm':'bkpt_as_param', 'cut_score':cut_score_names[cut_score_index],  'max_num_of_bkpts': trial_bkpts[number_breakpoints_index]}
	numb_cuts = trial_cuts[number_of_cuts_index]
	model = Model()
	model.readProblem(filename=model_path)
	model.setParam("limits/time", 3600)
	model.setSeparating(SCIP_PARAMSETTING.OFF)
    sepa = OptimalCut(write_cgf_data=True, cgp_kwds=cpg_kwds)
    model.includeSepa(sepa, "optimal_cut", "optimal cut over space of paramaterized cut generating functions", priority=10000, freq=0)
    model.setParam("separating/maxcutsroot", 1)
    model.setParam("separating/maxroundsroot", numb_cuts)
	model.setHeuristics(SCIP_PARAMSETTING.OFF)
    heuristic = OracleHeurisitc(sol_path)
    model.includeHeur(heuristic, "OracleHeurisitc", "for observing changes in dual bound from cuts", "Y", timingmask=SCIP_HEURTIMING.DURINGLPLOOP)
	model.setPresolve(SCIP_PARAMSETTING.OFF)
    model.setParam("limits/nodes", 1)
	model=record_data(model)
	model.hideOutput()	
    model.optimize()
	model.writeStatistics(filename=os.path.join(os.getenv("DATA"),f"{model_name}.bkpt_as_param.{cgp_kwds['cut_score']}.{cgp_kwds[max_num_of_bkpts]}.{numb_cuts}.out"))
    with open(os.path.join(os.getenv("DATA"),f"{model_name}.no_cuts.txt"), 'w') as data_file:
        json.dump(model.data, data_file, default=sage_rational_to_json)
	
elif exp_id > 255:
	exp_id = exp_id - 255
	bit_string =  '0'*(6-exp_id.bit_length()) + bin(exp_id)[2:]
	number_breakpoints_index =  int(bit_string[0:3], 2)
	number_of_cuts_index =  int(bit_string[3:], 2)
	cpg_kwds = {'algorithm':'value_poly_lp', max_bkpt=trial_bkpts[number_breakpoints_index]}
	numb_cuts =  trial_cuts[number_of_cuts_index]
	numb_cuts = trial_cuts[number_of_cuts_index]
	model = Model()
	model.readProblem(filename=model_path)
	model.setParam("limits/time", 3600)
	model.setSeparating(SCIP_PARAMSETTING.OFF)
    sepa = OptimalCut(write_cgf_data=True, cgp_kwds=cpg_kwds)
    model.includeSepa(sepa, "optimal_cut", "optimal cut over space of paramaterized cut generating functions", priority=10000, freq=0)
    model.setParam("separating/maxcutsroot", 1)
    model.setParam("separating/maxroundsroot", numb_cuts)
	model.setHeuristics(SCIP_PARAMSETTING.OFF)
    heuristic = OracleHeurisitc(sol_path)
    model.includeHeur(heuristic, "OracleHeurisitc", "for observing changes in dual bound from cuts", "Y", timingmask=SCIP_HEURTIMING.DURINGLPLOOP)
	model.setPresolve(SCIP_PARAMSETTING.OFF)
    model.setParam("limits/nodes", 1)
	model=record_data(model)
	model.hideOutput()
	model.optimize()
	model.writeStatistics(filename=os.path.join(os.getenv("DATA"),f"{model_name}.bkpt_as_param.{cgp_kwds['cut_score']}.{cgp_kwds[max_num_of_bkpts]}.{numb_cuts}.out"))
    with open(os.path.join(os.getenv("DATA"),f"{model_name}.no_cuts.txt"), 'w') as data_file:
        json.dump(model.data, data_file, default=sage_rational_to_json)
