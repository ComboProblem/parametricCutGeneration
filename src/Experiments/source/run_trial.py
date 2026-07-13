from cutgeneratingfunctionology.igp import *
from parametricCutGen.pyscipopt_optimal_cut_generation import OptimalCut
from parametricCutGen.experimental_utils.pyscipopt_data_collection_events import *
from pyscipopt import Model
import json, logging, signal, time, os

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
functions_logger = logging.getLogger("cutgeneratingfunctionology.igp.functions")
functions_logger.setLevel(logging.ERROR)

def slurm_handler(signum, frame):
    print('Signal handler called with signal', signum)
    raise Exception("Slurm Time Warning. Dump Data")

signal.signal(signal.SIGUSR1, slurm_handler)
signal.signal(signal.SIGTERM, slurm_handler)
signal.signal(signal.SIGINT, slurm_handler)


model_name = os.getenv("MODEL_NAME")
model_path = f"model_files/{model_name}.mps"
sol_path = f"solution_files/solutions/{model_name}/1/{model_name}.sol"
exp_id = int(os.getenv("SLURM_ARRAY_TASK_ID"))
scip_time = 2*60*60 # 2 hours in seconds
logger.info(f"Model:{model_name}\n Experiment:{exp_id}")

cut_score_names = ['parallelism', 'cut_off_distance', 'violation', 'realitive_violation']
trial_bkpts = [2**i for i in range(1,9)]
trial_cuts = [i for i in range(1,9)]

def sage_rational_to_json(obj):
    if isinstance(obj, sage.rings.rational.Rational):
        return {'__sage.rings.rational.Rational__': True, 'numerator': int(obj.numerator()), 'denominator': int(obj.denominator())}
    raise TypeError(f'Cannot serialize object of {type(obj)}')

#try:
if exp_id == 256:
    model = Model()
    model.readProblem(filename=model_path)
    model.setParam("limits/time", scip_time)
    model.setSeparating(SCIP_PARAMSETTING.OFF)
    model.setHeuristics(SCIP_PARAMSETTING.OFF)
    heuristic = OracleHeurisitc(sol_path)
    model.includeHeur(heuristic, "OracleHeurisitc", "for observing changes in dual bound from cuts", "Y", timingmask=SCIP_HEURTIMING.DURINGLPLOOP)
    model.setPresolve(SCIP_PARAMSETTING.OFF)
    model.setParam("limits/nodes", 1)
    model=record_data(model)
    model.redirectOutput()
    model.optimize()
    model.writeStatisticsJson(filename=f"data/{model_name}.no_cuts.stats.json")
    with open(f"data/{model_name}.no_cuts.txt", 'w') as data_file:
        json.dump(model.data, data_file, default=sage_rational_to_json)
    
elif 0 <= exp_id <= 255:
    bit_string =  '0'*(8-exp_id.bit_length()) + bin(exp_id)[2:]
    cut_score_index = int(bit_string[0:2], 2)
    number_breakpoints_index =  int(bit_string[2:5], 2)
    number_of_cuts_index =  int(bit_string[5:], 2)
    cpg_kwds = {'algorithm':'bkpt_as_param', 'cut_score':cut_score_names[cut_score_index],  'max_num_of_bkpts':trial_bkpts[number_breakpoints_index]}
    numb_cuts = trial_cuts[number_of_cuts_index]
    model = Model()
    model.readProblem(filename=model_path)
    model.setParam("limits/time", scip_time)
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
    model.redirectOutput()
    model.optimize()
    model.writeStatistics(f"data/{model_name}.bkpt_as_param.{cut_score_names[cut_score_index]}.{trial_bkpts[number_breakpoints_index]}.{numb_cuts}.stats.json")
    with open(f"data/{model_name}.bkpt_as_param.{cut_score_names[cut_score_index]}.{trial_bkpts[number_breakpoints_index]}.{numb_cuts}.txt", 'w') as data_file:
        json.dump(model.data, data_file, default=sage_rational_to_json)
    
#elif exp_id > 255:
#    exp_id = exp_id - 255
#    bit_string =  '0'*(6-exp_id.bit_length()) + bin(exp_id)[2:]
#    number_breakpoints_index =  int(bit_string[0:3], 2)
#    number_of_cuts_index =  int(bit_string[3:], 2)
#    cpg_kwds = {'algorithm':'value_poly_lp', "max_bkpt":trial_bkpts[number_breakpoints_index]}
#    numb_cuts =  trial_cuts[number_of_cuts_index]
#    numb_cuts = trial_cuts[number_of_cuts_index]
#    model = Model()
#    model.readProblem(filename=model_path)
#    model.setParam("limits/time", scip_time)
#    model.setSeparating(SCIP_PARAMSETTING.OFF)
#    sepa = OptimalCut(write_cgf_data=True, cgp_kwds=cpg_kwds)
#    model.includeSepa(sepa, "optimal_cut", "optimal cut over space of paramaterized cut generating functions", priority=10000, freq=0)
#    model.setParam("separating/maxcutsroot", 1)
#    model.setParam("separating/maxroundsroot", numb_cuts)
#    model.setHeuristics(SCIP_PARAMSETTING.OFF)
#    heuristic = OracleHeurisitc(sol_path)
#    model.includeHeur(heuristic, "OracleHeurisitc", "for observing changes in dual bound from cuts", "Y", timingmask=SCIP_HEURTIMING.DURINGLPLOOP)
#    model.setPresolve(SCIP_PARAMSETTING.OFF)
#    model.setParam("limits/nodes", 1)
#    model=record_data(model)
#    model.hideOutput()
#    model.optimize()
#    model.writeStatistics(filename=f"data/{model_name}.bkpt_as_param.value_poly_lp.bkpt_as_param.{cut_score_names[cut_score_index]}.{trial_bkpts[number_breakpoints_index]}.{numb_cuts}.out")
#    with open(f"data/{model_name}.value_poly_lp.bkpt_as_param.{cut_score_names[cut_score_index]}.{trial_bkpts[number_breakpoints_index]}.{numb_cuts}.txt", 'w') as data_file:
#        json.dump(model.data, data_file, default=sage_rational_to_json)
#
#except Exception as e:
#    logging.debug(e)
#    print()
#    try:
#        data = {"Exception": str(e)} | model.data
#    except Exception as e2:
#        data = {"Exception": str(e), "failure": str(e2)}
