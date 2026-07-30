from cutgeneratingfunctionology.igp import *
from parametricCutGen.experimental_utils.pyscipopt_data_collection_events import *
from pyscipopt import Model, quicksum, SCIP_PARAMSETTING, exp, log, sqrt, sin, SCIP_HEURTIMING
import os




model_name = "gen-ip002"
fun = gmic(92136398/212067779)
row = 1
exp_id = "test"
scip_time = 100
sol_path = f"solution_files/solutions/{model_name}/1/{model_name}.sol"

result = {}
model = Model()
model.readProblem(filename=f"model_files/{model_name}.mps")
model.setParam("limits/time", scip_time)
model.setSeparating(SCIP_PARAMSETTING.OFF)
sepa = backTest(exp_id, row, fun, model_name)
model.includeSepa(sepa, "optimal_cut", "optimal cut over space of paramaterized cut generating functions", priority=10000, freq=0)
model.setParam("separating/maxcutsroot", 1)
model.setParam("separating/maxroundsroot", 1)
model.setHeuristics(SCIP_PARAMSETTING.OFF)
heuristic = OracleHeurisitc(sol_path)
model.setPresolve(SCIP_PARAMSETTING.OFF)
model.includeHeur(heuristic, "OracleHeurisitc", "for observing changes in dual bound from cuts", "Y", timingmask=SCIP_HEURTIMING.DURINGLPLOOP)
model.setParam("limits/nodes", 1)
model.optimize()
model.writeMIP(filename=f"TEMP/{model_name}.{exp_id}.{row}.lp")
oracle_sol = model.getSols()
print(oracle_sol)
model_lp_with_cut = Model()
model_lp_with_cut.readProblem(filename=f"TEMP/{model_name}.{exp_id}.{row}.lp")
model_lp_with_cut.relax()
model_lp_with_cut.optimize()
sol = model_lp_with_cut.getSols()
dual_value = model_lp_with_cut.getDualbound()
result[exp_id] = (sol, dual_value)
print(result[exp_id])

