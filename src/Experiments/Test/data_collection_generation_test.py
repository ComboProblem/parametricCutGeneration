# from parametricCutGen.scip_data_collection_events import CutGapDataRecording
from parametricCutGen.optimal_cut_generation import OptimalCut
from parametricCutGen.logging_utils import * # cgp_written_log is here
from parametricCutGen.scip_data_collection_events import record_data
from parametricCutGen.cgf_specializations import *
from pyscipopt import Model
from pyscipopt import Model, quicksum, SCIP_PARAMSETTING, exp, log, sqrt, sin
from cvxpy import Variable, Problem, Minimize

import logging

test_logging = logging.getLogger(__name__)
test_logging.setLevel(logging.DEBUG)



#b = [0, 76/33672873, 7290191/119903003, 16965671/130581379, 12570730/62300759, 358957910381105/745783659626811, 58934854246777/77490406121583, 3825777230473888/4397050191231867, 3792006029277104/4037478592337619]
#n = len(b)
#x = Variable(n)
#value_cons = value_nnc_polyhedron_constraints(to_sage_rationals(b), n-1, x, coeff_type='float')
#binvc = [1.1 for i in range(n)]
#cut_score_weights = expression_of_steepest_direction_score(binvc, to_sage_rationals(b), n-1, x)
#obj = Minimize(cut_score_weights)
#prob = Problem(obj, value_cons)
#prob.solve()
#print(x.value)
#print(prob.value)
#

model = Model()
model.readProblem(filename="/home/acadia/Downloads/gen-ip016.mps")
model.setSeparating(SCIP_PARAMSETTING.OFF)
#logging.disable()
sepa = OptimalCut(write_cgf_data=True, cgp_kwds={"algorithm" : "value_poly_lp", "max_num_of_bkpts": 100, "backend":"pplite"})
model.includeSepa(sepa, "optimal_cut", "optimal_cut test", priority=10000, freq=0)
#model.setHeuristics(SCIP_PARAMSETTING.OFF)
model.setPresolve(SCIP_PARAMSETTING.OFF)
#model.disablePropagation()
model.setParam("separating/maxcutsroot", 1)
model.setParam("separating/maxroundsroot", 1)
#model.setParam("separating/poolfreq", 0)
#model.setParam("separating/maxroundsrootsubrun", -1)
#model.setParam("separating/maxcuts", 1)
#model.setParam("separating/maxcutsgenfactor", 1)
#model.setParam("separating/maxcutsrootgenfactor", -1)
#model.setParam("separating/maxrounds", 1)
#model.setParam("separating/maxruns", 1)
#model.setParam("branching/random/priority", 20000)
model.setParam("limits/nodes", 100000)
model=record_data(model)
#model.hideOutput()
try:
    model.optimize()
except Exception as e:
    print(e)
print(model.getStatus())
model.printStatistics()
#print(model.data)


