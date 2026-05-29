# from parametricCutGen.scip_data_collection_events import CutGapDataRecording
from parametricCutGen.optimal_cut_generation import OptimalCut
from parametricCutGen.logging_utils import * # cgp_written_log is here
from parametricCutGen.scip_data_collection_events import record_data
from pyscipopt import Model
from pyscipopt import Model, quicksum, SCIP_PARAMSETTING, exp, log, sqrt, sin

import logging

test_logging = logging.getLogger(__name__)
test_logging.setLevel(logging.DEBUG)


model = Model()
model.readProblem(filename="/home/acadia/Downloads/gen-ip016.mps")
model.setSeparating(SCIP_PARAMSETTING.OFF)
#logging.disable()
sepa = OptimalCut(write_cgf_data=True, cgp_kwds={"algorithm" : "value_poly_lp", "max_num_of_bkpts": 16, "backend":None})
model.includeSepa(sepa, "optimal_cut", "optimal_cut test", priority=10000, freq=0)
model.setHeuristics(SCIP_PARAMSETTING.OFF)
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
model.setParam("branching/random/priority", 20000)
model.setParam("limits/nodes", 10)
model=record_data(model)
try:
    model.optimize()
except Exception as e:
    print(e)
print(model.getStatus())
model.printStatistics()
print(model.data)


