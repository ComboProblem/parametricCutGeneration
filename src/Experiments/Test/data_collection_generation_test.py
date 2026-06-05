# from parametricCutGen.scip_data_collection_events import CutGapDataRecording
from parametricCutGen.pyscipopt_optimal_cut_generation import OptimalCut
from parametricCutGen.experimental_utils.pyscipopt_data_collection_events import *
from parametricCutGen.cgf_specializations import *
from pyscipopt import Model
from pyscipopt import Model, quicksum, SCIP_PARAMSETTING, exp, log, sqrt, sin, SCIP_HEURTIMING
from cvxpy import Variable, Problem, Minimize

import logging

test_logging = logging.getLogger(__name__)
test_logging.setLevel(logging.DEBUG)

def plot_primal_dual_evolution(model: Model, num_bkpt):
    try:
        from matplotlib import pyplot as plt
    except ImportError:
        raise ImportError("matplotlib is required to plot the solution. Try running `pip install matplotlib` in the command line.\
                          You may also need to install PyQt6 to show the plot.")

    assert model.data["primal_log"], "Could not find any feasible solutions"
    time_primal, val_primal = map(list,zip(*model.data["primal_log"]))
    time_dual, val_dual = map(list,zip(*model.data["dual_log"]))

    
#    if time_primal[-1] < time_dual[-1]:
#        time_primal.append(time_dual[-1])
#        val_primal.append(val_primal[-1])
#
    if time_primal[-1] > time_dual[-1]:
        time_dual.append(time_primal[-1])
        val_dual.append(val_dual[-1])
        
#    plt.plot(time_primal, val_primal, label="Primal bound")
    plt.plot(time_dual, val_dual, label="Dual bound")
    plt.title("dual evolution - gen-ip016.mps")
    plt.suptitle(f"At most k={num_bkpt} breakpoints, 100000 nodes" )
    plt.legend(loc="best")
    return plt

def check_for_gmics(paramed_cgfs):
    pass
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
# think above levels and which primal heursitics are allowed; like rounding ect to find feasible solutions for relaxations proving relaxations ects. 
# What is the parameters; 
data = {}
final_stats = {}
model = Model()
model.readProblem(filename="/home/acadia/Downloads/gen-ip016.mps")
model.setParam("limits/time", 300)
model.setSeparating(SCIP_PARAMSETTING.OFF)
#logging.disable()
#sepa = OptimalCut(write_cgf_data=True, cgp_kwds={"algorithm" : "value_poly_lp", "max_num_of_bkpts": k, "backend":None})
#model.includeSepa(sepa, "optimal_cut", "optimal_cut test", priority=10000, freq=0)
model.setHeuristics(SCIP_PARAMSETTING.OFF)
model.setPresolve(SCIP_PARAMSETTING.OFF)
model.disablePropagation()
#model.setParam("separating/maxcutsroot", 1)
#model.setParam("separating/maxroundsroot", l)
#model.setParam("separating/poolfreq", 0)
#model.setParam("separating/maxroundsrootsubrun", -1)
#model.setParam("separating/maxcuts", 1)
#model.setParam("separating/maxcutsgenfactor", 1)
#model.setParam("separating/maxcutsrootgenfactor", -1)
#model.setParam("separating/maxrounds", 1)
#model.setParam("separating/maxruns", 1)
#sol = model.readSolFile("/home/acadia/Downloads/markshare_4_0.sol.gz")
heuristic = OracleHeurisitc('/home/acadia/Downloads/gen-ip016.sol')
model.includeHeur(heuristic, "OracleHeurisitc", "for observing changes in dual bound from cuts", "Y", timingmask=SCIP_HEURTIMING.DURINGLPLOOP)
#model.setParam("branching/random/priority", 20000)
#model.setParam("limits/nodes", 100000)
model=record_data(model)
model.hideOutput()
try:
    model.optimize()
except Exception as e:
    print(e)
print(model.getStatus())
data["default"] = model.data
model.writeStatistics(filename=f"gen-ip016.bkpt_as_param.steepest_dir.default.out")
default_solve_time = model.getSolvingTime()
default_num_nodes = model.getNTotalNodes()
final_stats["default"] = (default_solve_time, default_num_nodes)

#for l in [1,2,3,4,5]:
#    for k in [2**i for i in range(1,5)]:
#        model = Model()
#        model.readProblem(filename="/home/acadia/Downloads/gen-ip016.mps")
#        model.setSeparating(SCIP_PARAMSETTING.OFF)
#        model.setParam("limits/time", 300)
#        #logging.disable()
#        sepa = OptimalCut(write_cgf_data=True, cgp_kwds={"algorithm" : "bkpt_as_param", "cut_score":"parallelism", "max_num_of_bkpts": k, "backend":None})
#        model.includeSepa(sepa, "optimal_cut", "optimal_cut test", priority=10000, freq=0)
#        model.setHeuristics(SCIP_PARAMSETTING.OFF)
#        model.setPresolve(SCIP_PARAMSETTING.OFF)
#        model.disablePropagation()
#        model.setParam("separating/maxcutsroot", 1)
#        model.setParam("separating/maxroundsroot", l)
#        #model.setParam("separating/poolfreq", 0)
#        #model.setParam("separating/maxroundsrootsubrun", -1)
#        #model.setParam("separating/maxcuts", 1)
#        #model.setParam("separating/maxcutsgenfactor", 1)
#        #model.setParam("separating/maxcutsrootgenfactor", -1)
#        #model.setParam("separating/maxrounds", 1)
#        #model.setParam("separating/maxruns", 1)
#        #sol = model.readSolFile("/home/acadia/Downloads/markshare_4_0.sol.gz")
#        heuristic = OracleHeurisitc('/home/acadia/Downloads/gen-ip016.sol')
#        model.includeHeur(heuristic, "OracleHeurisitc", "for observing changes in dual bound from cuts", "Y", timingmask=SCIP_HEURTIMING.DURINGLPLOOP)
#        #model.setParam("branching/random/priority", 20000)
#        #model.setParam("limits/nodes", 100000)
#        model=record_data(model)
#        model.hideOutput()
#        try:
#            model.optimize()
#        except Exception as e:
#            print(e)
#        print(model.getStatus())
#        trial_solve_time = model.getSolvingTime()
#        trial_num_nodes = model.getNTotalNodes()
#        # model.printStatistics()
#        #model.redirectOutput() # to python
#        model.writeStatistics(filename=f"gen-ip016.bkpt_as_param.steepest_dir.bkpts.{k}.cuts.{l}.out")
#        data[(k,l)] = model.data
#        final_stats[(k,l)] = (trial_solve_time, trial_num_nodes)
#
#        if trial_solve_time < default_solve_time < 300:
#            print(f"bkpts:{k},no_cuts:{l}, decrease solve time by {(default_solve_time-trial_solve_time)/default_solve_time *100}%")
#            if trial_num_nodes < default_num_nodes:
#                print(f"bkpts:{k},no_cuts:{l}, decrease nodes time by {float((default_num_nodes-trial_num_nodes))/default_num_nodes *100}%")
        #plt = plot_primal_dual_evolution(model, k)
        #plt.show()
#print(data)
print(final_stats)

