from pyscipopt import Model, SCIP_EVENTTYPE, SCIP_RESULT, Eventhdlr, SCIP_PARAMSETTING
from pyscipopt import Model, Heur, SCIP_HEURTIMING, SCIP_LPSOLSTAT
import os
import logging

data_collection_logger = logging.getLogger(__name__)
data_collection_logger.setLevel(logging.DEBUG)

def record_data(model: Model):
    """
    Attaches an event handler to a given SCIP model that collects primal and dual solutions,
    along with the solving time when they were found.
    The data is saved in model.data["primal_log"] and model.data["dual_log"]. They consist of
    a list of tuples, each tuple containing the solving time and the corresponding solution.

    A usage example can be found in examples/finished/plot_primal_dual_evolution.py. The
    example takes the information provided by this recipe and uses it to plot the evolution
    of the dual and primal bounds over time. 
    """
    class GapEventhdlr(Eventhdlr):
        def eventinit(self): # we want to collect best primal solutions and best dual solutions
            self.model.catchEvent(SCIP_EVENTTYPE.BESTSOLFOUND, self)
            self.model.catchEvent(SCIP_EVENTTYPE.DUALBOUNDIMPROVED, self)

        def eventexec(self, event):
            # if a new best primal solution was found, we save when it was found and also its objective
            if event.getType() == SCIP_EVENTTYPE.BESTSOLFOUND:
                self.model.data["primal_log"].append([self.model.getSolvingTime(), self.model.getPrimalbound()])
            
            if event.getType() == SCIP_EVENTTYPE.DUALBOUNDIMPROVED:
                self.model.data["dual_log"].append([self.model.getSolvingTime(), self.model.getDualbound()])
             

    if not hasattr(model, "data") or model.data==None:
        model.data = {}

    model.data.update({
            'primal_log': [],
            'dual_log': [],
            'cgf_log': [],
            })

    gap_hdlr = GapEventhdlr()
    model.includeEventhdlr(gap_hdlr, "gapEventHandler", "Event handler which collects primal and dual solution evolution")
    return model


# model.readSolFile()
class OracleHeurisitc(Heur):
    def __init__(self, path):
        self._sol = self.model.readSolFile(path)

    def heurexec(self, heurtiming, nodeinfeasible):

        scip = self.model
        result = SCIP_RESULT.DIDNOTRUN
 
        # This heuristic does not run if the LP status is not optimal
        lpsolstat = scip.getLPSolstat()
        if lpsolstat != SCIP_LPSOLSTAT.OPTIMAL:
            return {"result": result}

        # We haven't added handling of implicit integers to this heuristic
        if scip.getNImplVars() > 0:
            return {"result": result}

        # Get the current branching candidate, i.e., the current fractional variables with integer requirements
        branch_cands, branch_cand_sols, branch_cand_fracs, ncands, npriocands, nimplcands = scip.getLPBranchCands()

        # Ignore if there are no branching candidates
        if ncands == 0:
            return {"result": result}

        # Create a solution that is initialised to the LP values
 
        stored = scip.trySol(sol)

        if stored:
            return {"result": SCIP_RESULT.FOUNDSOL}
        else:
            return {"result": SCIP_RESULT.DIDNOTFIND}


import numpy as np

class GreedySolBranch(Branchrule):

    def __init__(self, scip, path):
        self.scip = scip
        self._sol = self.model.readSolFile(path)

    def branchexeclp(self, allowaddcons):

        # Get the branching candidates. Only consider the number of priority candidates (they are sorted to be first)
        # The implicit integer candidates in general shouldn't be branched on. Unless specified by the user
        # npriocands and ncands are the same (npriocands are variables that have been designated as priorities)
        branch_cands, branch_cand_sols, branch_cand_fracs, ncands, npriocands, nimplcands = self.scip.getLPBranchCands()

        # Find the variable that is most fractional
        best_cand_idx = 0
        best_dist = np.inf
        for i in range(npriocands):
            if abs(branch_cand_fracs[i] - 0.5) <= best_dist:
                best_dist = abs(branch_cand_fracs[i] - 0.5)
                best_cand_idx = i

        # Branch on the variable with the largest score
        down_child, eq_child, up_child = self.model.branchVarVal(
            branch_cands[best_cand_idx], branch_cand_sols[best_cand_idx])

        return {"result": SCIP_RESULT.BRANCHED}

#
#def OracleHeurisitc(Heur):
#
#    def heurexec(self, heurtiming, nodeinfeasible):
#
#        scip = self.model
#        result = SCIP_RESULT.DIDNOTRUN
#
#        # This heuristic does not run if the LP status is not optimal
#        lpsolstat = scip.getLPSolstat()
#        if lpsolstat != SCIP_LPSOLSTAT.OPTIMAL:
#            return {"result": result}
#
#        # We haven't added handling of implicit integers to this heuristic
#        if scip.getNImplVars() > 0:
#            return {"result": result}
#
#        # Get the current branching candidate, i.e., the current fractional variables with integer requirements
#        branch_cands, branch_cand_sols, branch_cand_fracs, ncands, npriocands, nimplcands = scip.getLPBranchCands()
#
#        # Ignore if there are no branching candidates
#        if ncands == 0:
#            return {"result": result}
#
#        # Create a solution that is initialised to the LP values
#        sol = scip.createSol(self, initlp=True)
#
#        # Now round the variables that can be rounded
#        for i in range(ncands):
#            old_sol_val = branch_cand_sols[i]
#            scip_var = branch_cands[i]
#            may_round_up = scip_var.varMayRound(direction="up")
#            may_round_down = scip_var.varMayRound(direction="down")
#            # If we can round in both directions then round in objective function direction
#            if may_round_up and may_round_down:
#                if scip_var.getObj() >= 0.0:
#                    new_sol_val = scip.feasFloor(old_sol_val)
#                else:
#                    new_sol_val = scip.feasCeil(old_sol_val)
#            elif may_round_down:
#                new_sol_val = scip.feasFloor(old_sol_val)
#            elif may_round_up:
#                new_sol_val = scip.feasCeil(old_sol_val)
#            else:
#                # The variable cannot be rounded. The heuristic will fail.
#                continue
#
#            # Set the rounded new solution value
#            scip.setSolVal(sol, scip_var, new_sol_val)
#
#        # Now try the solution. Note: This will free the solution afterwards by default.
#        stored = scip.trySol(sol)
#
#        if stored:
#            return {"result": SCIP_RESULT.FOUNDSOL}
#        else:
#            return {"result": SCIP_RESULT.DIDNOTFIND}

#heuristic = OracleHeurisitc()
#scip.includeHeur(heuristic, "OracleHeurisitc", "for observing changes in dual bound from cuts", "Y",
#                 timingmask=SCIP_HEURTIMING.BEFORENODE)
