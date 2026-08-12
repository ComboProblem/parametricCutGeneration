from pyscipopt import Model, SCIP_EVENTTYPE, SCIP_RESULT, Eventhdlr, SCIP_PARAMSETTING
from pyscipopt import Model, Heur, SCIP_HEURTIMING, SCIP_LPSOLSTAT
from pyscipopt import Model, Sepa, SCIP_RESULT
from cutgeneratingfunctionology.igp import *
import os
import logging
from contextlib import redirect_stdout
import sys
import io


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
                cuts_applied = []
                for row in self.model.getLPRowsData():
                    if row.getOrigintype() == 3: # created from a separator
                        with redirect_stdout(io.StringIO()) as f:
                            self.model.printRow(row)
                        s = f.getvalue()
                        cuts_applied.append(s[:-2]) # get rid of new line character.
                self.model.data["dual_log"].append([self.model.getSolvingTime(), self.model.getDualbound(), self.model.getNCutsApplied(), cuts_applied])
                
             

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
        self._path = path

    def heurexec(self, heurtiming, nodeinfeasible):

        scip = self.model
        result = SCIP_RESULT.DIDNOTRUN
        sol = scip.readSolFile(self._path)
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


class OracleSelector(Eventhdlr):

    def __init__(self, model,  path):
        Eventhdlr.__init__(model)
        self._path = path

    def eventinit(self):
        self.model.catchEvent(SCIP_EVENTTYPE.NODEFOCUSED, self)
    
    def eventexec(self, event):
        sol = self.model.readSolFile(self._path)
        oracle_val =  self.model.getSolVal(sol, self.model.getObjective())
        if  self.model.getCurrentNode().getLowerbound() > oracle_val:
            self.model.cutoffNode(self.model.getCurrentNode())




class backTest(Sepa):
    def __init__(self, count=None, row=None, fun=None, model_name=None):
        """
        TESTS::
        >>> from parametricCutGen.optimal_cut_generation import OptimalCut
        >>> from piscipopt import Model
        >>> model = Model()
        >>> sepa = OptimalCut() # gmic by default
        >>> model.includeSepa(sepa, "optimal_cut", "full space gmic", priority=1000, freq=1)
        >>> cgp_kwds = {"algorithm":"bkpt_as_param"}
        >>> other_sepa = OptimalCut(cgp_kwds=cgp_kwds)
        >>> model.includeSepa(other_sepa, "optima_cut", "bkpt as param gmic", priority=1000, freq=1)
        """
        self.ncuts = 0
        self.count = count
        self.row=row
        self.fun=fun
        self.model_name=model_name

    def genCutFromCGF(self, cols, rows, binvrow, binvarow, primsol, cgf=None):
        """ Given the row (binvarow, binvrow) of the tableau, computes optimized cut.

        :param primsol:  is the rhs of the tableau row
        :param cols:     are the variables
        :param rows:     are the slack variables
        :param binvrow:  components of the tableau row associated to the basis inverse
        :param binvarow: components of the tableau row associated to the basis inverse * A
        :param 1dPWL:    a minimal function, ab element of PiMin<=k with pi_p(primsol)=1
    
        The intersection cut is given by
         sum(pi_p(a_j) x_j, j in J_I) \geq 1
        where J_I are the integer non-basic variables and J_C are the continuous.
        f_0 is the fractional part of primsol
        a_j is the j-th coefficient of the row and f_j its fractional part
        Note: we create -% <= -f_0 !!
        Note: this formula is valid for a problem of the form Ax = b, x>= 0. Since we do not have
        such problem structure in general, we have to (implicitly) transform whatever we are given
        to that form. Specifically, non-basic variables at their lower bound are shifted so that the lower
        bound is 0 and non-basic at their upper bound are complemented.
        """
        # initialize
        cutcoefs_fun = [0] * len(cols)
        cutrhs_fun = 0
        if cgf is None:
            return cutcoefs_fun, cutrhs_fun
        pi_p = cgf
        # get scip
        scip = self.model
        f = fractional(QQ(primsol))
        def psi(x):
            if x < 0:
                slope = pi_p.functions()[-1]._slope
                neg_part = FastLinearFunction(slope, 0)
                return neg_part(x)
            else:
                return pi_p.functions()[0](x)
        # rhs of the cut is the fractional part of the LP solution for the basic variable
        cutrhs_fun = -f
        # Generate cut coefficients for the original variables
        for c in range(len(cols)):
            col = cols[c]
            assert col is not None
            status = col.getBasisStatus()

            # Get simplex tableau coefficient
            if status == "lower":
                # Take coefficient if nonbasic at lower bound
                rowelem = binvarow[c]
            elif status == "upper":
                # Flip coefficient if nonbasic at upper bound: x --> u - x
                rowelem = -binvarow[c]
            else:
                # variable is nonbasic free at zero -> cut coefficient is zero, skip OR
                # variable is basic, skip
                assert status == "zero" or status == "basic"
                continue

            # Integer variables
            if col.isIntegral():
                cutelem_fun = fractional(QQ(rowelem))
                cutelem_fun = float(-1*f*pi_p(cutelem_fun))
            else:
                # Continuous variables
                cutelem_fun = float(-1*f*psi(QQ(rowelem)))

            # cut is define when variables are in [0, infty). Translate to general bounds
            if not scip.isZero(cutelem_fun):
                if col.getBasisStatus() == "upper":
                    cutelem_fun = -cutelem_fun
                    cutrhs_fun += cutelem_fun * col.getUb()
                else:
                    cutrhs_fun += cutelem_fun * col.getLb()
                # Add coefficient to cut in dense form
                cutcoefs_fun[col.getLPPos()] = cutelem_fun

        # Generate cut coefficients for the slack variables; skip basic ones
        for c in range(len(rows)):
            row = rows[c]
            assert row != None
            status = row.getBasisStatus()

            # free slack variable shouldn't appear
            assert status != "zero"

            # Get simplex tableau coefficient
            if status == "lower":
                # Take coefficient if nonbasic at lower bound
                rowelem = binvrow[row.getLPPos()]
                # But if this is a >= or ranged constraint at the lower bound, we have to flip the row element
                if not scip.isInfinity(-row.getLhs()):
                    rowelem = -rowelem
            elif status == "upper":
                # Take element if nonbasic at upper bound - see notes at beginning of file: only nonpositive slack variables
                # can be nonbasic at upper, therefore they should be flipped twice and we can take the element directly.
                rowelem = binvrow[row.getLPPos()]
            else:
                assert status == "basic"
                continue

            # if row is integral we can strengthen the cut coefficient
            if row.isIntegral() and not row.isModifiable():
                # warning: because of numerics cutelem < 0 is possible (though the fractional part is, mathematically, always positive)
                # However, when cutelem < 0 it is also very close to 0, enough that isZero(cutelem) is true (see later)
                cutelem_fun = fractional(QQ(rowelem))
                cutelem_fun = float(-1*f*pi_p(cutelem_fun))
            else:
                # Continuous variables
                cutelem_fun = float(-1*f*psi(QQ(rowelem)))

            # cut is define in original variables, so we replace slack by its definition
            if not scip.isZero(cutelem_fun):
                # get lhs/rhs
                rlhs = row.getLhs()
                rrhs = row.getRhs()
                assert scip.isLE(rlhs, rrhs)
                assert not scip.isInfinity(rlhs) or not scip.isInfinity(rrhs)

                # If the slack variable is fixed, we can ignore this cut coefficient
                if scip.isFeasZero(rrhs - rlhs):
                  continue

                # Unflip slack variable and adjust rhs if necessary: row at lower means the slack variable is at its upper bound.
                # Since SCIP adds +1 slacks, this can only happen when constraints have a finite lhs
                if row.getBasisStatus() == "lower":
                    assert not scip.isInfinity(-rlhs)
                    cutelem_fun = -cutelem_fun

                rowcols = row.getCols()
                rowvals = row.getVals()

                assert len(rowcols) == len(rowvals)

                # Eliminate slack variable: rowcols is sorted: [columns in LP, columns not in LP]
                for i in range(row.getNLPNonz()):
                    cutcoefs_fun[rowcols[i].getLPPos()] -= cutelem_fun * rowvals[i]

                act = scip.getRowLPActivity(row)
                rhsslack = rrhs - act
                if scip.isFeasZero(rhsslack):
                    assert row.getBasisStatus() == "upper" # cutelem != 0 and row active at upper bound -> slack at lower, row at upper
                    cutrhs_fun -= cutelem_fun * (rrhs - row.getConstant())
                else:
                    assert scip.isFeasZero(act - rlhs)
                    cutrhs_fun -= cutelem_fun * (rlhs - row.getConstant())

        return cutcoefs_fun, cutrhs_fun

    def sepaexeclp(self):
        result = SCIP_RESULT.DIDNOTRUN
        model = self.model
        if not model.isLPSolBasic():
            return {"result": result}

        # get LP data
        cols = model.getLPColsData()
        rows = model.getLPRowsData()

        # exit if LP is trivial
        if len(cols) == 0 or len(rows) == 0:
            return {"result": result}

        result = SCIP_RESULT.DIDNOTFIND

        # get basis indices
        basisind = model.getLPBasisInd()
        # For all basic columns (not slacks) belonging to integer variables, decide if the row should be tried. 

        if self.row is None:
            data_collection_logger.debug(f"None Branch")
            model.writeLP(filename=f"TEMP/{self.model_name}/{self.count}.lp")
            model.interruptSolve()  
            return {"result": result}
        else:
            cgf = self.fun

            #data_collection_logger.debug(f"{fractional(QQ(primsol))}, {cgf}")
            i = None
            for ind in range(len(rows)):
                c = basisind[ind]
                if c == self.row:
                    i = ind
                    break
                elif -c-1 == self.row:
                    i = ind
                    break
            primsol = cols[c].getPrimsol()
            if i is None:
                raise ValueError(f"row:{self.row} is not in LP relaxation")

            # i = [basisind[i] if basisind[i] >= 0 else -basisind[i] -1 for i in range(len(rows))].index(c)	


            # get the row of B^-1 for this basic integer variable with fractional solution value
            binvrow = model.getLPBInvRow(i)

            # get the tableau row for this basic integer variable with fractional solution value
            binvarow = model.getLPBInvARow(i)

            # get current reduced costs for objective evaluation.
            costs = [model.getColRedCost(j) for j in cols if j not in basisind]

            cutcoefs, cutrhs = self.genCutFromCGF(cols, rows, binvrow, binvarow, primsol, cgf)

            cut = model.createEmptyRowSepa(self, "optimal_cut%d_x%d"%(self.ncuts, c if c >= 0 else -c-1), lhs = None, rhs = cutrhs)
            model.cacheRowExtensions(cut)

            for j in range(len(cutcoefs)):
                if model.isZero(cutcoefs[j]): # maybe here we need isFeasZero
                    continue
                model.addVarToRow(cut, cols[j].getVar(), cutcoefs[j])
            model.flushRowExtensions(cut)
            model.addCut(cut, forcecut=False)
            infeasible = model.addCut(cut, forcecut=False)
            self.ncuts += 1

            if infeasible:
               result = SCIP_RESULT.CUTOFF
            else:
               result = SCIP_RESULT.SEPARATED
            # model.relax()
            # model.writeMIP(filename=f"TEMP/{self.model_name}.{self.exp_id}.{self.row}.lp")
            # model.interruptSolve()
            return {"result": result}
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
