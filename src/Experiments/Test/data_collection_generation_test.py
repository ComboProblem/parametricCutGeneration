# from parametricCutGen.scip_data_collection_events import CutGapDataRecording
from parametricCutGen.optimal_cut_generation import OptimalCut, GMI
from parametricCutGen.logging_utils import * # cgp_written_log is here
from parametricCutGen.scip_data_collection_events import record_data
from pyscipopt import Model
from pyscipopt import Model, quicksum, SCIP_PARAMSETTING, exp, log, sqrt, sin
from pyscipopt.scip import Cutsel
from typing import List

import logging

test_logging = logging.getLogger(__name__)
test_logging.setLevel(logging.DEBUG)


model = Model() # from https://github.com/scipopt/PySCIPOpt/blob/master/tests/helpers/utils.py\
model.readProblem(filename="/home/acadia/Downloads/gen-ip016.mps")
#logging.disable()
sepa = OptimalCut(cgp_kwds={"algorithm" : "value_poly_lp", "max_num_of_bkpts": 2}) # From the example of seperator; https://pyscipopt.readthedocs.io/en/latest/tutorials/separator.html]
#sepa = GMI()
#model.readParams("src/Experiments/paramFiles/scip_disable_other_cuts.set") # see link; https://github.com/ComboProblem/parametricCutGeneration/blob/main/src/Experiments/paramFiles/scip_disable_other_cuts.set
model.setSeparating(SCIP_PARAMSETTING.OFF)
model.includeSepa(sepa, "optimal_cut", "optimal_cut test", priority=10000, freq=0)
#model.includeSepa(sepa, "gmi", "gmi test", priority=10000, freq=0)
model.setHeuristics(SCIP_PARAMSETTING.OFF)
model.setPresolve(SCIP_PARAMSETTING.OFF)
#model.disablePropagation()
#cutsel = TestCutsel()
#model.includeCutsel(cutsel, 'test_cut_sel', 'test', 5000000)
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
model.setParam("limits/nodes", 100)
model=record_data(model, 10) # redirecting write pathissue?
#add_cut_at_root = addCutsAtRoot(model)
# disable = disableCuts(model)
# model.includeEventhdlr(disable, "disable_other_cuts", "disable other cuts aside from the targeted cut at depths not at the root")
#model.includeEventhdlr(add_cut_at_root, "add cuts", "cuts_added")
model.optimize()
model.printStatistics()
print(model.data)


from pyscipopt import Model, SCIP_EVENTTYPE, SCIP_RESULT, Eventhdlr, SCIP_PARAMSETTING

class checkCutsAdded(Eventhdlr):
    def __init__(self, model):
        Eventhdlr.__init__(model)
        self.count = 0 
        
    def eventinit(self):
        self.model.catchEvent(SCIP_EVENTTYPE.LPSOLVED, self)

    def eventexec(self, event):  
        if self.model.getCurrentNode().getDepth() > 1: # not at the root; we should have (afik) at most 1 cut added to this node
            self.count += 1
            if self.count % 50 == 0:
                test_logging.debug(f"Node depth: {self.model.getCurrentNode().getDepth()}")
                test_logging.debug(f"Cuts: {self.model.getNCutsApplied()}")
                test_logging.debug(f"sepa_rounds: {self.model.getNSepaRounds()}")
                self.model.writeMIP(f"node_lp_{self.count}.lp")
                test_logging.debug(f"Node type: {self.model.getCurrentNode().getType()}")
                test_logging.debug(f"Node parent: {self.model.getCurrentNode().getParent().getDepth()}")
        if self.model.getCurrentNode().getDepth() == 0:
            test_logging.debug(f"Node depth: {self.model.getCurrentNode().getDepth()}")
            test_logging.debug(f"Node added conss: {self.model.getCurrentNode().getNAddedConss()}")
            test_logging.debug(f"Node type: {self.model.getCurrentNode().getType()}")
        if self.count == 200:
            self.model.interruptSolve()


class addCutsAtRoot(Eventhdlr):
    def __init__(self, model):
        Eventhdlr.__init__(model)
        
    def eventinit(self):
        self.model.catchEvent(SCIP_EVENTTYPE.FIRSTLPSOLVED, self)

    def eventexit(self):
        self.model.dropEvent(SCIP_EVENTTYPE.FIRSTLPSOLVED, self)

    def eventexec(self, event):
        if not model.isLPSolBasic():
            print("boo")
        else:
            if self.model.getCurrentNode().getDepth() == 0:
                print("whoo")

class TestCutsel(Cutsel):

    def cutselselect(self, cuts, forcedcuts, root, maxnselectedcuts):
        """
        Selects the 10 cuts with largest efficacy.
        """

        scip = self.model

        scores = [0] * len(cuts)
        for i in range(len(scores)):
            scores[i] = scip.getCutEfficacy(cuts[i])

        rankings = sorted(range(len(cuts)), key=lambda x: scores[x], reverse=True)

        sorted_cuts = [cuts[rank] for rank in rankings]
        return {'cuts': sorted_cuts, 'nselectedcuts': min(maxnselectedcuts, len(cuts), 1),
            'result': SCIP_RESULT.SUCCESS}

#class disableCuts(Eventhdlr):
#    def __init__(self, model):
#        self.model = model
#    def eventinit(self):
#        self.model.catchEvent(SCIP_EVENTTYPE.NODEFOCUSED, self) # the event is called whenever a node is about to be solved
#
#    def eventexec(self, event):
#        if self.model.getCurrentNode().getDepth() > 1: # if we aren't in the root node
#            self.model.setSeparating(SCIP_PARAMSETTING.OFF) # disable separators
#        else:
#            self.model.setSeparating(SCIP_PARAMSETTING.DEFAULT)
#            self.model.readParams("src/Experiments/paramFiles/scip_disable_other_cuts.set") # see link


#model.printStatistics()


def add_at_most_k_cuts_to_model(model, cgp, k):
    """
    Add the top k scoring cuts to the root based on the inital LP and specificed cut generation problem (as an ).  
    """
    
    result = SCIP_RESULT.DIDNOTRUN
    if not model.isLPSolBasic():
        return {"result": result}

    # get LP data
    cols = model.getLPColsData()
    rows = model.getLPRowsData()

    # exit if LP is trivial
    if len(cols) == 0 or len(rows) == 0:
        return {"result": result}

    result = SCIP_RESULT.DIDNOTFIND
    cuts_to_add = []
    # get basis indices
    basisind = model.getLPBasisInd()

    # For all basic columns (not slacks) belonging to integer variables, try to generate an optimal cut
    for i in range(len(rows)):
        tryrow = False
        c = basisind[i]

        if c >= 0:
            assert c < len(cols)
            var = cols[c].getVar()

            if var.vtype() != "CONTINUOUS":
                primsol = cols[c].getPrimsol()
                assert model.getSolVal(None, var) == primsol

                if cgp._espilon <= model.frac(primsol) <= 1 - cgp._espilon: # use cgp notion of 0/1
                    tryrow = True

        # generate the cut!
        if tryrow:
            # get the row of B^-1 for this basic integer variable with fractional solution value
            binvrow = model.getLPBInvRow(i)

            # get the tableau row for this basic integer variable with fractional solution value
            binvarow = model.getLPBInvARow(i)

            # get current reduced costs for objective evaluation.
            costs = [model.getColRedCost(j) for j in cols if j not in basisind]

            cgf, cut_score = cgp.solve(binvarow, costs, primsol) # produce an optimal cgf
            if len(cuts_to_add) < k:
                cuts_to_add.append((i, cut_score, cgf))
                cuts_to_add = sorted(cuts_to_add, key=lambda elem : elem[1])
            elif cuts_to_add[-1][1] < cut_score:
                cuts_to_add.pop()
                cuts_to_add.append((i, cut_score, cgf))
                cuts_to_add = sorted(cuts_to_add, key=lambda elem : elem[1])
                
            

            
    for cut in cuts_to_add:
        i = cut[0]
        cgf = cut[2]
        c = basisind[i]
        primsol = cols[c].getPrimsol()
        binvrow = model.getLPBInvRow(i)
        # get the tableau row for this basic integer variable with fractional solution value
        binvarow = model.getLPBInvARow(i)
        cutcoefs, cutrhs = getOptimalCutFromRow(cols, rows, binvrow, binvarow, primsol, cgf)
        # model.


    return {"result": result}

def getOptimalCutFromRow(model, cols, rows, binvrow, binvarow, primsol, pi_p):
    """ Given the row (binvarow, binvrow) of the tableau, computes optimized cut.

    :param primsol:  is the rhs of the tableau row
    :param cols:     are the variables
    :param rows:     are the slack variables
    :param binvrow:  components of the tableau row associated to the basis inverse
    :param binvarow: components of the tableau row associated to the basis inverse * A
    :param 1dPWL:    a minimal function, ab element of PiMin<=k with pi_p(primsol)=1

    The intersection cut is given by
     sum(pi_p(a_j) x_j, j in J_I) geq 1
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
    cutcoefs = [0] * len(cols)
    cutrhs = 0

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
            # warning: because of numerics cutelem < 0 is possible (though the fractional part is, mathematically, always positive)
            # However, when cutelem < 0 it is also very close to 0, enough that isZero(cutelem) is true, so we ignore
            # the coefficient (see below)
            cutelem = float(pi_p(fractional(QQ(rowelem)))) #keep types correct: data is always kept as "rational" as interfaces with inexact problems i.e. rational approximations. 
        else:
            # Continuous variables
            # From Matthias, the super additive portion of the function about 0
            def psi(x):
                if x < 0:
                    return pi_p.functions()[-1](x)
                else:
                    return pi_p.functions()[0](x)
            cutelem = float(psi(fractional(QQ(rowelem))))
        # cut is define when variables are in [0, infty). Translate to general bounds
        if not model.isZero(cutelem):
            if col.getBasisStatus() == "upper":
                cutelem = -cutelem
                cutrhs += cutelem * col.getUb()
            else:
                cutrhs += cutelem * col.getLb()
            # Add coefficient to cut in dense form
            cutcoefs[col.getLPPos()] = cutelem

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
            if not model.isInfinity(-row.getLhs()):
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
            cutelem = float(pi_p(fractional(QQ(rowelem))))
        else:
            # Continuous variables
            def psi(x):
                if x < 0:
                    return pi_p.functions()[-1](x)
                else:
                    return pi_p.functions()[0](x)
            cutelem = float(psi(fractional(QQ(rowelem))))

        # cut is define in original variables, so we replace slack by its definition
        if not model.isZero(cutelem):
            # get lhs/rhs
            rlhs = row.getLhs()
            rrhs = row.getRhs()
            assert model.isLE(rlhs, rrhs)
            assert not model.isInfinity(rlhs) or not model.isInfinity(rrhs)

            # If the slack variable is fixed, we can ignore this cut coefficient
            if model.isFeasZero(rrhs - rlhs):
              continue

            # Unflip slack variable and adjust rhs if necessary: row at lower means the slack variable is at its upper bound.
            # Since model adds +1 slacks, this can only happen when constraints have a finite lhs
            if row.getBasisStatus() == "lower":
                assert not model.isInfinity(-rlhs)
                cutelem = -cutelem

            rowcols = row.getCols()
            rowvals = row.getVals()

            assert len(rowcols) == len(rowvals)

            # Eliminate slack variable: rowcols is sorted: [columns in LP, columns not in LP]
            for i in range(row.getNLPNonz()):
                cutcoefs[rowcols[i].getLPPos()] -= cutelem * rowvals[i]

            act = model.getRowLPActivity(row)
            rhsslack = rrhs - act
            if model.isFeasZero(rhsslack):
                assert row.getBasisStatus() == "upper" # cutelem != 0 and row active at upper bound -> slack at lower, row at upper
                cutrhs -= cutelem * (rrhs - row.getConstant())
            else:
                assert model.isFeasZero(act - rlhs)
                cutrhs -= cutelem * (rlhs - row.getConstant())
