from cutgeneratingfunctionology.igp import *
from minimalFunctionCache.utils import minimal_function_cache_info
# the rational conversion QQ is imported from cutgeneratingfunctionology.
from scipy.optimize import minimize, LinearConstraint, NonlinearConstraint
from pyscipopt import Model, Sepa, SCIP_RESULT
import logging
import time


class OptimalCut(Sepa):
    def __init__(self, algorithm=None, backend=None, cut_score=None,  epsilon=10**-7, M=10**7, max_cgf_solver_time=None, max_num_of_bkpts=None, multithread=False,
        paramaterized_solver=None, prove_seperator=False, rel_tol=10**-6, show_proof=False):
        # the signature here needs to be the same as the cgp.
        # TODO: Change this to a keywords type intitalization.
        self.ncuts = 0
        self.cgp = cutGenerationProblem(algorithm=algorithm, backend=backend, cut_score=cut_score,  epsilon=epsilon, M=M, max_cgf_solver_time=max_cgf_solver_time, max_num_of_bkpts=max_num_of_bkpts, multithread=multithread,
        paramaterized_solver=paramaterized_solver, prove_seperator=prove_seperator, rel_tol=rel_tol, show_proof=show_proof)
# Adapted from the example in 
    def getOptimalCutFromRow(self, cols, rows, binvrow, binvarow, primsol, pi_p):
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
        cutcoefs = [0] * len(cols)
        cutrhs = 0

        # get scip
        scip = self.model

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
                cutelem = float(pi_p(fractional(QQ(rowelem)))) #keep types correct
            else:
                # Continuous variables
                # From Matthias, the super additive portion of the function about 0
                def psi(x):
                    if x < 0:
                        return pi_p.functions()[-1](x)
                    else:
                        return pi_p.functions()[0](x)
                cutelem = float(psi(fractional(QQ)(rowelem)))
            # cut is define when variables are in [0, infty). Translate to general bounds
            if not scip.isZero(cutelem):
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
                cutelem = float(pi_p(fractional(QQ(rowelem))))
            else:
                # Continuous variables
                def psi(x):
                    if x < 0:
                        return pi_p.functions()[-1](x)
                    else:
                        return pi_p.functions()[0](x)
                cutelem = float(psi(fractional(QQ)(rowelem)))

            # cut is define in original variables, so we replace slack by its definition
            if not scip.isZero(cutelem):
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
                    cutelem = -cutelem

                rowcols = row.getCols()
                rowvals = row.getVals()

                assert len(rowcols) == len(rowvals)

                # Eliminate slack variable: rowcols is sorted: [columns in LP, columns not in LP]
                for i in range(row.getNLPNonz()):
                    cutcoefs[rowcols[i].getLPPos()] -= cutelem * rowvals[i]

                act = scip.getRowLPActivity(row)
                rhsslack = rrhs - act
                if scip.isFeasZero(rhsslack):
                    assert row.getBasisStatus() == "upper" # cutelem != 0 and row active at upper bound -> slack at lower, row at upper
                    cutrhs -= cutelem * (rrhs - row.getConstant())
                else:
                    assert scip.isFeasZero(act - rlhs)
                    cutrhs -= cutelem * (rlhs - row.getConstant())

        return cutcoefs, cutrhs

    def sepaexeclp(self):
        result = SCIP_RESULT.DIDNOTRUN
        scip = self.model

        if not scip.isLPSolBasic():
            return {"result": result}

        # get LP data
        cols = scip.getLPColsData()
        rows = scip.getLPRowsData()

        # exit if LP is trivial
        if len(cols) == 0 or len(rows) == 0:
            return {"result": result}

        result = SCIP_RESULT.DIDNOTFIND

        # get basis indices
        basisind = scip.getLPBasisInd()

        # For all basic columns (not slacks) belonging to integer variables, try to generate an optimal cut
        for i in range(len(rows)):
            tryrow = False
            c = basisind[i]

            if c >= 0:
                assert c < len(cols)
                var = cols[c].getVar()

                if var.vtype() != "CONTINUOUS":
                    primsol = cols[c].getPrimsol()
                    assert scip.getSolVal(None, var) == primsol

                    if 0.005 <= scip.frac(primsol) <= 1 - 0.005:
                        tryrow = True

            # generate the cut!
            if tryrow:
                # get the row of B^-1 for this basic integer variable with fractional solution value
                binvrow = scip.getLPBInvRow(i)

                # get the tableau row for this basic integer variable with fractional solution value
                binvarow = scip.getLPBInvARow(i)

                # get current reduced costs for objective evaluation.
                costs = [scip.getColRedCost(j) for j in cols if j not in basisind]

                cgf = self.cgp.solve(binvarow, costs, primsol) # produce an optimal cgf

                cutcoefs, cutrhs = self.getOptimalCutFromRow(cols, rows, binvrow, binvarow, primsol, cgf)

                # add cut
                cut = scip.createEmptyRowSepa(self, "optimalcut%d_x%d"%(self.ncuts,c if c >= 0 else -c-1), lhs = None, rhs = cutrhs)
                scip.cacheRowExtensions(cut)

                for j in range(len(cutcoefs)):
                    if scip.isZero(cutcoefs[j]): # maybe here we need isFeasZero
                        continue
                    scip.addVarToRow(cut, cols[j].getVar(), cutcoefs[j])

                if cut.getNNonz() == 0:
                    assert scip.isFeasNegative(cutrhs)
                    return {"result": SCIP_RESULT.CUTOFF}


                # Only take efficacious cuts, except for cuts with one non-zero coefficient (= bound changes)
                # the latter cuts will be handled internally in sepastore.
                if cut.getNNonz() == 1 or scip.isCutEfficacious(cut):

                    # flush all changes before adding the cut
                    scip.flushRowExtensions(cut)

                    infeasible = scip.addCut(cut, forcecut=True)
                    self.ncuts += 1

                    if infeasible:
                       result = SCIP_RESULT.CUTOFF
                    else:
                       result = SCIP_RESULT.SEPARATED
                scip.releaseRow(cut)

        return {"result": result}
