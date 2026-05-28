from cutgeneratingfunctionology.igp import *
from minimalFunctionCache.utils import minimal_function_cache_info
from .cut_generation_problem import *
# the rational conversion QQ is imported from cutgeneratingfunctionology.
#from scipy.optimize import minimize, LinearConstraint, NonlinearConstraint
from pyscipopt import Model, Sepa, SCIP_RESULT
import json
import logging
import os
import time

optimal_cut_logger = logging.getLogger(__name__)
optimal_cut_logger.setLevel(logging.DEBUG)
# Adapted from the example in the docs. https://pymodelopt.readthedocs.io/en/latest/tutorials/separator.html

# Note; pplitepy has numeric issues. Use backend=None when large vaues may be encountered.
class OptimalCut(Sepa):
    def __init__(self, *, write_mip_and_cut=False, file_name_base="OptimalCutData", max_number_of_data_records=None, cgp_kwds=None, paths=None):
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
        self.max_number_of_data_records = max_number_of_data_records
        self.write_mip_and_cut = write_mip_and_cut # writes mip, cut, and metadata used to generate/assess the cut such as row, optimal cut initalization prameters, paths to problem and cut, .
        self.file_name_base = file_name_base
        if paths is None:
            self.mip_and_cut_write_path = ""
            self.metadata_write_path = ""
        else:
            try:
                if self.write_mip_and_cut:
                    self.mip_and_cut_write_path = paths["mip_and_cut_write_path"]
                    self.metadata_write_path = paths["metadata_write_path"]
            except KeyError:
                optimal_cut_logger.debug("paths does not contain keyword \"cgp_mip_data_dir\". Writing to current directory.")
                self.mip_and_cut_write_path = "" # write to current directory
        if cgp_kwds is None:
            self.cgp = cutGenerationProblem()
        else:    
            self.cgp = cutGenerationProblem(**cgp_kwds)

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
        cutcoefs_fun = [0] * len(cols)
        cutrhs_fun = 0

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

                    if .001 <= model.frac(primsol) <= 1 - .001: # use cgp notion of 0/1
                        tryrow = True

            # generate the cut!
            if tryrow:
                optimal_cut_logger.debug(f"c={c}")
                # get the row of B^-1 for this basic integer variable with fractional solution value
                binvrow = model.getLPBInvRow(i)

                # get the tableau row for this basic integer variable with fractional solution value
                binvarow = model.getLPBInvARow(i)

                # get current reduced costs for objective evaluation.
                costs = [model.getColRedCost(j) for j in cols if j not in basisind]

                cgf, cut_score = self.cgp.solve(binvarow, costs, primsol) # produce an optimal cgf
                optimal_cut_logger.debug(f"b={cgf.end_points()}\nv={cgf.values_at_end_points()}")
                cutcoefs, cutrhs = self.getOptimalCutFromRow(cols, rows, binvrow, binvarow, primsol, cgf)

                cut = model.createEmptyRowSepa(self, "optimal_cut%d_x%d"%(self.ncuts,c if c >= 0 else -c-1), lhs = None, rhs = cutrhs)
                model.cacheRowExtensions(cut)

                for j in range(len(cutcoefs)):
                    if model.isZero(cutcoefs[j]): # maybe here we need isFeasZero
                        continue
                    model.addVarToRow(cut, cols[j].getVar(), cutcoefs[j])

                if cut.getNNonz() == 0:
                    assert model.isFeasNegative(cutrhs)
                    return {"result": SCIP_RESULT.CUTOFF}


                # Only take efficacious cuts, except for cuts with one non-zero coefficient (= bound changes)
                # the latter cuts will be handled internally in sepastore.
                if cut.getNNonz() == 1 or model.isCutEfficacious(cut):

                    # flush all changes before adding the cut
                    model.flushRowExtensions(cut)

                    infeasible = model.addCut(cut, forcecut=False)
                    self.ncuts += 1

                    if infeasible:
                       result = SCIP_RESULT.CUTOFF
                    else:
                       result = SCIP_RESULT.SEPARATED
                model.releaseRow(cut)

        return {"result": result}

class GMI(Sepa):

    def __init__(self):
        self.ncuts = 0

    def getGMIFromRow(self, cols, rows, binvrow, binvarow, primsol):
        """ Given the row (binvarow, binvrow) of the tableau, computes gomory cut

        :param primsol:  is the rhs of the tableau row.
        :param cols:     are the variables
        :param rows:     are the slack variables
        :param binvrow:  components of the tableau row associated to the basis inverse
        :param binvarow: components of the tableau row associated to the basis inverse * A

        The GMI is given by
         sum(f_j x_j                  , j in J_I s.t. f_j <= f_0) +
         sum((1-f_j)*f_0/(1 - f_0) x_j, j in J_I s.t. f_j  > f_0) +
         sum(a_j x_j,                 , j in J_C s.t. a_j >=   0) -
         sum(a_j*f_0/(1-f_0) x_j      , j in J_C s.t. a_j  <   0) >= f_0.
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

        cutcoefs_fun = [0] * len(cols)
        cutrhs_fun = 0

        # get scip
        scip = self.model

        # Compute cut fractionality f0 and f0/(1-f0)
        f0 = scip.frac(primsol)
        f = fractional(QQ(primsol))
        ratiof0compl = f0/(1-f0)
        pi_p = gmic(f)
        def psi(x):
            if x < 0:
                slope = pi_p.functions()[-1]._slope
                neg_part = FastLinearFunction(slope, 0)
                return neg_part(x)
            else:
                return pi_p.functions()[0](x)
        # rhs of the cut is the fractional part of the LP solution for the basic variable
        cutrhs = -f0
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
                # warning: because of numerics cutelem < 0 is possible (though the fractional part is, mathematically, always positive)
                # However, when cutelem < 0 it is also very close to 0, enough that isZero(cutelem) is true, so we ignore
                # the coefficient (see below)
                cutelem = scip.frac(rowelem)
                cutelem_fun = fractional(QQ(rowelem))
                cutelem_fun = float(-1*f*pi_p(cutelem_fun))
                if cutelem > f0:
                    # sum((1-f_j)*f_0/(1 - f_0) x_j, j in J_I s.t. f_j  > f_0) +
                    cutelem = -((1.0 - cutelem) * ratiof0compl)
                else:
                    #  sum(f_j x_j                  , j in J_I s.t. f_j <= f_0) +
                    cutelem = -cutelem

            else:
                # Continuous variables
                cutelem_fun = float(-1*f*psi(QQ(rowelem)))
                if rowelem < 0.0:
                    # -sum(a_j*f_0/(1-f_0) x_j      , j in J_C s.t. a_j  <   0) >= f_0.
                    cutelem = rowelem * ratiof0compl
                else:
                    #  sum(a_j x_j,                 , j in J_C s.t. a_j >=   0) -
                    cutelem = -rowelem

            # cut is define when variables are in [0, infty). Translate to general bounds
            if not scip.isZero(cutelem):
                if col.getBasisStatus() == "upper":
                    cutelem = -cutelem
                    cutelem_fun = -cutelem_fun
                    cutrhs += cutelem * col.getUb()
                    cutrhs_fun += cutelem_fun * col.getUb()
                else:
                    cutrhs += cutelem * col.getLb()
                    cutrhs_fun += cutelem_fun * col.getLb()
                # Add coefficient to cut in dense form
                cutcoefs[col.getLPPos()] = cutelem
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
                cutelem = scip.frac(rowelem)
                cutelem_fun = fractional(QQ(rowelem))
                cutelem_fun = float(-1*f*pi_p(cutelem_fun))
                if cutelem > f0:
                    #  sum((1-f_j)*f_0/(1 - f_0) x_j, j in J_I s.t. f_j  > f_0) +
                    cutelem = -((1.0 - cutelem) * ratiof0compl)
                else:
                    #  sum(f_j x_j                  , j in J_I s.t. f_j <= f_0) +
                    cutelem = -cutelem
            else:
                # Continuous variables
                cutelem_fun = float(-1*f*psi(QQ(rowelem)))
                if rowelem < 0.0:
                    # -sum(a_j*f_0/(1-f_0) x_j      , j in J_C s.t. a_j  <   0) >= f_0.
                    cutelem = rowelem * ratiof0compl
                else:
                    #  sum(a_j x_j,                 , j in J_C s.t. a_j >=   0) -
                    cutelem = -rowelem
            if abs(cutelem_fun - cutelem) > .000001:
                optimal_cut_logger.debug(f"cut strengthening \n ratio={cutelem_fun/cutelem}")
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
                    cutcoefs_fun[rowcols[i].getLPPos()] -= cutelem_fun * rowvals[i]

                act = scip.getRowLPActivity(row)
                rhsslack = rrhs - act
                if scip.isFeasZero(rhsslack):
                    assert row.getBasisStatus() == "upper" # cutelem != 0 and row active at upper bound -> slack at lower, row at upper
                    cutrhs -= cutelem * (rrhs - row.getConstant())
                    cutrhs_fun -= cutelem_fun * (rrhs - row.getConstant())
                else:
                    assert scip.isFeasZero(act - rlhs)
                    cutrhs -= cutelem * (rlhs - row.getConstant())
                    cutrhs_fun -= cutelem_fun * (rlhs - row.getConstant())
        optimal_cut_logger.debug(f"f0={f0}")
        optimal_cut_logger.debug(f"cutcoefs={cutcoefs} \n cutcoefs_fun={cutcoefs_fun}")
        optimal_cut_logger.debug(f"cutrhs={cutrhs}")
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

        # For all basic columns (not slacks) belonging to integer variables, try to generate a gomory cut
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
                optimal_cut_logger.debug(f"row={c}")
                # get the row of B^-1 for this basic integer variable with fractional solution value
                binvrow = scip.getLPBInvRow(i)

                # get the tableau row for this basic integer variable with fractional solution value
                binvarow = scip.getLPBInvARow(i)

                # get cut's coefficients
                cutcoefs, cutrhs = self.getGMIFromRow(cols, rows, binvrow, binvarow, primsol)

                # add cut
                cut = scip.createEmptyRowSepa(self, "gmi%d_x%d"%(self.ncuts,c if c >= 0 else -c-1), lhs = None, rhs = cutrhs)
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

                    infeasible = scip.addCut(cut, forcecut=False)
                    self.ncuts += 1

                    if infeasible:
                       result = SCIP_RESULT.CUTOFF
                    else:
                       result = SCIP_RESULT.SEPARATED
                scip.releaseRow(cut)

        return {"result": result}
