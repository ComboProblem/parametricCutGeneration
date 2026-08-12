from cutgeneratingfunctionology.igp import *
from minimalFunctionCache.utils import minimal_function_cache_info
from .cut_generation_problem import *
from pyscipopt import Model, Sepa, SCIP_RESULT
import logging

# QQ is imported implicitly from cutgeneratingfunctionology.igp

optimal_cut_logger = logging.getLogger(__name__)
optimal_cut_logger.setLevel(logging.ERROR)

# Adapted from the example in the docs. hhttps://pyscipopt.readthedocs.io/en/latest/tutorials/separator.html

class OptimalCut(Sepa):
    """
    Impements an interface to pyscipopt for a ``cutGenerationProblem``.

    :Parameters:
    \'cgp_kwds\' - Keywords for the :class:``cutGenerationProblem``.
    \'max_cuts\' - ``None`` or positive integer. Maximum number of cuts to gengerate.
    \'cgp_timing\' - ``None``, \'scip\' or real number. Bounds the number of cuts by cgp_timing 
    
    """
    def __init__(self, *, cgp_kwds={}, cgp_timing='scip', max_cuts=None, write_cgf_data=False):
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
        self.write_cgf_data = write_cgf_data
        self.max_cuts = max_cuts
        if cgp_kwds is None:
            self.cgp = cutGenerationProblem()
        else:    
            self.cgp = cutGenerationProblem(**cgp_kwds)
        # For help keeping within scip limits/time. There is no timing enforced in while solving a cgp at the moment. 
        self.cgp_timing = cgp_timing
        if self.cgp_timing:
            self.cgp_solve_time = 0
            self.ncgp_solves = 0
            self.mean_cgp_solve_time = 0

    def getOptimalCutFromRow(self, cols, rows, binvrow, binvarow, primsol, pi_p):
        """ Given the row (binvarow, binvrow) of the tableau, computes optimized cut generated
        from the cut generating problem as specified by ``cgp_kwds``.

        :param primsol:  is the rhs of the tableau row
        :param cols:     are the variables
        :param rows:     are the slack variables
        :param binvrow:  components of the tableau row associated to the basis inverse
        :param binvarow: components of the tableau row associated to the basis inverse * A
        :param 1dPWL:    a minimal function, ab element of PiMin<=k with pi_p(primsol)=1
    
        The cut generated is an intersection cut.
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
        # For all basic columns (not slacks) belonging to integer variables, decide if the row should be tried. 
        to_try = []
        for i in range(len(rows)):
            c = basisind[i]
            if c >= 0:
                assert c < len(cols)
                var = cols[c].getVar()

                if var.vtype() != "CONTINUOUS":
                    primsol = cols[c].getPrimsol()
                    assert model.getSolVal(None, var) == primsol
                    # while cut generating problems can address small value, we don't konw the effects on the cuts numerical stablility.
                    if .01 <= model.frac(primsol) <= 1 - .01:
                        to_try.append(i)
        # pick the order to process the rows we will try to generate cuts
        process_order = sorted(to_try, key = lambda i: abs(model.frac(cols[basisind[i]].getPrimsol())-1/2))
        # work on generating an optimal cut. 
        for i in process_order:
            if self.cgp_timing is not None:
                # this is not 
                if self.cgp_timing == "scip":
                    max_time = model.getParam("limits/time")
                else:
                    max_time = self.cgp_timing
                if max_time - model.getTotalTime()- self.mean_cgp_solve_time < 0:
                    break
            c = basisind[i]
            primsol = cols[c].getPrimsol()
            optimal_cut_logger.debug(f"c={c}")
            # get the row of B^-1 for this basic integer variable with fractional solution value
            binvrow = model.getLPBInvRow(i)

            # get the tableau row for this basic integer variable with fractional solution value
            binvarow = model.getLPBInvARow(i)

            # get current reduced costs for objective evaluation.
            costs = [model.getColRedCost(j) for j in cols if j not in basisind]

            # start clock
            if self.cgp_timing is not None:
                cgp_start_time = model.getTotalTime()

            cgf, cut_score, prob_dim = self.cgp.solve(binvrow, binvarow, primsol, costs, cols, rows, model) # produce an optimal cgf

            # keep running track of average solve time.
            if self.cgp_timing is not None:
                self.cgp_solve_time += model.getTotalTime()-cgp_start_time
                self.ncgp_solves += 1
                self.mean_cgp_solve_time = self.cgp_solve_time/self.ncgp_solves

        
            if self.write_cgf_data:
                if not hasattr(model, "data") or model.data==None:
                    model.data = {}
                    model.data["cgf_log"] = []
                b = cgf.end_points()
                v = cgf.values_at_end_points()
                model.data["cgf_log"].append(((b,v), cut_score, prob_dim, self.ncuts, c if c >= 0 else -c-1)) # can be used to map cgfs to cuts applied. 

            optimal_cut_logger.debug(f"b={cgf.end_points()}\nv={cgf.values_at_end_points()}")

            cutcoefs, cutrhs = self.getOptimalCutFromRow(cols, rows, binvrow, binvarow, primsol, cgf)

            cut = model.createEmptyRowSepa(self, "optimal_cut%d_x%d"%(self.ncuts, c if c >= 0 else -c-1), lhs = None, rhs = cutrhs)
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
