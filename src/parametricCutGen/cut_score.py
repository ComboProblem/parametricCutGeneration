"""
Defines cut scores available to the cut generation problem.
"""

from cutgeneratingfunctionology.igp import *
from .generic_solvers import cvxpyCutGenProblemSolverInterface
from .execptions import *
import logging
import time

cut_score_logger  = logging.getLogger(__name__)
cut_score_logger.setLevel(logging.INFO)

def pwl_with_value_parameters_and_bkpts_fixed(bkpt, f_index, log_paramateric_real_field=False, log_pw_functions=False):
    n = len(bkpt)
    if not log_paramateric_real_field:
        parametric_logging_level = logging.getLogger("cutgeneratingfunctionology.igp.parametric").getEffectiveLevel()
        logging.getLogger("cutgeneratingfunctionology.igp.parametric").setLevel(logging.ERROR)
    if not log_pw_functions:
        pw_logging_level = logging.getLogger("cutgeneratingfunctionology.igp.functions").getEffectiveLevel()
        logging.getLogger("cutgeneratingfunctionology.igp.functions").setLevel(logging.ERROR)
    assert(n >= 2)
    assert(f_index >= 1)
    assert(f_index <= n - 1)
    if not isinstance(bkpt, list):
        bkpt = list(bkpt)
    coord_names = ['gamma'+str(i) for i in range(n)]
    K = PolynomialRing(QQ, names=coord_names, order='lex')
    vals = [0] + [K.gens()[i] if i != f_index else 1  for i in range(1,n)]
    return piecewise_function_from_breakpoints_and_values(bkpt + [1],  vals + [0], merge=False)

class abstractCutScore:
    r"""
    Abstract class for cut optimization objective functions aka cut scores.
    Named after heuristics used to evaluate a cuts effectiveness.
    """
    def __init__(self):
        pass
        
    def cut_score(cut, mip_obj):
        r"""
        A(n) (assumed to be) smooth function from {(\pi(bar a_ij))_{j\in N} : pi in PiMin} to RR.

        Suppose that R is a ring such that either QQ subseteq R subseteq RR or  R is a ring that can
        coercied to a ring R' with QQ subseteq R' subseteq RR.

        cut_score should use sagemath types to ensure generating a separating cut.

        input: cut, mip_obj
        cut: A list of length n-m of elements of R representing a proposed cut to a given MIP.
        mip_obj: A list of elements of R representing the MIPs objective function.

        output: an element of R.
        """
        raise NotImplementedError

    def cut_score_grad(cut, mip_obj):
        r"""
        The gradient of the cut score function.

        input: cut, mip_obj
        cut: A list of length n-m of elements of R representing a proposed cut to a given MIP.
        mip_obj: A list of elements of R representing the MIPs objective function.

        output: A vector of length n-m of elements of R.
        """
        raise NotImplementedError

    def cut_score_hess(cut, mip_obj):
        r"""
        The hessian of the cut score function.

        input: cut, mip_obj
        cut: A list of length n-m of elements of R representing a proposed cut to a given MIP.
        mip_obj: A list of elements of R representing the MIPs objective function.

        output: An n-m by n-m matrix of elements of R.
        """
        raise NotImplementedError

    def is_linear():
        r"""
        Returns True if the cut score (when defined on the parameterized cut space) is linear.
        """
        raise NotImplementedError
    
    def wrap_cut_score_to_solver_linear_objective(solver, **kwrds):
        r"""
        Returns a valid input linear function for solver to use in an LP problem.
        """
        raise NotImplementedError

class cutScore:
    """
    cutScore is objective function used in the cutOptimzationProblem.

    cutScore is a(n) (assumed to be) smooth function from cutSpace to RR where cutSpace = {(\pi_p(bar a_ij))_{j\in N} : pi in PiMin<=k}.

    cutScore's domain is the cutSpace written in terms of the parameterization of PiMin<=k.

    cutScore comes with optional methods to provide first and second order information to solvers.
    
    TESTS::
    >>> from parametricCutGen.cut_score import *
    >>> cs = cutScore("parallelism")
    >>> cs
    cut score Parallelism
    >>> class TestCutScore(AbstractCutScore):
            pass
    >>> cutScore(TestCutScore)
    cut score TestCutScore
    """
    def __init__(self, *, cut_score=None, objective_sense="maximize"):
        r"""
        Initialize the paramatrized cut scoring function.

        Data used with cutScore is managed by the methods that call cutScore.
        """
        if cut_score is None:
            self._cut_score = SCIP_STANDARD
        elif cut_score is "parallelism":
            self._cut_score = Parallelism
        elif cut_score is "steepest_direction":
            self._cut_score = SteepestDirection
        elif cut_score is "2*steepest_direction":
            self._cut_score = SteepestDirection2
        elif issubclass(cut_score, abstractCutScore):
            self._cut_score = cut_score
        else:
            raise TypeError("Use a predefined cut scoring method or use custom instance of abstractCutScore.")
        self._MIP_objective = None
        self._MIP_row = None
        self._sage_to_solver_type = None
        self._timer = None
        self._feasible_point = None
        self._rel_tol = 10**-9
        self._prev_result = None
        self._objective_sense = objective_sense
        if self._objective_sense not in ["maximize", "minimize"]:
            raise ValueError("Objective senese should be either maximize or minimize.")

    def __repr__(self):
        return f"cut score {self._cut_score.__name__}"

    def __call__(self, parameters):
        r"""
        Evaluate the cutScore.

        parameters is a list like object of real numbers with even length of at most 2k.
        parameters = (bkpt, val) represents a parameterized element of Pimin<=k by the breakpoint and value paramaterization.
        """
        # It is necessary for frac_f to be converted to exact rational from the MIP.
        if self._MIP_objective is None:
            raise UnsetData("Set MIP_objective before use of CutScore.")
        if self._MIP_row is None:
            raise UnsetData("Set MIP_row before use of CutScore.")
        if self._sage_to_solver_type is None:
            raise UnsetData("Set sage_to_solver_type before use of CutScore.")
        # To generate a separator, we have to ensure, the function pi_parameters
        # is minimal.
        # parameters is given from some non linear solver and might not meet the conditions
        # of minimality.
        # validate_point will either give a point p = b,v which
        # we believe to up to rounding and L.C. that pi_p is a minimal function in the current cell
        # and satisfies the conditions of the model or will raise an error.
        b, v = self.validate_point(parameters)
        self.set_feasible_point(b+v)
        pi = piecewise_function_from_breakpoints_and_values(b + [1], v + [0])
        row_data = self.get_MIP_row()
        sage_cut = [pi(fractional(QQ(bar_a_ij))) for bar_a_ij in row_data]
        sage_mip_obj =  [QQ(bar_cj) for bar_cj in self._MIP_objective]
        sage_result = self._cut_score.cut_score(sage_cut, sage_mip_obj)
        if self._objective_sense == "minimize":
            sage_result = -1*sage_result
        self._sage_cut = sage_cut
        self._sage_mip_obj = sage_mip_obj
        if self.get_prev_result() is not None and sage_result != 0:
            if abs(sage_result - self.get_prev_result())/sage_result < self._rel_tol:
                cut_score_logger.debug(f"cutScore.__call__: Relative distance between successive solutions is less than {self._rel_tol}. Stopping non-linear solver.")
                self.set_prev_result(sage_result)
                raise SolverRelTolReached(f"cutScore.__call__: Relative distance between successive solutions is less than {self._rel_tol}. Stopping non-linear solver.")
            else:
                self.set_prev_result(sage_result)
        else:
            self.set_prev_result(sage_result)
        if self._timer is not None:
            if self._timer.solver_time_out():
                raise SolverTimeOut
        return self.get_sage_to_solver_type()(sage_result)

    @staticmethod
    def grad(self):
        """
        Gradient of cutScore. This method should only be called after a __call__ has been made to cutScore.
        """
        return self._cut_score.cut_score_grad(self._sage_cut, self._sage_mip_obj)

    @staticmethod
    def hess(self):
        """
        Hessian of cutScore. This method should only be called after a __call__ has been made to cutScore.
        """
        return self._cut_score.cut_score_hess(self._sage_cut, self._sage_mip_obj)
        

### Get and set methods for communicating data between solvers.

    def get_current_cell(self):
        return self._cell

    def get_espilon(self):
        return self._espilon

    def get_f_index(self):
        return self._f_index

    def get_f_trust(self):
        return self._f_trust

    def get_feasible_point(self):
        return self._feasible_point

    def get_lipschitz_constant(self):
        return self._M

    def get_objective_sense(self):
        return self._objective_sense

    def get_MIP_obj(self):
        return self._MIP_objective

    def get_MIP_row(self):
        """
        For fixed row i of the corner polyhedron;  bar_i - (bar a_i)^T x_N
        """
        return self._MIP_row

    def get_prev_result(self):
        return self._prev_result

    def get_sage_to_solver_type(self):
        """
        Defined in the cutGenerationSolver. This method is only intended to be set and unset by solving
        routines in cutGenerationSolver.
        """
        return self._sage_to_solver_type

    def set_current_cell(self, cell):
        self._cell = cell

    def set_espilon(self, espilon):
        self._espilon = espilon

    def set_f_index(self, f_index):
        self._f_index = f_index

    def set_f_trust(self, f_trust):
        self._f_trust = f_trust

    def set_feasible_point(self, point):
        self._feasible_point = point

    def set_lipschitz_constant(self, M):
        self._M = M

    def set_objective_sense(self, objective_sense):
        self._objective_sense = objective_sense

    def set_MIP_row(self, new_row):
        self._MIP_row = new_row

    def set_MIP_obj(self, new_objective):
        """
        Use reduced costs of the basis relaxation here.

        This should be a vector of length n-m for the problem.
        """
        self._MIP_objective = new_objective

    def set_prev_result(self, prev_result):
        self._prev_result = prev_result

    def set_rel_tol(self, rel_tol):
        self._rel_tol = rel_tol

    def set_sage_to_solver_type(self, new_conversion):
        self._sage_to_solver_type = new_conversion

    def set_timer(self, timer):
        """
        Timer is a class which keeps track of the solvers time spent solving the cgf problem.
        """
        self._timer = timer

    def objective_sense(self):
         self._objective_sense

    def validate_point(self, point, poly_check=False):
        """
        Take parameters from the non linear solver and returns point that satisfies the minimal function model.
        """
        # This enforces the "manifoldness" of PWL.
        epsilon = self._espilon
        M = self._M
        f_index = self._f_index
        f_trust = fractional(QQ(self._f_trust))
        cell = self._cell
        sage_point = [QQ(x) for x in point]
        n = int(len(point)/2)
        b, v  = [b for b in  sage_point[:n]], [v for v in sage_point[n:]]
        # model assumptions
        # we log all errors
        # breakpoints should be in [0,1)
        # values should be in [0,1]
        for i in range(n):
            if (b[i] < 0 or b[i] + epsilon >= 1) and i != f_index:
                raise ModelViolation(f"validate_point: breakpoint lambda_{i} < 0 or lambda_{i}>= 1; lambda_{i}=={b[i]}")
            if (v[i] < 0 or v[i] + epsilon > 1) and i != f_index:
                raise ModelViolation(f"validate_point: breakpoint gamma_{i} < 0 or gamma_{i}>= 1; gamma_{i}=={v[i]}")
        # pi(0) = 0, pi(f) = 1
        if abs(b[0]-0) <= epsilon:
            b[0] = 0
        else:
            cut_score_logger.debug(f"validate_point: breakpoint lambda_0 >0, point: {point}, cell: {cell}")
            raise ModelViolation("breakpoint lambda_0 >0")
        if abs(v[0]-0) <= epsilon:
            v[0] = 0
        else:
            cut_score_logger.debug(f"validate_point: breakpoint gamma_0 >0, point: {point}, cell: {cell}")
            raise ModelViolation("value gamma_0 >0")
        # bkpt[f_index] == f
        if abs(b[f_index] - f_trust) <= epsilon:
            b[f_index] = f_trust
        else:
            cut_score_logger.debug(f"validate_point: breakpoint lambda_{f_index} != {f_trust}: point: {point}, cell: {cell}")
            raise ModelViolation(f"breakpoint lambda_{f_index} != {f_trust}")
        # pi_p(f) = 1
        if abs(v[f_index] - 1) <= epsilon:
            v[f_index] = 1
        else:
            cut_score_logger.debug(f"validate_point: breakpoint lambda_{f_index} !=1: point: {point}, cell: {cell}")
            raise ModelViolation(f"value gamma_{f_index} != 1")
        # lipschitz constant and continuity.
        for i in range(n-1):
            if 0 < b[i+1]-b[i]<= epsilon:
                if abs(v[i+1]-v[i]) >= epsilon*M:
                    # potential discontunity
                    # not in (epsilon_i, M) - charts.
                    cut_score_logger.debug(f"validate_point: Solution does not have lipschitz constant {M}: point: {point}, cell: {cell}")
                    raise ModelViolation(f"Solution does not have lipschitz constant {M}")
                # lambda_i+1 = lambda_i
                b[i+1] = b[i]
                # continuity, gamma_i =gamma_i+1
                if abs(v[i+1]-v[i]) > epsilon:
                    raise ModelViolation
                v[i+1] = v[i]
                # when epsilon <= v[i+1]-v[i] < epsilon*M
                # the solution exists in the intersection of the epsilon,M
                # can assume the values are correct/as intended
        # the last breakpoint should be distinct from 1 to enforce a breakpoint sequence.
        if 1-b[n-1] <= epsilon:
            cut_score_logger.debug(f"validate_point: breakpoint lambda_{n-1} >= 1: point: {point}, cell: {cell}")
            raise ModelViolation(f"breakpoint lambda_{n-1} >= 1")
        # b,v are rounded values.
        # ensure constraints hold
        # For a polynomial constraint, poly, assume poly(b,v) = 0. Then poly(point) = poly((b,v) + O(epsilon))
        # implies  poly(point) = poly((b,v) + O(epsilon)) = poly(b,v) + total_deriv(poly)( point)+O(epspsilon)) + O(epsilon^2)
        # We treat episolon^2 -> 0. Hence poly(point) = total_deriv(poly)(O(epsilon)) <= total_deriv(poly)(epsilon)
        # the above does not hold, then poly(b,v) != 0.
        # note the polynomial parents should have variable names lambda_0,...,lambda_n-1, gamma_0,...,gamma_n-1 in order.
        # hence we can lazily evaluate polys from the cell.
        # this does not seem to make a difference so it should be safe to ignore.
        if poly_check:
            # Set poly_check to true to prove within the chart used, the constraints hold.
            for poly in cell.lt_poly():
                poly_val = poly(b+v)
                if poly_val < -1*epsilon:
                    pass
                else:
                    if poly_val < 0:
                        # see if poly_val == 0 or not.
                        # assume poly(b+v) == 0.
                        if abs(poly(sage_point)) <= sum(abs(grad(sage_point)) for grad in poly.gradient())*epsilon:
                            cut_score_logger.debug(f"validate_point: {poly} evaluated at {b+v} == 0 when {poly} evaluated at {b+v} should be < 0: point: {point}, cell: {cell}")
                            raise ModelViolation(f"{poly} evaluated at {b+v} == 0 when {poly} evaluated at {b+v} should be < 0")
                        else:
                            pass
                    else:
                        cut_score_logger.debug(f"validate_point: {poly} evaluated at {b+v} >= 0 when {poly} evaluated at {b+v} should be < 0: point: {point}, cell: {cell}")
                        raise ModelViolation(f"{poly} evaluated at {b+v} >= 0 when {poly} evaluated at {b+v} should be < 0")
            for poly in cell.le_poly():
                if poly(b+v) > 0:
                    cut_score_logger.debug(f"validate_point: {poly} evaluated at {b+v} >0 when {poly} evaluated at {b+v} should be <=  0: point: {point}, cell: {cell}")
                    raise ModelViolation(f"{poly} evaluated at {b+v} >0 when {poly} evaluated at {b+v} should be <= 0")

            for poly in cell.eq_poly():
                if abs(poly(b+v)) >= epsilon or abs(poly(sage_point)) > sum(abs(grad(sage_point)) for grad in poly.gradient())*epsilon:
                    cut_score_logger.debug(f"validate_point:{poly} evaluated at {b+v} != 0 when {poly} evaluated at {b+v} should be == 0: point: {point}, cell: {cell}")
                    raise ModelViolation(f"{poly} evaluated at {b+v} != 0 when {poly} evaluated at {b+v} should be == 0")
        return b,v


class Parallelism(abstractCutScore):
    """
    Normalized cut parallelism score.
    """
    def cut_score(self, cut, mip_obj):
        obj_norm = vector(mip_obj).norm()
        cut_norm = vector(cut).norm()
        dot_product = vector(mip_obj).row()*vector(cut).column()
        return (dot_product[0]/(obj_norm*cut_norm))[0]
    
    def is_linear(self):
        r"""
        Returns True if the cut score (when defined on the parameterized cut space) is linear.
        """
        return False


class SteepestDirection(abstractCutScore):
    """
    Steepest direction score.
    """
    def cut_score(cut, mip_obj):
        dot_product = vector(mip_obj).row()*vector(cut).column()
        return dot_product[0][0]

    def is_linear():
        r"""
        Returns True if the cut score (when defined on the parameterized cut space) is linear.
        """
        return True
    
    def wrap_cut_score_to_solver_linear_objective(solver, **kwds):
        r"""
        Returns a valid input linear function for solver to use in an LP problem.
        """
        if issubclass(solver, cvxpyCutGenProblemSolverInterface):
            from cvxpy import Minimize, Maximize
            x = kwds['x']
            bkpt = kwds['bkpt']
            f_index = kwds['f_index']
            mip_obj = kwds['mip_obj']
            pi = pwl_with_value_parameters_and_bkpts_fixed(bkpt, f_index)
            cut_score_in_value_params = sum(pi(fractional(QQ(c))) for c in mip_obj)
            coord_names = ['gamma'+str(i) for i in range(len(bkpt))]
            param_obj = np.array([cut_score_in_value_params.coefficient(cut_score_in_value_params.parent().gens_dict()[name]) for name in coord_names])
            cut_score_logger.debug(f"Objective inputs... {param_obj, x}")
            if kwds['objective_sense'] == "maximize":
                return Maximize(param_obj @ x)
            else:
                return Minimize(param_obj @ x)

class SteepestDirection2(abstractCutScore):
    """
    Steepest direction score.
    """
    def cut_score(cut, mip_obj):
        dot_product = 2*vector(mip_obj).row()*vector(cut).column()
        return dot_product[0][0]

    def is_linear():
        r"""
        Returns True if the cut score (when defined on the parameterized cut space) is linear.
        """
        return True
    
    def wrap_cut_score_to_solver_linear_objective(solver, **kwds):
        r"""
        Returns a valid input linear function for solver to use in an LP problem.
        """
        if issubclass(solver, cvxpyCutGenProblemSolverInterface):
            from cvxpy import Minimize, Maximize
            x = kwds['x']
            bkpt = kwds['bkpt']
            f_index = kwds['f_index']
            mip_obj = kwds['mip_obj']
            pi = pwl_with_value_parameters_and_bkpts_fixed(bkpt, f_index)
            cut_score_in_value_params = sum(pi(fractional(QQ(c))) for c in mip_obj)
            coord_names = ['gamma'+str(i) for i in range(len(bkpt))]
            param_obj = np.array([cut_score_in_value_params.coefficient(cut_score_in_value_params.parent().gens_dict()[name]) for name in coord_names])
            if kwds['objective_sense'] == "maximize":
                return Maximize(2*param_obj @ x)
            else:
                return Minimize(2*param_obj @ x)


class SCIP_STANDARD(abstractCutScore):
    """
    
    """
    def cut_score(cut, mip_obj):
        raise NotImplementedError
