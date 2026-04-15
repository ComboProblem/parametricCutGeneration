from cutgeneratingfunctionology.igp import *
# the rational conversion QQ is imported from cut generating functionology.
from scipy.optimize import minimize, LinearConstraint, NonlinearConstraint
from pyscipopt import Model, Sepa, SCIP_RESULT
import logging
import time

cut_gen_logger = logging.getLogger(__name__)
cut_gen_logger.setLevel(logging.WARNING)

minimal_function_cashe_logging = True

max_bkpts = 4


def find_f_index(min_pwl):
    r"""
    Assume a model function, pi_p, with a fininte nubmer of breakpoints is given.
    Finds the index i such that pi_p(lambda_i) = 1.

    INPUT:
    - igp piecewise linear function

    OUTPUT:
    - integer
    """
    return min_pwl.end_points().index(find_f(min_pwl))


def sparse_enough_breakpoints(bkpt_old, epsilon):
    r"""
    Considers the space PWL(*<=n) and finds a point breakpoint seqeunce bkpt
    such that  ||bkpt_old - bkpt||_infty < epsilon and that  bkpt[i+1]-bkpt[i] > epsilon xor bkpt[i] = bkpt[i+1].
    This imples that d((bkpt_old,f), (bktp,f))<espilon and bkpt is "sparse enough".

    INPUT:
    - breakpoint sequence

    OUTPUT:
    - breakpoint sequence such that bkpt[i+1]-bkpt[i] > epsilon xor bkpt[i] = bkpt[i+1].
    """
    bkpt = list(tuple(bkpt_old)) # cheap way of deep copying lists
    for i in range(len(bkpt)-1):
        if abs(bkpt[i]) < epsilon:
            bkpt[i] = 0
        # 1 equiv 0 mod 1
        elif abs(bkpt[i] - 1) < epsilon:
            bkpt[i] = 0
        elif abs(bkpt[i] - bkpt[i+1]) < epsilon:
            bkpt[i+1] = bkpt[i]
    if abs(bkpt[len(bkpt)-1] -1) < epsilon or abs(bkpt[len(bkpt)-1])< epsilon:
        bkpt[len(bkpt)-1] = 0
    return bkpt


def log_problem_result(bkpt, val, binvarow, binvc, f):
    cut_gen_logger.info(f"Cut generation problem solved.")
    cut_gen_logger.info(f"breakpoitns: {bkpt}")
    cut_gen_logger.info(f"values: {val}")
    cut_gen_logger.info(f"row: {binvarow}")
    cut_gen_logger.info(f"objective:{binvc}")
    cut_gen_logger.info(f"f: {f}")


class UnsetData(Exception):
    pass


class SolverError(Exception):
    pass


class SolverHalt(Exception):
    pass


class SolverTimeOut(Exception):
    pass


class cgfTimer:
    def __init__(self, max_time):
        if max_time is None:
            max_time = 2**63-1
        self._max_time = max_time
        self._start_time  =  time.process_time()

    def solver_time_out(self):
        if time.process_time() - self._start_time >= self._max_time:
            cut_gen_logger.warning(f"Solver timed out.")
            return True
        return False

class abstractCutScore:
    r"""
    Abstract class for cut optimization objective functions aka cut scores.
    Named after huersticis used to evaluate a cuts effectiveness.
    """
    @classmethod
    def __init__(cls, **kwrds):
        pass

    def cut_score(cls, cut, mip_obj):
        r"""
        A(n) (assumed to be) smooth function from cutSpace to RR where cutSpace = {(\pi(bar a_ij))_{j\in N} : pi in PiMin}.

        Suppose that R is a ring such that either QQ subseteq R subseteq RR or  R is a ring that can
        coercied to a ring R' wiht QQ subseteq R' subseteq RR.

        cut_score should use sagemath types to ensure generating a separating cut.

        input: cut, mip_obj
        cut: A list of length n-m of elements of R representing a proposed cut to a given MIP.
        mip_obj: A list of elements of R representing the MIPs objective function.

        output: an element of R.

        EXAMPLES:
        """
        raise NotImplementedError

    def cut_score_grad(cls, cut, mip_obj):
        r"""
        The gradient of the cut score function.

        input: cut, mip_obj
        cut: A list of length n-m of elements of R representing a proposed cut to a given MIP.
        mip_obj: A list of elements of R representing the MIPs objective function.

        output: A vector of length n-m of elements of R.
        """
        raise NotImplementedError

    def cut_score_hess(cls, cut, mip_obj):
        r"""
        The hessian of the cut score function.

        input: cut, mip_obj
        cut: A list of length n-m of elements of R representing a proposed cut to a given MIP.
        mip_obj: A list of elements of R representing the MIPs objective function.

        output: An n-m by n-m matrix of elements of R.
        """
        raise NotImplementedError


class cutScore:
    """
    cutScore is objective function used in the cutOptimzationProblem.

    cutScore is a(n) (assumed to be) smooth function from cutSpace to RR where cutSpace = {(\pi_p(bar a_ij))_{j\in N} : pi in PiMin<=k}.

    cutScore's domain is the cutSpace written in terms of the parameterization of PiMin<=k.

    cutScore comes with optional methods to provide first and second order information to solvers.
    """
    @staticmethod
    def __classcall__(cls, name=None, **kwrds):
        r"""
        Input normalization of cutScore class.
        """
        if name == "parallelism" or name is None:
            return super().__classcall__(cls, cut_score=Parallelism)
        if name == "steepest_direction" or name is None:
            return super().__classcall__(cls, cut_score=SteepestDirection)
        if issubclass(name, abstractCutScore):
            return super().__classcall__(cls, cut_score=name, **kwrds)
        else:
            raise TypeError("Use a predefined cut scoring method or use custom instance of abstractCutScore.")

    def __init__(self, cutscore, **kwrds):
        r"""
        Initialize the paramatrized cut scoring function.

        Data used with cutScore is managed by the methods that call cutScore.
        """
        self._cut_score = cutscore(**kwrds)
        super().__init__()
        self._MIP_objective = None
        self._MIP_row = None
        self._sage_to_solver_type = None
        self._timer = None
        self._feasible_point = None
        self._rel_tol = 10**-6
        self._prev_result = None
        if "obj_type" in kwrds.keys():
            self._cut_obj_type = kwrds.keys()["obj_type"]
        else:
            self._cut_obj_type = "max"

    def __call__(self, parameters):
        r"""
        Evaluate the cutScore.

        parameters is a list like object of real numbers with even length of at most 2k.
        parameters = (bkpt, val) represents a parameterized element of Pimin<=k by the breakpoint and value paramaterization.

        EXAMPLES::

        """
        # It is necessary for frac_f to be converted to exact rational from the MIP.
        if self._MIP_objective is None:
            raise UnsetData("Set MIP_objective before use of CutScore.")
        if self._MIP_row is None:
            raise UnsetData("Set MIP_row before use of CutScore.")
        if self._sage_to_solver_type is None:
            raise UnsetData("Set sage_to_solver_type before use of CutScore.")
        # To generate a seperator, we have to ensure, the function pi_parameters
        # is minimal.
        # parameters is given from some non linear solver and might not meet the conditions
        # of minimality.
        # validate_point will either give a point p = b,v which
        # we believe to up to rounding and L.C. that pi_p is a minimal function in the current cell
        # and satisfies the conditions of the model.
        # or will raise an error
        b, v = self.validate_point(parameters)
        self.set_feasible_point(b+v)
        pi = piecewise_function_from_breakpoints_and_values(b + [1], v + [0])
        row_data = self.get_MIP_row()
        sage_cut = [pi(fractional(QQ(bar_a_ij))) for bar_a_ij in row_data]
        sage_mip_obj =  [QQ(bar_cj) for bar_cj in self._MIP_objective]
        sage_result = self._cut_score.cut_score(sage_cut, sage_mip_obj)
        if self.get_prev_result() is not None:
            if abs(sage_result - self.get_prev_result())/sage_result < self._rel_tol:
                cut_gen_logger.debug(f"cutScore.__call__: Relalitve difference betweeen successive solutions is less than {self._rel_tol}. Stopping non-linear solver.")
                raise SolverHalt
            else:
                self.set_prev_result(sage_result)
        else:
            self.set_prev_result(sage_result)
        if self._timer is not None:
            if self._timer.solver_time_out():
                raise SolverTimeOut
        return self.get_sage_to_solver_type()(sage_result)

    @staticmethod
    def grad(parameters):
        """
        Graident of cutScore.
        """
        raise NotImplementedError
        # return self._cut_score.cut_score_grad(sage_cut, sage_mip_obj)

    @staticmethod
    def hess(parameters):
        raise NotImplementedError
        # return self._cut_score.cut_score_hess(sage_cut, sage_mip_obj)

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
        Defined in the cutGenerationSolver. This method is only indeded to be set and unset by solving
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

    def cut_obj_type(self):
         self._cut_obj_type

    def validate_point(self, point, poly_check=False):
        """
        Take parameters from the non linear solver and returns point that satasfies the minimal function model.
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
                raise SolverHalt(f"validate_point: breakpoint lambda_{i} < 0 or lambda_{i}>= 1; lambda_{i}=={b[i]}")
            if (v[i] < 0 or v[i] + epsilon > 1) and i != f_index:
                raise SolverHalt(f"validate_point: breakpoint gamma_{i} < 0 or gamma_{i}>= 1; gamma_{i}=={v[i]}")
        # pi(0) = 0, pi(f) = 1
        if abs(b[0]-0) <= epsilon:
            b[0] = 0
        else:
            cut_gen_logger.debug(f"validate_point: breakpoint lambda_0 >0, point: {point}, cell: {cell}")
            raise SolverHalt("breakpoint lambda_0 >0")
        if abs(v[0]-0) <= epsilon:
            v[0] = 0
        else:
            cut_gen_logger.debug(f"validate_point: breakpoint gamma_0 >0, point: {point}, cell: {cell}")
            raise SolverHalt("value gamma_0 >0")
        # bkpt[f_index] == f
        if abs(b[f_index] - f_trust) <= epsilon:
            b[f_index] = f_trust
        else:
            cut_gen_logger.debug(f"validate_point: breakpoint lambda_{f_index} != {f_trust}: point: {point}, cell: {cell}")
            raise SolverHalt(f"breakpoint lambda_{f_index} != {f_trust}")
        # pi_p(f) = 1
        if abs(v[f_index] - 1) <= epsilon:
            v[f_index] = 1
        else:
            cut_gen_logger.debug(f"validate_point: breakpoint lambda_{f_index} !=1: point: {point}, cell: {cell}")
            raise SolverHalt(f"value gamma_{f_index} != 1")
        # lipschitz constant and contunity.
        for i in range(n-1):
            if 0 < b[i+1]-b[i]<= epsilon:
                if abs(v[i+1]-v[i]) >= epsilon*M:
                    # potential discontunity
                    # not in (epsilon_i, M) - charts.
                    cut_gen_logger.debug(f"validate_point: Solution does not have lipschitz constant {M}: point: {point}, cell: {cell}")
                    raise SolverHalt(f"Solution does not have lipschitz constant {M}")
                # lambda_i+1 = lambda_i
                b[i+1] = b[i]
                # contunity, gamma_i =gamma_i+1
                if abs(v[i+1]-v[i]) > epsilon:
                    raise SolverHalt
                v[i+1] = v[i]
                # when epsilon <= v[i+1]-v[i] < epsilon*M
                # the solution exists in the intersection of the epsilon,M
                # can assume the values are correct/as intended
        # the last breakpoint should be distinct from 1 to enforce a breakpoint sequence.
        if 1-b[n-1] <= epsilon:
            cut_gen_logger.debug(f"validate_point: breakpoint lambda_{n-1} >= 1: point: {point}, cell: {cell}")
            raise SolverHalt(f"breakpoint lambda_{n-1} >= 1")
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
                            cut_gen_logger.debug(f"validate_point: {poly} evaluated at {b+v} == 0 when {poly} evaluated at {b+v} should be < 0: point: {point}, cell: {cell}")
                            raise SolverHalt(f"{poly} evaluated at {b+v} == 0 when {poly} evaluated at {b+v} should be < 0")
                        else:
                            pass
                    else:
                        cut_gen_logger.debug(f"validate_point: {poly} evaluated at {b+v} >= 0 when {poly} evaluated at {b+v} should be < 0: point: {point}, cell: {cell}")
                        raise SolverHalt(f"{poly} evaluated at {b+v} >= 0 when {poly} evaluated at {b+v} should be < 0")
            for poly in cell.le_poly():
                if poly(b+v) > 0:
                    cut_gen_logger.debug(f"validate_point: {poly} evaluated at {b+v} >0 when {poly} evaluated at {b+v} should be <=  0: point: {point}, cell: {cell}")
                    raise SolverHalt(f"{poly} evaluated at {b+v} >0 when {poly} evaluated at {b+v} should be <= 0")

            for poly in cell.eq_poly():
                if abs(poly(b+v)) >= epsilon or abs(poly(sage_point)) > sum(abs(grad(sage_point)) for grad in poly.gradient())*epsilon:
                    cut_gen_logger.debug(f"validate_point:{poly} evaluated at {b+v} != 0 when {poly} evaluated at {b+v} should be == 0: point: {point}, cell: {cell}")
                    raise SolverHalt(f"{poly} evaluated at {b+v} != 0 when {poly} evaluated at {b+v} should be == 0")
        return b,v


class Parallelism(abstractCutScore):
    """
    Normalized cut parallelism score.
    """
    def cut_score(cls, cut, mip_obj):
        obj_norm = vector(mip_obj).norm()
        cut_norm = vector(cut).norm()
        dot_product = vector(mip_obj).row()*vector(cut).column()
        return (dot_product[0]/(obj_norm*cut_norm))[0]


class SteepestDirection(abstractCutScore):
    """
    Steepest direction score.
    """
    def cut_score(cls, cut, mip_obj):
        dot_product = vector(mip_obj).row()*vector(cut).column()
        return dot_product[0][0]

# cutGenProblemParametersNames = ["algorithm", "cut_score", "max_num_bkpt", "multithread", "prove_seperator", "epsilon", "M", "max_cgf_solver_time", "max_cgf_solver_iter"]

class cutGenerationProblem:
    r"""
    A base class for interfacing with solvers and solving a cut generation problem.

    The cut generation problem is defined as max cutScore(cut) s.t. cut in cutSpace.

    The cutGenerationProblem options are listed below.

    Option: row algorithm - full; bkpt_as_param
    Option: num_bkpt - full: k <= max_bkpts
    Option: cut_score - parallelism, steepestdirection, scip, or custom
    Option: multithread - notImplemented
    Option: prove_seperator - bool
    Option: epsilon - value to det
    """
    def __init__(self, algorithm_name=None, cut_score=None, num_bkpt=None, multithread=False, prove_seperator=False, show_proof=False,
        epsilon=10**-7, M = 10**7, rel_tol=10**-6, max_cgf_solver_time=None, paramaterized_problem_solver=None,
        minimal_function_cache_path = None, backend=None):
        if algorithm_name is None or algorithm_name.lower() == "full":
            self._algorithm_name = "full"
            if num_bkpt is None or num_bkpt < 1 or num_bkpt > max_bkpts:
                raise ValueError(f"Incorrect number of breakpoints defined for full algorithm. 2 <= num_bkpt <= {max_bkpts}.")
            self._num_bkpt = num_bkpt
            self._cut_space = None
            self._solver = scipyCutGenProbelmSolverInterface
        elif algorithm_name.lower() == "bkpt_as_param":
            self._algorithm_name = "bkpt_as_param"
            self._solver = scipyCutGenProbelmSolverInterface
        elif algorithm_name.lower() == "value_poly_lp":
            self._algorithm_name = "value_poly_lp"
        else:
            raise ValueError("No other algorithms are supported at this time.")
        if cut_score is None:
            self._cut_score = cutScore(SteepestDirection)
        else:
            try:
                self._cut_score = cutScore(cut_score)
            except NameError:
                raise ValueError("Please provided a valid defined cutscore.")
        self._cut_score.set_sage_to_solver_type(self._solver.sage_to_solver_type)
        self._cut_space = None
        self._prove_seperator = prove_seperator
        self._show_proof = show_proof
        self._backend = backend
        self._espilon = epsilon
        self._M = M
        self._rel_tol = rel_tol
        self._cut_score.set_rel_tol(self._rel_tol)
        self._cut_score.set_espilon(self._espilon)
        self._cut_score.set_lipschitz_constant(self._M)
        self._max_cgf_solver_time = max_cgf_solver_time
        self._minimal_function_cache_path = minimal_function_cache_path

    def solve(self, binvarow, binvc, f):
        r"""Solves the paramaterized problem.

        Interperts the options and calls the correct solving algorithm.

        Passes any instructions to the underlying solver.

        """
        # assume MIP is a scip model; really we should be passing in and LP relaxation with variable information here.
        # The cut generation problem
        if self._algorithm_name == "full":
            cgf = self._algorithm_full_space(binvarow, binvc, f)
        elif self._algorithm_name == "bkpt_as_param":
            cgf = self._algorithm_bkpt_as_param(binvarow, binvc, f)
        elif self._algorithm_full_space == "value_poly_lp":
            cgf = self._algorithm_value_poly_lp(binvarow, binvc, f)
        return cgf

    def _algorithm_full_space(self, binvarow, binvc, f):
        r"""
        Solves the problem given a row of B^-1A and the reduced costs
        """
        self._cut_score.set_MIP_row(binvarow)
        self._cut_score.set_MIP_obj(binvc)
        self._cut_score.set_f_trust(f)
        frac_f = fractional(QQ(f))
        def cut_score(params):
            return self._cut_score(params)
        # if max_or_min == "max":
        best_value = -1*np.inf
        # else: #To do implement minimze options
        #     raise NotImplementedError
        solution_for_best_result = None
        if self._cut_space is None: # load the semi algebraic descriptions.
             if self._num_bkpt > max_bkpts:
                raise ValueError("The Minimal Functions Cache for {} breakpoints requested has not been computed.".format(self._num_bkpt))
             self._cut_space = PiMinContContainer(self._num_bkpt, backend=self._backend)
        # start the clock when the actual portion of the solving processs starts.
        problem_timer = cgfTimer(self._max_cgf_solver_time)
        self._cut_score.set_timer(problem_timer)
        for b, v in self._cut_space.get_rep_elems():
            # f is a bkpt when pi has a finite number of bkpts.
            # start by finding a bkpt sequence in the same cell
            # such that lambda_f_index = f holds.
            # pi(f) = 1;
            pi_test = piecewise_function_from_breakpoints_and_values(b+[1], v+[0])
            bsa_f_index = find_f_index(pi_test)
            self._cut_score.set_f_index(bsa_f_index)
            bkpt_bsa = nnc_poly_from_bkpt_sequence(b, backend=self._backend)
            lambda_f_index = bkpt_bsa.polynomial_map()[0].parent().gens()[bsa_f_index]
            bkpt_bsa.add_polynomial_constraint(lambda_f_index - frac_f, operator.eq)
            try:
                if not bkpt_bsa.upstairs().is_empty():
                    b0 = list(bkpt_bsa.find_point())
                    v0 = list(value_nnc_polyhedron_value_cords(b0, bsa_f_index).find_point())
                    # A feasible solution for cell problem has been found.
                    point = b0+v0
                    self._cut_score.set_feasible_point(point)
                    subdomain_with_f_constraint = bsa_of_rep_element(b0, v0)
                    lambda_f_index = subdomain_with_f_constraint.polynomial_map()[0].parent().gens()[bsa_f_index]
                    lhs =  lambda_f_index - frac_f
                    subdomain_with_f_constraint.add_polynomial_constraint(lhs, operator.eq)
                    self._cut_score.set_current_cell(subdomain_with_f_constraint)
                    subdomain_solver_constraints = self._solver.write_nonlinear_constraints_from_bsa(subdomain_with_f_constraint)
                    # Call a NL solver, attempt to solve the cell optimization problem.
                    try:
                        self._solver.nonlinear_solve(cut_score, point, subdomain_solver_constraints)
                    except SolverHalt:
                        pass
                    except SolverTimeOut:
                        # use the last good point to see if we have improvement before timing out our compuation
                        # and sending the result to the MIP
                        self._cut_score.set_timer(None)
                        point = self._cut_score.get_feasible_point()
                        if point is not None:
                            value_for_cell = self._cut_score(point)
                            if solution_for_best_result is None:
                                best_result = value_for_cell
                                solution_for_best_result = point
                                rep_elem_of_best_cell = b+v
                            if best_value < value_for_cell:
                                best_value = value_for_cell
                                solution_for_best_result = point
                                rep_elem_of_best_cell = b+v
                        break
                    # When a SolverHalt is encountered
                    # the NL solver has violated a constraint of the model or minimality
                    # within the cell. Use the last known feasible point from the
                    # solver.
                    point = self._cut_score.get_feasible_point()
                    try:
                        value_for_cell = self._cut_score(point)
                        continue_solving =  True
                    except SolverTimeOut:
                        self._cut_score.set_timer(None)
                        value_for_cell = self._cut_score(point)
                        continue_solving =  False
                    if solution_for_best_result is None:
                        best_result = value_for_cell
                        solution_for_best_result = point
                        rep_elem_of_best_cell = b+v
                    if best_value < value_for_cell:
                        best_value = value_for_cell
                        solution_for_best_result = point
                        rep_elem_of_best_cell = b+v
                    if not continue_solving:
                        break

            except EmptyBSA:
                pass
        # If result is None, the solver has failed to find any meaningful result or the computation has timed out.
        # There should always be a result and the SolverError should never be raised unless the computaion has timed out.
        if best_result is None:
            raise SolverError("The solver has failed, we should always get a result from the computation. Try increasing the time allowed for the solver to run.")
        val_result = [QQ(gamma_i) for gamma_i in solution_for_best_result[self._num_bkpt:]]
        bkpt_result = [QQ(lambda_i) for lambda_i in solution_for_best_result[:self._num_bkpt]]
        pi_p = piecewise_function_from_breakpoints_and_values(bkpt_result+[1],val_result+[0])
        cut_gen_logger.info(f"{pi_p} is the found function for the row: {binvarow}, objective:{binvc},and f{f}")
        if self._prove_seperator:
            res = minimality_test(pi_p) # add someway to log certificates.
            if not res:
                cut_gen_logger.error(f"minimality of  {pi_p}: {res}")
        return pi_p

    def _algorithm_bkpt_as_param(self, binvarow, binvc, f):
        """
        Solves the problem given a row of B^-1A and the reduced costs
        """
        self._cut_score.set_MIP_row(binvarow)
        self._cut_score.set_MIP_obj(binvc)
        self._cut_score.set_f_trust(f)
        problem_timer = cgfTimer(self._max_cgf_solver_time)
        self._cut_score.set_timer(problem_timer)
        frac_f = fractional(QQ(f))
        def cut_score(params):
            return self._cut_score(params)
        symmetrized_bkpts = [0, frac_f]
        # symmertized breakpoints should all be in [0,1)
        for b in binvarow:
            sage_b = fractional(QQ(b))
            b_sym = frac_f - sage_b
            if b_sym > 0:
                symmetrized_bkpts += [sage_b, b_sym]
            elif b_sym < 0:
                symmetrized_bkpts += [sage_b, 1+b_sym]
        symmetrized_bkpts = unique_list(symmetrized_bkpts)
        symmetrized_bkpts.sort()
        # it might be worth while to ensure if we have sufficient difference between breakpoints.
        sparse_bkpt = unique_list(sparse_enough_breakpoints(symmetrized_bkpts, self._espilon))
        if frac_f not in sparse_bkpt:
            sparse_bkpt.append(frac_f)
        num_bkpt = len(sparse_bkpt)
        if num_bkpt == 2:
            cut_gen_logger.info("Parsed row data suggests to use GMIC.")
            cut_gen_logger.info("Dim of value polyhedron: 0")
            pi_p =  gmic(frac_f)
            log_problem_result(sparse_bkpt, [0, 1], binvarow, binvc, f)
            if self._prove_seperator:
                # we always have a seperator here.
                cut_gen_logger.debug(f"The minimality of the found cgf is {True}")
            # return gmic, the feasible set for the optimization problem is a single point which corresponds to gmic.
            return pi_p
        # ensure a breakpoint sequence is given
        sparse_bkpt.sort()
        f_index = sparse_bkpt.index(frac_f)
        self._cut_score.set_f_index(f_index)
        value_polyhedron =  value_nnc_polyhedron(sparse_bkpt, f_index, backend=self._backend)
        cut_gen_logger.info(f"Dim of value polyhedron :{value_polyhedron.upstairs().ambient_dim()}")
        point = list(value_polyhedron.find_point())
        # initialize a feasible point for the cut scoring function to remember.
        self._cut_score.set_feasible_point(point)
        self._cut_score.set_current_cell(value_polyhedron)
        value_polyhedron_constraints = self._solver.write_linear_constraints_from_bsa(value_polyhedron)
        try:
            result = self._solver.nonlinear_solve(cut_score, point, value_polyhedron_constraints)
        except SolverHalt:
            pass
        except SolverTimeOut:
            self._cut_score.set_timer(None)
            if self._cut_score.get_feasible_point() is None:
                raise SolverError("The solver has failed to find a feasible point, allocate more solver time")
        self._cut_score.set_timer(None)
        result_point = self._cut_score.get_feasible_point()  # this should always be defined.
        val_result = [QQ(gamma_i) for gamma_i in point[num_bkpt:]]
        bkpt_result = [QQ(lambda_i) for lambda_i in point[:num_bkpt]]
        pi_p = piecewise_function_from_breakpoints_and_values(bkpt_result+[1], val_result+[0])
        log_problem_result(bkpt_result, val_result, binvarow, binvc, f)
        if self._prove_seperator:
            res = minimality_test(pi_p, self._show_proof) # add someway to log certificates.
            cut_gen_logger.info(f"Minimality of cgf: {res}")
        return pi_p

    def _algorithm_value_poly_lp(self, binvarow, binvc, f):
        raise NotImplementedError
        frac_f = fractional(QQ(f))
        symmetrized_bkpts = [0, frac_f]
        # symmertized breakpoints should all be in [0,1)
        for b in binvarow:
            sage_b = fractional(QQ(b))
            b_sym = frac_f - sage_b
            if b_sym > 0:
                symmetrized_bkpts += [sage_b, b_sym]
            elif b_sym < 0:
                symmetrized_bkpts += [sage_b, 1+b_sym]
        symmetrized_bkpts = unique_list(symmetrized_bkpts)
        num_bkpt = len(symmetrized_bkpts)
        if num_bkpt == 2:
            # put a logging step here.
            # return gmic, the feasible set for the optimization problem is a singlton which corresponds to gmic.
            return gmic(frac_f)
        symmetrized_bkpts.sort()
        f_index = symmetrized_bkpts.index(frac_f)
        # write value names
        coord_names = []
        val = [None]*(num_bkpt)
        for i in range(num_bkpt):
            coord_names.append('gamma'+str(i))
        K.gens()[0] == 0
        for i in range(1, num_bkpt):
            K.gens()[i] <=1
            K.gens()[i] > 0
        pi = piecewise_function_from_breakpoints_and_values(bkpt + [1], K.gens() + [0], merge=False)
        paramaterized_objective_function = sum(pi(QQ(b))*QQ(c) for b, c in zip(binvarow, binvc))
        # objective as a linear function of value parameters
        value_paramater_coeffs = paramaterized_objective_function.coefficients()
        lp_obj_in_solver_type =  [self._solver.sage_to_solver_type(obj_coeff) for obj_coeff in value_paramater_coeffs]
        val_poly = value_nnc_polyhedron_value_cords(bkpt, f_index)
        # now construct LP
        cons = self._solver.write_linear_constraints_from_bsa(val_poly)
        result = self._solver.lp_solve(cons, lp_obj_in_solver_type)
        vals = [QQ(v) for v in result]
        pi_p = piecewise_function_from_breakpoints_and_values(bkpt + [1], vals+ [0])
        log_problem_result(bkpt_result, val_result, binvarow, binvc, f)
        if self._prove_seperator:
            cut_gen_logger.info(f"Minimality of cgf: {res}")

        return pi_p

    def _algorithm_custom(self, binvarow, binvc, f):
        """
        Input: row data and cost data.
        Output: minimal function
        """
        # everthing that is written in one of this
        raise NotImplementedError


class abstractCutGenProblemSolverInterface:
    r"""
    Interfaces types from ``cutgeratingfunctionolgy`` to a specified solver.
    """
    def __init__():
        pass
    @staticmethod
    def write_linear_constraints_from_bsa(bsa):
        r"""
        Given a BSA with only linear constraints, converts the bsa object into a format that the underlying solver can use.
        """
        raise NotImplementedError

    @staticmethod
    def write_nonlinear_constraints_from_bsa(bsa):
        r"""
        Given a BSA with non linear constraints, converts the bsa object into a format that the underlying solver can use.
        """
        raise NotImplementedError


    @staticmethod
    def lp_solve(constraints, objective,  **solver_options):
        r"""
        Interface to solver's min/max f(x) s.t. Ax<=b.

        Should return optimal objective value, optimal objective solution, solver success, and solver_output
        """
        raise NotImplementedError

    @staticmethod
    def nonlinear_solve(constraints, objective, **solver_options):
        r"""
        Interface to solver's min/max f(x) s.t. p_i(x) <= b_i, where p_i is a polynomial and at least 1 p_i has degree larger than 1.
        """
        raise NotImplementedError


    @staticmethod
    def sage_to_solver_type(sage_ring_element):
        r"""

        """
        raise NotImplementedError


class scipyCutGenProbelmSolverInterface(abstractCutGenProblemSolverInterface):
    """
    Interfaces types and objects from ``cutgeneratingfunctiology`` to scipy.


    """
    @staticmethod
    def write_linear_constraints_from_bsa(bsa, epsilon=10**-9):
        r"""
        Given a BSA with only linear constraints, converts the bsa object into a format that the underlying solver can use.
        """
        #
        lb = []
        ub = []
        eq = []
        A_ieq = []
        A_eq = []
        for polynomial in bsa.eq_poly():
            if polynomial.degree() != 1:
                raise ValueError(f"Constraint {polynomial} == 0 is not linear.")
            linear_coeffs = [polynomial.coefficient(i) for i in polynomial.parent().gens()]
            A_eq.append(linear_coeffs)
            # polys in bsa are written poly op 0
            # rewrite to correct scipy notation.
            constant = -1*polynomial.constant_coefficient()
            eq.append(constant)
        for polynomial in bsa.lt_poly():
            if polynomial.degree() != 1:
                raise ValueError(f"Constraint {polynomial} < 0 is not linear.")
            linear_coeffs = [polynomial.coefficient(i) for i in polynomial.parent().gens()]
            A_ieq.append(linear_coeffs)
            # mimic in a non-rigrous way <
            constant = -1*polynomial.constant_coefficient() - epsilon
            lb.append(-np.inf)
            ub.append(constant)
        for polynomial in bsa.le_poly():
            if polynomial.degree() != 1:
                raise ValueError(f"Constraint {polynomial} < 0 is not linear.")
            linear_coeffs = [polynomial.coefficient(i) for i in polynomial.parent().gens()]
            A_ieq.append(linear_coeffs)
            constant = -1*polynomial.constant_coefficient()
            lb.append(-np.inf)
            ub.append(constant)
        lb = np.array(lb)
        ub = np.array(ub)
        if len(lb) == 0 and len(eq) > 0 :
            LinearConstraint(np.array(A_eq), eq, eq)
        elif len(lb) >0  and len(eq) > 0:
            return [LinearConstraint(np.array(A_eq), eq, eq), LinearConstraint(np.array(A_ieq), lb, ub)]
        elif len(lb) >0  and len(eq) == 0:
            return LinearConstraint(np.array(A_ieq), lb, ub)
        else:
            raise ValueError(f"No constraints have been written.")
    @staticmethod
    def write_nonlinear_constraints_from_bsa(bsa, epsilon=10**-9):
        r"""
        Given a BSA with nonlinear constraints, converts into an equivlent set of nonlinear constraints for scipy.

        Treats p(x) < c as p(x) + epsilon <= c for all epsilon>0.
        """
        nonlinear_constraints = []
        # All variables are implicitly bounded between 0 and 1.
        # We should establish using a lower bound.
        # This section can be improved. Hessians need to be rewritten to have the right signature.
        for polynomial in bsa.eq_poly():
            def poly(array_like):
                # map coordinates names in BSA to coordinates of solvers
                input_map = {polynomial.parent().gens()[i]: array_like[i] for i in range(polynomial.parent().ngens())}
                return np.array([polynomial.subs(input_map)])
            def poly_grad(array_like):
                input_map = {polynomial.parent().gens()[i]: array_like[i] for i in range(polynomial.parent().ngens())}
                return np.array([partial.subs(input_map) for partial in polynomial.gradient()])
            # def poly_hess(array_like):
                # input_map = {polynomial.parent().gens()[i]: array_like[i] for i in range(polynomial.parent().ngens())}
                # return np.array([[second_partial.subs(input_map) for second_partial in partial.gradient()]  for partial in polynomial.gradient()]])
            nonlinear_constraints.append(NonlinearConstraint(poly, 0, 0, jac=poly_grad))
        for polynomial in bsa.le_poly():
            def poly(array_like):
                input_map = {polynomial.parent().gens()[i]: array_like[i] for i in range(polynomial.parent().ngens())}
                return np.array([polynomial.subs(input_map)])
            def poly_grad(array_like):
                input_map = {polynomial.parent().gens()[i]: array_like[i] for i in range(polynomial.parent().ngens())}
                return np.array([partial.subs(input_map) for partial in polynomial.gradient()])
            # def poly_hess(array_like):
                # input_map = {polynomial.parent().gens()[i]: array_like[i] for i in range(polynomial.parent().ngens())}
                # return np.array([[second_partial.subs(input_map) for second_partial in partial.gradient()]  for partial in polynomial.gradient()]])
            nonlinear_constraints.append(NonlinearConstraint(poly, -np.inf, 0,  jac=poly_grad))
        for polynomial in bsa.lt_poly():
            def poly(array_like):
                input_map = {polynomial.parent().gens()[i]: array_like[i] for i in range(polynomial.parent().ngens())}
                return np.array([polynomial.subs(input_map)+epsilon])
            def poly_grad(array_like):
                input_map = {polynomial.parent().gens()[i]: array_like[i] for i in range(polynomial.parent().ngens())}
                return np.array([partial.subs(input_map) for partial in polynomial.gradient()])
            # def poly_hess(array_like):
                # input_map = {polynomial.parent().gens()[i]: array_like[i] for i in range(polynomial.parent().ngens())}
                # return np.array([[second_partial.subs(input_map) for second_partial in partial.gradient()]  for partial in polynomial.gradient()]])
            nonlinear_constraints.append(NonlinearConstraint(poly, -np.inf, 0,  jac=poly_grad))
        return nonlinear_constraints


    @staticmethod
    def lp_solve(objective, constraints, **solver_options):
        r"""
        Interface to solver's min/max f(x) s.t. Ax<=b.

        Should return optimal objective value, optimal objective solution.
        """
        raise NotImplementedError

    @staticmethod
    def nonlinear_solve(objective, x0, cons, jac=None, hess=None,  **solver_options):
        r"""
        Given converted constraints and an objective function that is compitable with the solver,
        use scipy minimize to ...
        """
        if jac is not None:
            if hess is not None:
                result =  minimize(objective, x0, constraints=cons, jac=jac, hess=hess)
            else:
                result = minimize(objective, x0, constraints=cons, jac=jac)
        else:
            result = minimize(objective, x0, constraints=cons, jac=jac, hess=hess)
        return result.success, result.fun, result.x, result

    @staticmethod
    def sage_to_solver_type(sage_ring_element):
        """
        scipy supports inputs of floats, convert sage field element to its equivlant numerical (python) floating point value.
        """

        # be lazy and assume the_ring_element is something the converts to a rational number (or can be put into a floating point approximation).
        # this is a point where we lose the exactness of sage.
        return float(sage_ring_element)


class OptimalCut(Sepa):
    def __init__(self, algorithm_name=None, cut_score=None, num_bkpt=None, multithread=False, prove_seperator=True, show_proof = False, epsilon=10**-7, M = 10**7):
        self.ncuts = 0
        self.cgp = cutGenerationProblem(algorithm_name=algorithm_name, cut_score=cut_score, num_bkpt=num_bkpt, multithread=multithread, prove_seperator=prove_seperator, show_proof = show_proof, epsilon=epsilon, M = epsilon)
    def getOptimalCutFromRow(self, cols, rows, binvrow, binvarow, primsol, pi_p):
        """ Given the row (binvarow, binvrow) of the tableau, computes optimized cut.

        :param primsol:  is the rhs of the tableau row
        :param cols:     are the variables
        :param rows:     are the slack variables
        :param binvrow:  components of the tableau row associated to the basis inverse
        :param binvarow: components of the tableau row associated to the basis inverse * A
        :param 1dPWL:    a minimal funciton element of PiMin<=k with pi_p(primsol)=1

        The interesection cut is given by
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
                # From matthias, do the supperattitive portion of the function about 0.
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
                # get the row of B^-1 for this basic integer variable with fractional solution value
                binvrow = scip.getLPBInvRow(i)

                # get the tableau row for this basic integer variable with fractional solution value
                binvarow = scip.getLPBInvARow(i)

                # get current reduced costs for objective evaluation.
                costs = [scip.getColRedCost(j) for j in cols if j not in basisind]

                cgf = self.cgp.solve(binvarow, costs, primsol) # produce an optimal cgf

                cutcoefs, cutrhs = self.getOptimalCutFromRow(cols, rows, binvrow, binvarow, primsol, cgf)

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

                    infeasible = scip.addCut(cut, forcecut=True)
                    self.ncuts += 1

                    if infeasible:
                       result = SCIP_RESULT.CUTOFF
                    else:
                       result = SCIP_RESULT.SEPARATED
                scip.releaseRow(cut)

        return {"result": result}
