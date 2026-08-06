"""
Define and solve cut generation problems.
"""

from cutgeneratingfunctionology.igp import *
from .execptions import *
from .scipy_cgp_interface import *
from .cgf_specializations import * # imports cut generating functionololgy and specalizations
from scipy.optimize import minimize as scipy_minimize
from cvxpy import Variable as cvxpy_Variable, Maximize as cvxpy_Maximize, Minimize as cvxpy_Minimize, Problem as cvxpy_Problem
import logging
import time


cut_generation_problem_logger = logging.getLogger(__name__)
cut_generation_problem_logger.setLevel(logging.ERROR)

# minimal_exprction_cashe_logging = True

#avail_rep_elems = minimal_exprction_cache_info()["avail_rep_elems"]

#class breakpoint_sequence(list):
#    pass

# Helper functions/model enforcement

def find_f_index(min_pwl):
    r"""
    Assume a model function, pi_p, with a finite number of breakpoints is given.
    Finds the index i such that pi_p(lambda_i) = 1.

    INPUT:
    - igp piecewise linear function

    OUTPUT:
    - integer
    """
    return min_pwl.end_points().index(find_f(min_pwl))


def infer_mutliplicity_vector_and_proper_breakpoints(bkpt, epsilon):
    r"""
    Deduces a multiplicity vector sum_{j=0}^{k-1} m_j of length k such that  in which |b_j - sum_{l=1}^m_j b_{j, epsilon_l}| < m_j * epsilon. 
    
    This gives the possiblity of at most 2k propper breakpoints in the solution of the cut generation problem. 

    INPUT:
    - breakpoint sequence

    OUTPUT:
    - breakpoint sequence such that bkpt[i+1]-bkpt[i] > epsilon xor bkpt[i] = bkpt[i+1].

    TESTS::
    >>> from parametricCutGen.cut_generation_problem import *
    >>> infer_mutliplicity_vector_and_proper_breakpoints([0, 10**-7, 4/5], 10**-6)
    [0, 4/5]
    """
    pbkpt = list(tuple(bkpt)) # cheap way of deep copying lists
    pbkpt.sort()
    for i in range(len(pbkpt)-1):
        if abs(pbkpt[i]) < epsilon:
            pbkpt[i] = 0
        # 1 equiv 0 mod 1
        elif abs(pbkpt[i] - 1) < epsilon:
            pbkpt[i] = 0
        elif abs(pbkpt[i] - pbkpt[i+1]) < epsilon:
            pbkpt[i+1] = pbkpt[i]
    if abs(pbkpt[len(pbkpt)-1] -1) < epsilon or abs(pbkpt[len(pbkpt)-1])< epsilon:
        pbkpt[len(pbkpt)-1] = 0
    return pbkpt


def symmetrize_about_f_mod_1(bkpt, f):
    """
    Add symmerties about f modulo 1.
    Retun a list such for all i in bkpt, there exist a j with sym_bkpt[j] equiv bkpt[i] - f mod 1.
    """
    frac_f = fractional(QQ(f))
    symmetrized_bkpts = [QQ(0), frac_f]
    # symmertized breakpoints should all be in [0,1)
    for b in bkpt:
        sage_b = fractional(QQ(b))
        b_sym = frac_f - sage_b
        if b_sym > 0:
            symmetrized_bkpts += [sage_b, b_sym]
        elif b_sym < 0:
            symmetrized_bkpts += [sage_b, 1+b_sym]
    symmetrized_bkpts = unique_list(symmetrized_bkpts)
    symmetrized_bkpts.sort()
    return symmetrized_bkpts

def PWL_with_bkpts_manifold_chart_constraints(bkpt, M, values=None, *, coeff_type=None, backend=None, feasiblity_problem=False, ring=QQ, backend_kwds={"sage_mip":{"solver":"GLPK", "maximization":True}}, backend_args={"sage_mip":{"mip":None}}):
    """Defines constraints for U_{m, M}, the domain of the selected chart from the cut generation problem.

    :bkpt: - list, a breakpoint sequence, entries are elements of QQ.
    
    :M: - ring(M) in QQ. 

    :ring:
    - `sage.Ring`
    - `None` -  
    
    :backend:
    Speicfy the polyhedral backend. `None` is assumed to be ppl. 
    - `None` - Var is ppl.Variable
    - \"pplite\" - Var is pplite.Variable
    - \"sage_mip\" - Var is sage.numerical.mip.MIPVariable
    
    :coeff_type:
    :backend:'s' constraint coeffiecent type. 
    Possible values:
    - \'int\' - python.int
    - \'float\' - python.float
    - `None` - sage.Ring, default sage.Ring.Rings.Rational=QQ, not checked.

    :feasiblity_problem: 
    Write constriants for solving the feasiblity problem min  x_{n} st. x in set((U_{m,M},  x_{n}), x_{n} >=0).  
    Possible Values:
    - bool
    """
    n = len(bkpt)
    if values is None:
        if backend == 'pplite':
            Var = pplite_Variable
            coeff_type='int'
            values = [Var(i) for i in range(n)]
        elif backend == 'sage_mip':
            if backend_args["sage_mip"]["mip"] is not None:
                backend_args["sage_mip"]["mip"] = mip
                values = mip.default_variable()
                coeff_type='ring'
            else:
                raise ValueError("please provide a MixedIntegerProgram")
        else:
            Var = Variable
            coeff_type='int'
            values = [Var(i) for i in range(n)]
    R=ring
    
    cons = []
    if not feasiblity_problem:
        # assert len(values) == n
        if coeff_type=='int':
            M = QQ(M)
            for i in range(n-1):
                lcd = lcm([bkpt[i+1].denominator(), bkpt[i].denominator(), M.denominator()])
                cons.append( int(lcd ) * (values[i+1] - values[i]) - int(lcd * M * (bkpt[i+1]-bkpt[i])) <= 0 )
                cons.append( int(-1*lcd * M * (bkpt[i+1]-bkpt[i])) +  (-1*int(lcd) * (values[i+1] - values[i]))  <= 0 )
            lcd = lcm(bkpt[n-1].denominator(),  M.denominator())
            cons.append( int(lcd) * (1 - values[n-1]) - int(lcd * M* (1 - bkpt[n-1])) <= 0 )
            cons.append( int(-1*lcd * M * (1 - bkpt[n-1])) + (-1* int(lcd )) * (1 - values[n-1]) <= 0 ) # correct
        elif coeff_type=='float':
            M=float(M)
            for i in range(n-1):
                cons.append( (values[i+1] - values[i]) - M* float(bkpt[i+1]-bkpt[i]) <= 0 )
                cons.append( -1*M* float(bkpt[i+1]-bkpt[i]) -  (values[i+1] - values[i])  <= 0 )
            cons.append(  (1 - values[n-1]) - M * float(1 - bkpt[n-1]) <= 0 )
            cons.append( -1*M* float(1 - bkpt[n-1]) +  (-1 * (1 - values[n-1])) <= 0 )    
        else:
            M=R(M)
            # bkpt is assumed to be sage rational type,
            for i in range(n-1):
                cons.append( (values[i+1] - values[i]) + (-1*M)* (bkpt[i+1]-bkpt[i]) <= 0 )
                cons.append( M * (bkpt[i+1]-bkpt[i]) -  (values[i+1] - values[i])  <= 0 )
            cons.append( (1 - values[n-1]) - M * (1 - bkpt[n-1]) <= 0 )
            cons.append( M * float(1 - bkpt[n-1]) - (1 - values[n-1]) <= 0 )
    else:
        # assert len(values) == (n+1)
        cons.append(values[n] >= 0)
        if coeff_type=='int':
            M = QQ(M)
            for i in range(n-1):
                lcd = lcm([bkpt[i+1].denominator(), bkpt[i].denominator(), M.denominator()])
                cons.append( int(lcd ) * (values[i+1] - values[i]) - int(lcd * M * (bkpt[i+1]-bkpt[i])) - values[n] <= 0 )
                cons.append( int(-1*lcd * M * (bkpt[i+1]-bkpt[i])) +  (-1*int(lcd) * (values[i+1] - values[i])) - values[n] <= 0 )
            lcd = lcm(bkpt[n-1].denominator(),  M.denominator())
            cons.append( int(lcd) * (1 - values[n-1]) - int(lcd * M* (1 - bkpt[n-1])) - values[n] <= 0 )
            cons.append( int(-1*lcd * M * (1 - bkpt[n-1])) + (-1* int(lcd )) * (1 - values[n-1]) - values[n] <= 0 ) # correct
        elif coeff_type=='float':
            M=float(M)
            for i in range(n-1):
                cons.append( (values[i+1] - values[i]) - M* float(bkpt[i+1]-bkpt[i]) - values[n] <= 0 )
                cons.append( -1*M* float(bkpt[i+1]-bkpt[i]) -  (values[i+1] - values[i]) - values[n]  <= 0 )
            cons.append(  (1 - values[n-1]) - M * float(1 - bkpt[n-1]) - values[n] <= 0 )
            cons.append(-1* M* float(1 - bkpt[n-1]) +  (- 1 + values[n-1]) - values[n] <= 0 )    
        else:
            M=R(M)
            for i in range(n-1):
                cons.append( (values[i+1] - values[i]) + (-1*M)* (bkpt[i+1]-bkpt[i]) - values[n] <= 0 )
                cons.append( -1*M * (bkpt[i+1]-bkpt[i]) -  (values[i+1] - values[i]) - values[n] <= 0 )
            cons.append( (1 - values[n-1]) - M * (1 - bkpt[n-1]) - values[n] <= 0 )
            cons.append(-1* M * (1 - bkpt[n-1]) + (-1 + values[n-1]) - values[n] <= 0 )
    return cons

def find_feasible_point(bkpt, M, f_index,  *, epsilon=1e-9, backend=None, \
    ring=QQ, backend_kwds={"sage_mip":{"solver":"GLPK", "maximization":True}}, \
    backend_args={"sage_mip":{"mip":None}}): 
    """
    Find x0 in VP_bkpt,f_index cap U_{m, M} where m is the multiplcity vector of bkpt.

    See doc strings of value_nnc_polyhedron_constraints and PWL_with_bkpts_manifold_chart_constraints for parameter infromation. 
    
    returns: 
    `x0` - np.array of dtype=float of shape (n,)
    """
    n = len(bkpt)
    if backend == 'pplite' or backend == None:
        if backend == 'pplite':
            Value_Poly = pplite_NNC_Polyhedron(dim_type = n, spec_elem = "universe", topology = "nnc")
            bsa = BasicSemialgebraicSet_polyhedral_pplite_NNC_Polyhedron(n, Value_Poly)
            values = [pplite_Variable(i) for i in range(n)]
        else:
            Value_Poly = NNC_Polyhedron(n)
            bsa = BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(n, Value_Poly)
            values = [Variable(i) for i in range(n)]
        minimality_cons = value_nnc_polyhedron_constraints(bkpt, f_index, values, epsilon=epsilon, coeff_type='int')
        manifold_cons = PWL_with_bkpts_manifold_chart_constraints(bkpt, M, values, coeff_type='int')
        for con in minimality_cons+manifold_cons:
            try:
                Value_Poly.add_constraint(con)
            except Exception as e:
                cut_generation_problem_logger.warning(f"constraint:{con} could not be added to value Polyhedron. Reason: {e}")
                continue
        x0 = bsa.find_point()
    elif backend == 'sage_mip' or backend=='cvxpy':
        if backend == 'sage_mip':
            # try with highs. 
            bsa = BasicSemialgebraicSet_polyhedral_MixedIntegerLinearProgram(QQ, n+1, solver=backend_kwds["sage_mip"]["solver"])
            values = list(bsa.mip_gens())
            lp = bsa.mip() # has method.add_constraint
            lp.set_objective(-1*values[n]) # The value polyhedron is always feasible. 
            minimality_cons_feasiblity_problem = value_nnc_polyhedron_constraints(bkpt, f_index, values, epsilon=epsilon, coeff_type='ring', feasiblity_problem=True)
            manifold_cons_feasiblity_problem = PWL_with_bkpts_manifold_chart_constraints(bkpt, M, values, coeff_type='ring', feasiblity_problem=True)
            for con in minimality_cons_feasiblity_problem+manifold_cons_feasiblity_problem:
                try:
                    lp.add_constraint(con)
                except Exception as e:
                    cut_generation_problem_logger.warning(f"constraint:{con} could not be added to value Polyhedron. Reason: {e}")
                    continue
            cut_generation_problem_logger.debug(f"constraints:{minimality_cons_feasiblity_problem+manifold_cons_feasiblity_problem}")
            lp.solve()
            x0 = lp.get_values(values)[:-1]
        elif backend == 'cvxpy':
            values = cvxpy_Variable(n+1)
            minimality_cons = value_nnc_polyhedron_constraints(bkpt, f_index, values, epsilon=epsilon, coeff_type='float', feasiblity_problem=True)
            manifold_cons = PWL_with_bkpts_manifold_chart_constraints(bkpt, M, values, coeff_type='float', feasiblity_problem=True)
            obj = cvxpy_Minimize(values[n])
            prob = cvxpy_Problem(obj, minimality_cons+manifold_cons)
            result=prob.solve()
            print(result)
            x0 = list(values.value)[:-1]
    else:
        raise ValueError("specify a valid backend for the feasiblity problem")
    return np.array([float(xi) for xi in x0])

def write_cgp_constraints(bkpt, M, f_index,  *, epsilon=1e-9, backend='pplite'):
    """
    Write the constraints of a cut genration problem for the value polyhedron and the manifold constraints.

    :bkpt:, :M:, :f_index: are defined in :function:``value_nnc_polyhedron_constraints`` and :function:`PWL_with_bkpts_manifold_chart_constraints.
    
    :epsilon: real number >= 0, models > as >=.

    :backend: Polyhedral backend; `None` or \'pplite\'. Default is \'pplite\' as its implemenation is faster. None uses pplpy as a baceknd.
    
    Outputs 

    :bsa: - :cutgeneratingfunctiology.spam:``BasicSemialgebraicSet`` with the given constraints. 
    
    """
    n = len(bkpt)
    if backend == 'pplite':
        Value_Poly = pplite_NNC_Polyhedron(dim_type = n, spec_elem = "universe", topology = "nnc")
        bsa = BasicSemialgebraicSet_polyhedral_pplite_NNC_Polyhedron(n, Value_Poly)
        values = [pplite_Variable(i) for i in range(n)]
    else:
        Value_Poly = NNC_Polyhedron(n)
        bsa = BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(n, Value_Poly)
        values = [Variable(i) for i in range(n)]
    minimality_cons = value_nnc_polyhedron_constraints(bkpt, f_index, values, epsilon=epsilon, coeff_type='int')
    manifold_cons = PWL_with_bkpts_manifold_chart_constraints(bkpt, M, values, coeff_type='int')
    for con in minimality_cons+manifold_cons:
        try:
            Value_Poly.add_constraint(con)
        except Exception as e:
            cut_generation_problem_logger.warning(f"constraint:{con} could not be added to value Polyhedron. Reason: {e}")
            continue
    return bsa

def inf_norm_of_cont_pwl(f, g):
    """
    ||f-g||_infty where f,g are PW1D functions. 
    """
    return max([abs(v) for v in (f-g).values_at_end_points()])


def log_problem_result(bkpt, val, binvarow, binvc, f):
    cut_generation_problem_logger.info(f"Cut generation problem solved.")
    cut_generation_problem_logger.info(f"breakpoitns: {bkpt}")
    cut_generation_problem_logger.info(f"values: {val}")
    cut_generation_problem_logger.info(f"row: {binvarow}")
    cut_generation_problem_logger.info(f"objective:{binvc}")
    cut_generation_problem_logger.info(f"f: {f}")
            

# cutGenProblemParametersNames = ["algorithm", "cut_score", "max_num_bkpt", "multithread", "prove_seperator", "epsilon", "M", "max_cgp_solver_time", "max_cgf_solver_iter"]

class cutGenerationProblem:
    r"""
    A class for solving a cut generation problem (cgp). 
    
    This class combines the manifold charts of the space of piecewise linear functions with the constraints on minimal functions to optimize cuts over the subspaces of minimal functions. 

    A chart, phi: U -> PWL_{<=n} has a multiplicity vector, m, and a liptishitz constant, M, associated with it. 
    The multiplicity vector is infrered by the parameter :epsilon:. 
    The liptishitz constant M is defined by :M:.

    The chart domain is intersected with the constraints of minimality.
    
    This gives a subset F of RR^n in which phi(F) <= PiMin_{<=} (containment). 

    Assume MIP lp relaxation data is provided.
    Define C_F as the set of cuts which there exists pi in phi(F) such that cut is defined by pi (as an interesction cut). 
    F is specified by :algorithm:. 

    Let s: C_F -> RR be a cut scoring hueristic and assume that s is C^2 on C_F. s is defined by :cut_score:.

    A cut generation problem is defined as max s(c) s.t. c in C_F. 

    The solution to a cut generating problem is reported as f in PWL_{<=n}. 
    
    Spaces of minimal functions:

    VP_{b, f_index} - Value Polyhedron_{b, f_index} 
    
    M - bigcup_(alpha in min function cache) S_alpha = {p in RR^k : pi_p is minimal if and only if pi_alpha is minial} 
    
    The cgp optional keywords are listed below.

    :algorithm: - string
    :backend: - `None` or string
    :cut_score: - string
    :epsilon: - types convertable to sage rational 
    :M: -  types convertable to sage rational
    :max_cgp_solver_time: - NotImplemented
    :max_num_of_bkpts: - integer
    :multithread: - NotImplemented
    :prove_seperator: - bool
    :rel_tol: NotImplemted
    :show_proof: - bool
    :solver_args: - dict
    :solver_kwds: - dict
    :enable_profiling: - bool

    PARAMETER DETAILS

    :algorithm: 

    - \"full\" : (pre-alpha) Takes F=M. 
    
    - \"bkpt\_as\_param\" : (alpha) Takes F = U cap VP_{b, b.index(f)}. 
    The chart :epsilon: and the row data are used to detemine a viable multiplicity vector.
    Symmetries about f=fractional(overline(b_i))) to ensure the VP_{b, b.index(f)} has dimension higher than 0.
    
    - \"value\_poly\_vert" (alpha) Takes F = vert(U cap VP_{b, b.index(f)}).

    :backend: 

    - `None` : (stable) use ppl as a polyhedral backend for computations.
    
    - \"pplite\" : (alpha) use pplite as a polyhedral backend for computations. 

    - \"sage_mip\" : (stable) use sagemaths `MixedIntegerProgram` class a polyhedral backend for computations.

    :cut_score:

    - \"parallelisim\" : s(cut) = <cut,obj>/||cut - obj||

    - \" \" 
    """
    def __init__(self, *, algorithm=None, backend=None, cut_score=None, 
        epsilon=1e-2, numerical_epsilon=1e-9, M = 1e6, 
        max_cgp_solver_time=None, max_num_of_bkpts=4, multithread=False, 
        objective_sense="maximize", prove_seperator=False, show_proof=False, 
        solver_args={"cvxpy": None, "scipy": None}, solver_kwds={"cvxpy:" : None, "scipy": None}, 
        backend_kwds={"sage_mip":{"solver":"GLPK"}},  enable_profiling=False):
        """
        TESTS::
        >>> from parametricCutGen.cut_generation_problem import *
        >>> cgp = cutGenerationProblem()
        False
        >>> cgp.get_cgp_input_parameters()
        
        """
        # To recall parameters for debugging puposes. 
        self._cgp_input_parameters = {**locals()} # copy input parameters; only inital inputs here.
        self._cgp_input_parameters.pop("self") # This is required to ensure that when parameters are reused for initlaziation, problems doens't arise.
        # TODO: Implement a way so that cgps with the same parameter are the same object; maybe.
        self._objective_sense = objective_sense
        if self._objective_sense not in ["maximize", "minimize"]:
            raise ValueError("Objective senese should be either maximize or minimize.")
        if algorithm is None or algorithm.lower() == "full":
            self._algorithm = "full"
            self._max_num_of_bktps = max_num_of_bkpts
            if max_num_of_bkpts not in avail_rep_elems:
                cut_generation_problem_logger.ERROR(f"Function cache {max_num_of_bkpts} is not available; generation of functions at runtime will have long runtime.")
        elif algorithm.lower() == "bkpt_as_param":
            self._algorithm = "bkpt_as_param"
        elif algorithm.lower() == "value_poly_vert":
            self._algorithm = "value_poly_vert"
        else:
            raise ValueError("No other algorithms are supported at this time.")
        self._cut_score_generator = scipyCutScoreForBkptAsParam(cut_score, epsilon, M)
        self._prove_seperator = prove_seperator
        self._show_proof = show_proof
        self._backend = backend
        self._epsilon = epsilon
        self._numerical_epsilon = numerical_epsilon 
        self._M = M
        self._max_cgp_solver_time = max_cgp_solver_time
        self._max_num_of_bkpts = max_num_of_bkpts
        self._solver_args = solver_args
        self._solver_kwds = solver_kwds
        self._backend_kwds= backend_kwds
        self._enable_profiling = enable_profiling
        if self._enable_profiling:
            import cProfile
            self._pr = cProfile.Profile()
    def __repr__(self):
        return f"Cut generation problem using algorithm {self._algorithm}."

    def solve(self, binvrow, binvarow, f, costs, cols=None, rows=None, scip=None):
        r"""Solves the parameterized problem.

        Interprets the options and calls the correct solving algorithm.

        Passes any instructions to the underlying solver.
        """
        if self._enable_profiling:
            self._pr.enable()
        if self._algorithm == "full":
            result = self._algorithm_full_space(binvrow, binvarow, f, costs, cols, rows, scip)
        elif self._algorithm == "bkpt_as_param":
            result = self._algorithm_bkpt_as_param(binvrow, binvarow, f, costs, cols, rows, scip)
        elif self._algorithm == "value_poly_vert":
            result = self._algorithm_value_poly_lp(binvrow, binvarow, f, costs, cols, rows, scip)
        if self._enable_profiling:
            self._pr.disable()
        return result

    def _algorithm_full_space(self, binvrow, binvarow, f, costs, cols=None, rows=None, scip=None):
        r"""
        Solves the problem given a row of B^-1A and the reduced costs
        """
        raise NotImplementedError("Full space algorithm is under development.")

    def _algorithm_bkpt_as_param(self, binvrow, binvarow, f, costs, cols=None, rows=None, scip=None):
        """
        Solves the problem given a row of B^-1A and the reduced costs
        """
        # Ensure numerical difference btween breakpoints
        sparse_bkpt = unique_list(infer_mutliplicity_vector_and_proper_breakpoints([fractional(QQ(b)) for b in binvrow], self._epsilon))
        sparse_bkpt.sort()
        bkpt = symmetrize_about_f_mod_1(sparse_bkpt, f)
        n = len(bkpt)
        # ensure a breakpoint sequence is given
        f_index = bkpt.index(fractional(QQ(f)))
        cut_generation_problem_logger.debug(f"bkpt={bkpt}\nsparce_bkpt={sparse_bkpt}\nf_index={f_index}")    
        #cut_generation_problem_logger.debug(f"gen poly")
        bsa = write_cgp_constraints(bkpt, self._M, f_index, backend=self._backend)
        if self._backend == "sage_mip":
            x0 = find_feasible_point(bkpt, self._M, f_index, epsilon=self._numerical_epsilon, backend=self._backend) # solves feasiblity lp, rather than find an interor point of a polyehdron.
        else:
            x0 = np.array([float(x) for x in bsa.find_point()]) # picks inital point on relitive interor of the value polyhedron.
        scipy_cons = map_polyhedral_bsa_to_scipy_LinearConstraint(bsa)
        if scip is not None:
            cutcoefs_expr, cutrhs_expr, integral_indices, lp_soln = self.getParamaterizedCutExpr(scip, cols, rows, binvrow, binvarow, f, bkpt, f_index)
            s = self._cut_score_generator.gen_cut_score(f_index,  cutcoefs_expr, cutrhs_expr, lp_soln, costs, integral_indices)
            if self._objective_sense == "maximize":
                def cut_score(value_parameters):
                    return -1* s(value_parameters)
            else:
                 def cut_score(value_parameters):
                    return s(value_parameters)
        else:
            raise NotImplementedError
        result = scipy_minimize(cut_score, x0, constraints=scipy_cons)
        #cut_generation_problem_logger.debug(f"fin")
        score = result.fun
        b = bkpt
        v = to_sage_rationals(result.x)
        pi_p = piecewise_function_from_breakpoints_and_values(b+[1], v+[0])
        log_problem_result(b, v, binvarow, costs, f)
        if self._prove_seperator:
            res = minimality_test(pi_p, self._show_proof) # add someway to log certificates.
            cut_generation_problem_logger.info(f"Minimality of cgf: {res}")
        return pi_p, score, n

#        if self._backend == 'cvxpy':
#            raise NotImplementedError("cvxpy backed not impemented")
#
    def _algorithm_value_poly_lp(self, binvrow, binvarow, f, costs, cols=None, rows=None, scip=None):
        problem_timer = cgpTimer(self._max_cgp_solver_time)
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
        symmetrized_bkpts.sort()
        # it might be worth while to ensure if we have sufficient difference between breakpoints.
        model_sparsity = max(float(1/self._max_num_of_bkpts), self._epsilon)
        sparse_bkpt = unique_list(infer_mutliplicity_vector_and_proper_breakpoints(symmetrized_bkpts, model_sparsity))
        if frac_f not in sparse_bkpt:
            sparse_bkpt = [x for x in sparse_bkpt if abs(frac_f - x) > self._epsilon]
            sparse_bkpt.append(frac_f)
        num_bkpt = len(sparse_bkpt)
        # Ensure a breakpoint sequence is given.
        sparse_bkpt.sort()
        # Data bkpt sequences need to be sage types.
        sparse_bkpt = to_sage_rationals(sparse_bkpt)
        f_index = sparse_bkpt.index(frac_f)
        # lets use cvxpy to solve the problem directly
        cvxpy_vals = cvxpy_Variable(len(sparse_bkpt)) # aka the values
        minimality_cons = value_nnc_polyhedron_constraints(sparse_bkpt, f_index, cvxpy_vals, coeff_type='float', epsilon=self._epsilon)
        
        map_sage_expr_to_cvpxy_expr = lambda expr : sum( float(expr.coefficient(expr.parent().gens_dict()['gamma'+str(i)]))*cvxpy_vals[i] + float(expr.constant_coefficient()) for i in range(num_bkpt) )
        # data from scip is being passed in
        if scip is not None:
            cutcoefs_expr, cutrhs_expr, integral_indices, lp_soln = self.getParamaterizedCutExpr(scip, cols, rows, binvrow, binvarow, f, sparse_bkpt, f_index)
            cut_generation_problem_logger.debug(f"sage:{cutcoefs_expr, cutrhs_expr}")
            # if the rhs gets used in teh furture; make sure to include the constant terms.
            steepest_dir_expr = sum( map_sage_expr_to_cvpxy_expr(cutcoefs_expr[i])*costs[i] for i in range(len(costs)) )
        else: # assume pure IP; testing and verifcation purposes.
            cut_generation_problem_logger.warning("pure IP problem is assumed.") #TODO: more useful warning.
            pi = pwl_with_value_parameters_and_bkpts_fixed(sparse_bkpt, f_index)
            steepest_dir_expr = sum( map_sage_expr_to_cvpxy_expr(pi(binvarow[i]))*costs[i] for i in range(len(costs)) )
        if self._objective_sense == "maximize":
            obj = cvxpy_Maximize(steepest_dir_expr)
        else:
            obj = cvxpy_Minimize(steepest_dir_expr)
        prob = cvxpy_Problem(obj, minimality_cons)
        if self._cvxpy_args is None and self._cvxpy_kwds is None:
            prob.solve()
        elif self._cvxpy_args is None and self._cvxpy_kwds is not None:
            prob.solve(**self._cvxpy_kwds)
        elif self._cvxpy_args is not None and self._cvxpy_kwds is not None:
            prob.solve(*self._cvxpy_args, **self._cvxpy_kwds)
        else:
            prob.solve(*self._cvxpy_args)
        score = prob.value
        values = to_sage_rationals(cvxpy_vals.value)
        b = sparse_bkpt
        v = values
        pi_p = piecewise_function_from_breakpoints_and_values(b+[1], v+[0])
        log_problem_result(b, v, binvarow, costs, f)
        if self._prove_seperator: # note; in general the minimality test will be false unless using an exact rational solver for the lp problem.
            res = minimality_test(pi_p, self._show_proof) # add someway to log certificates.
            cut_generation_problem_logger.info(f"Minimality of cgf: {res}")
        return pi_p, score, len(sparse_bkpt)

    def getParamaterizedCutExpr(self, scip, cols, rows, binvrow, binvarow, primsol, bkpt, f_index):
        """ Given the row (binvarow, binvrow) of the tableau, computes
        an expression for the cut in the value parameters (gammai).
        
        :param scip: scip model
        :param primsol:  is the rhs of the tableau row
        :param cols:     are the variables
        :param rows:     are the slack variables
        :param binvrow:  components of the tableau row associated to the basis inverse
        :param binvarow: components of the tableau row associated to the basis inverse * A
        :param bkpt:     breakpoint sequence
        :f_index:        for the value polyhedron

        Output
        parametarized cut expression and MIP relaxation lp information for cgp.
        """
        coord_names = ['gamma'+str(i) for i in range(len(bkpt))]
        K = PolynomialRing(QQ, names=coord_names, order='lex')
        # initialize
        cutcoefs_expr = [K.zero()] * len(cols)
        cutrhs_expr = K.zero()
        integral_indices = []
        lp_soln = []
        pi_p = pwl_with_value_parameters_and_bkpts_fixed(bkpt, f_index, poly_ring=K)
        f = K(fractional(QQ(primsol)))
        def psi(x):
            if x < 0:
                slope = pi_p.functions()[-1]._slope
                neg_part = FastLinearFunction(slope, K.zero())
                return neg_part(x)
            else:
                return pi_p.functions()[0](x)
        # rhs of the cut is the fractional part of the LP solution for the basic variable
        cutrhs_expr = -f
        # Generate cut coefficients for the original variables
        for c in range(len(cols)):
            col = cols[c]
            var = col.getVar()
            lp_soln.append(var.getLPSol())
            if var.vtype() != "CONTINUOUS":
                integral_indices.append(c)
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
                cutelem_expr = fractional(QQ(rowelem))
                cutelem_expr = -1*f*pi_p(cutelem_expr)
            else:
                # Continuous variables
                cutelem_expr = -1*f*psi(QQ(rowelem))

            # cut is define when variables are in [0, infty). Translate to general bounds
            if not cutelem_expr is 0: # could be 0 polynomial
                if col.getBasisStatus() == "upper":
                    cutelem_expr = -cutelem_expr
                    cutrhs_expr += cutelem_expr * col.getUb()
                else:
                    cutrhs_expr += cutelem_expr * col.getLb()
                # Add coefficient to cut in dense form
                cutcoefs_expr[col.getLPPos()] = cutelem_expr

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
                cutelem_expr = fractional(QQ(rowelem))
                cutelem_expr = -1*f*pi_p(cutelem_expr)
            else:
                # Continuous variables
                cutelem_expr = -1*f*psi(QQ(rowelem))

            # cut is define in original variables, so we replace slack by its definition
            if cutelem_expr is 0: # sage type; 
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
                    cutelem_expr = -cutelem_expr

                rowcols = row.getCols()
                rowvals = row.getVals()

                assert len(rowcols) == len(rowvals)

                # Eliminate slack variable: rowcols is sorted: [columns in LP, columns not in LP]
                for i in range(row.getNLPNonz()):
                    cutcoefs_expr[rowcols[i].getLPPos()] -= cutelem_expr * rowvals[i]

                act = scip.getRowLPActivity(row)
                rhsslack = rrhs - act
                if scip.isFeasZero(rhsslack):
                    assert row.getBasisStatus() == "upper" # cutelem != 0 and row active at upper bound -> slack at lower, row at upper
                    cutrhs_expr -= cutelem_expr * (rrhs - row.getConstant())
                else:
                    assert scip.isFeasZero(act - rlhs)
                    cutrhs_expr -= cutelem_expr * (rlhs - row.getConstant())
        cut_generation_problem_logger.debug(f"{lp_soln}")
        
        return cutcoefs_expr, cutrhs_expr, integral_indices, lp_soln
        # apply algebraic maps to produce cuts in terms of val; val is non-> return list
#        n = len(bkpt)
#        coord_names = ['gamma'+str(i) for i in range(n)]
#        map_to_vector = lambda expr : [expr.coefficient(expr.parent().gens_dict()[name]) for name in coord_names]
#        list_of_vecs_in_gamma_params = [ map_to_vector(expr) for expr in cutcoefs_expr ]
#        cutrhs_expr [ sum( vec[i]*val[i] for i in range(n)) for vec in list_of_vecs_in_gamma_params ]
#

    def _algorithm_custom(self, binvarow, binvc, f):
        # For future use.
        raise NotImplementedError

    
    def get_cgp_input_parameters(self):
        """
        Returns a dictionary parameters used to initalize the cut generation problem.
        """
        return self._cgp_input_parameters
        
    
