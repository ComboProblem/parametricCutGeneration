from cutgeneratingfunctionology.igp import *
from .cgf_specializations import *
from scipy.optimize import minimize, LinearConstraint, NonlinearConstraint
from .execptions import *
import logging

r"""
Inteface from :scipy: non-linear solvers interface to a :class:``cutGenerationProblem``.
"""

scipy_interface_logger = logging.getLogger(__name__)
scipy_interface_logger.setLevel(logging.ERROR)


def map_polyhedral_bsa_to_scipy_LinearConstraint(bsa, epsilon=10**-9):
    """
    Maps objects from :cutgeneratingfunctionlogy::spam::basic_semialgebraic:
    with linear constraints to :scipy: constraints as matrices :np.array:
    with shape (n,m) where n is the number of variables and m is the
    number of polynoimal constraints of the `bsa`. 
    """
    lb = []
    ub = []
    eq = []
    A_ieq = []
    A_eq = []
    for polynomial in bsa.eq_poly():
        if polynomial.degree() != 1:
            raise ValueError(f"Constraint {polynomial} == 0 is not linear.")
        coeffs = [polynomial.coefficient(i) for i in polynomial.parent().gens()]+[polynomial.constant_coefficient()]
        normalization_term = max([abs(int(x)) for x in coeffs])
        if normalization_term != 0:
            coeffs = np.array(coeffs)/normalization_term
        # polys in bsa are written poly op 0
        # rewrite to correct scipy notation.
        A_eq.append(coeffs[:-1])
        eq.append( -1*coeffs[-1])
    for polynomial in bsa.lt_poly():
        if polynomial.degree() != 1:
            raise ValueError(f"Constraint {polynomial} < 0 is not linear.")
        coeffs = [polynomial.coefficient(i) for i in polynomial.parent().gens()]+[polynomial.constant_coefficient()]
        normalization_term = max([abs(int(x)) for x in coeffs])
        if normalization_term != 0:
            coeffs = np.array(coeffs)/normalization_term
        # polys in bsa are written poly op 0
        # rewrite to correct scipy notation.
        A_ieq.append(coeffs[:-1])
        lb.append(-np.inf)
        ub.append(-1*coeffs[-1])
    for polynomial in bsa.le_poly():
        if polynomial.degree() != 1:
            raise ValueError(f"Constraint {polynomial} < 0 is not linear.")
        coeffs = [polynomial.coefficient(i) for i in polynomial.parent().gens()]+[polynomial.constant_coefficient()]
        normalization_term = max([abs(int(x)) for x in coeffs])
        if normalization_term != 0:
            coeffs = np.array(coeffs)/normalization_term
        # polys in bsa are written poly op 0
        # rewrite to correct scipy notation.
        A_ieq.append(coeffs[:-1])
        lb.append(-np.inf)
        ub.append(-1*coeffs[-1])
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

def map_bsa_to_scipy_NonLinearConstraint(bsa, epsilon=10**-9):
    r"""
    Maps objects from  :cutgeneratingfunctionlogy::spam::basic_semialgebraic: with polynomial constraints
    to :scipy:``NonLinearConstraint`` constraints by defining polynomial functions on np.array with shape
    (n,) with output of np.array with shap (m,) where n is the number of varialbes of the bsa and m is   
    the number of constraints. Graidents are computed algebraically and provided to :scipy:. 
    
    To model < as <=, we treat c_i < 0 if and only if c_i + epislon <= 0 since :scipy: does not support 
    strict inequalities in modeling. 
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

class scipyCutScoreForBkptAsParam:
    """
    For defining paramaterized objective functions for
    :class:`cutGenenratingProblem`when the `algorithm`
    selected is \'bkpt_as_param\'.

    Options include...
    \'parallelism\' - for objective fucntion parallelism.
    \'violation\' - to optimize for violation.
    \'realitive_violation\' - to optimize for realitve violation.
    \'cut_off_distance\' - to optimize for cut off distance 
    """
    def __init__(self, cut_score=None):
        """
        :cut_score: name or function with signature (cut_lhs, cut_rhs, binvc, lp_solution). 
        """
        if cut_score is None or cut_score == "scip":
            self.cut_score = self.scipy_scip_standard_cut_score
            # TODO: self.grad_cut_score = 
        elif cut_score == "parallelism":
            self.cut_score = self.scipy_objective_function_parallelism
        elif cut_score == "violation":
            self.cut_score = self.scipy_violation
        elif cut_score == "realitive_violation":
            self.cut_score = self.scipy_realitive_violation
        elif cut_score == "cut_off_distance":
            self.cut_score = self.scipy_cut_off_distance
        else:
            self.cut_score = cut_score
        self._current_point = None

    def gen_cut_score(self, f_index,  cutcoefs_expr, cutrhs_expr, current_feasible_soln, binvc, integral_indices):
        """
        Defines paramaterized objective function of  score based on fixed problem data.

        Input:
        :f_index: - 
        :cutcoefs_expr: - list of polynomial expression defined on (b,v) for the paramaterized cut's coeffients (lhs) based on the current node, row selecte, ect. 
        :cutrhs_expr: - polynomial experssions defined on (b,v) for the  cut's rhs in for the given node.
        :current_feasible_soln: - LP solution of the current B&B node.
        :binvc: - the reduced costs of the LP of teh current B&B node.
        """
        def s(value_parameters):
            val =  to_sage_rationals(value_parameters)
            coord_names = ['gamma'+str(i) for i in range(len(value_parameters))] # same as num bkpts
            map_to_vector = lambda expr : [expr.coefficient(expr.parent().gens_dict()[name]) for name in coord_names]
            vectors = [ map_to_vector(expr) for expr in cutcoefs_expr ] # inner lists are coeffics of gammis, outer lists are for orignal coordinates
            consts = [float(expr.constant_coefficient()) for expr in cutcoefs_expr ]
            #scipy_interface_logger.debug(f"val:{val} cutcoeffs_expr:{cutcoefs_expr} vectors:{vectors} consts:{consts}")
            cut_lhs = [float(sum(v * u for v, u in zip(val,vec))) + cst for vec, cst in zip(vectors, consts)]
            cut_rhs = float(sum(v*u for v, u in zip(val, map_to_vector(cutrhs_expr))))+ float(cutrhs_expr.constant_coefficient())
            #scipy_interface_logger.debug(f"{cut_lhs, cut_rhs}")
            return self.cut_score(cut_lhs, cut_rhs, binvc, current_feasible_soln, integral_indices)
        return s

    def scipy_objective_function_parallelism(self, cut_lhs, cut_rhs, binvc, lp_solution, integral_indices):
        """
        Paramaterized version of parallelism cut score.
        
        See `scipyCutScoreForBkptAsParam.gen_cut_score` for parameter arguments.
        """
        # eval cutcoefs_expr at val to get cut. 
        cut = np.array(cut_lhs)/np.linalg.norm(cut_lhs)
        c = np.array(binvc)/np.linalg.norm(binvc)
        # scipy_interface_logger.debug(f"{cut @ c}")
        return cut @ c

    def scipy_violation(self, cut_lhs, cut_rhs, binvc, lp_solution, integral_indices):
        """
        Paramaterized version of violation cut score.
        
        See `scipyCutScoreForBkptAsParam.gen_cut_score` for parameter arguments.
        """
        scipy_interface_logger.debug(f"{np.array(cut_lhs)@ np.array(lp_solution) - cut_rhs}")
        return np.array(cut_lhs)@ np.array(lp_solution) - cut_rhs
        
    def scipy_realitive_violation(self, cut_lhs, cut_rhs, binvc, lp_solution, integral_indices):
        """
        Paramaterized version of realitive violation cut score.
        
        See `scipyCutScoreForBkptAsParam.gen_cut_score` for parameter arguments.
        """
        if cut_rhs != 0:
            return scipy_violation(cut_lhs, cut_rhs, binvc, lp_solution, integral_indices)/abs(cut_rhs)
        else:
            return scipy_violation(cut_lhs, cut_rhs, binvc, lp_solution, integral_indices)

    def scipy_cut_off_distance(self, cut_lhs, cut_rhs, binvc, lp_solution, integral_indices):
        """
        Paramaterized version of cut off distance cut score.
        
        See `scipyCutScoreForBkptAsParam.gen_cut_score` for parameter arguments.
        """
        return scipy_violation(cut_lhs, cut_rhs, binvc, lp_solution, integral_indices)/np.linalg.norm(np.array(cut_lhs))

    def integer_support(self, cut_lhs, cut_rhs, binvc, lp_solution, integral_indices):
        """
        Paramaterized version of integer support cut score.
        
        See `scipyCutScoreForBkptAsParam.gen_cut_score` for parameter arguments.
        """
        raise NotImplementedError
        support_count = 0
        for c in integral_indices:
            if abs(cut_lhs[c]) > 1e-9:
                support_count += 1
        return support_count
        


