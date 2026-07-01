from cutgeneratingfunctionology.igp import *
from .cgf_specializations import *
from scipy.optimize import minimize, LinearConstraint, NonlinearConstraint
from .execptions import *
import logging

scipy_interface_logger = logging.getLogger(__name__)
scipy_interface_logger.setLevel(logging.ERROR)


def map_polyhedral_bsa_to_scipy_LinearConstraint(bsa, epsilon=10**-9):
    """
    Maps a cutgeneratingfunctionlogy ``basic_semialgebraic_set`` polynomial constraints to scipy constraint objects. 
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
    Given a BSA with nonlinear constraints, converts into an equivalent set of nonlinear constraints for scipy.

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

class scipyCutScoreForBkptAsParam:
    def __init__(self, cut_score=None, epsilon=1e-8, M=1e8):
        """
        cut_score: name or function with signature (cut_lhs, cut_rhs, binvc, lp_solution). 
        """
        if cut_score is None or cut_score == "scip":
            self.cut_score = scipy_scip_standard_cut_score
            # TODO: self.grad_cut_score = 
        elif cut_score == "parallelism":
            self.cut_score = scipy_objective_function_parallelism
        elif cut_score == "violation":
            self.cut_score = scipy_violation
        elif cut_score == "realitive_violation":
            self.cut_score = scipy_realitive_violation
        elif cut_score == "cut_off_distance":
            self.cut_score = scipy_cut_off_distance
        else:
            self.cut_score = cut_score
        self._epsilon = epsilon
        self._current_point = None

    def gen_cut_score(self, f_index,  cutcoefs_expr, cutrhs_expr, current_feasible_soln, binvc, integral_indices):
        """
        Defines cut score based on fixed problem data.
        """
        def s(value_parameters):
            val =  to_sage_rationals(value_parameters)
            # val = self.validate_point(val, f_index)
            # self._current_point = val
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

    def validate_point(self, value_parameters, f_index):
        # This enforces the "manifoldness" of PWL.
        # we just need to ensure that our values output are for sure correct. 
        epsilon = self._epsilon
        n = len(value_parameters)
        v = to_sage_rationals(value_parameters)
        # values should be in [0,1]
        for i in range(n):
            if (v[i] < 0 or v[i] + epsilon> 1) and i != f_index:
                raise ModelViolation(f"validate_point: breakpoint gamma_{i} < 0 or gamma_{i}>= 1; gamma_{i}=={v[i]}")
        # pi(0) = 0, pi(f) = 1
        if abs(v[0]-0) <= epsilon:
            v[0] = 0
        else:
            raise ModelViolation("value gamma_0 >0")
        # pi_p(f) = 1
        if abs(v[f_index] - 1) <= epsilon:
            v[f_index] = 1
        else:
            scipy_interface_logger.debug(f"validate_point: breakpoint lambda_{f_index} !=1: point: {point}, cell: {cell}")
            raise ModelViolation(f"value gamma_{f_index} != 1")
        return v
        # We are assuming breakpoints are model multiplicity free.
        # sprace enough breakpoints ensure no two breakpoitns are close enough to need an epsilon chart when solving.

def scipy_objective_function_parallelism(cut_lhs, cut_rhs, binvc, lp_solution, integral_indices):
    # eval cutcoefs_expr at val to get cut. 
    cut = np.array(cut_lhs)/np.linalg.norm(cut_lhs)
    c = np.array(binvc)/np.linalg.norm(binvc)
    # scipy_interface_logger.debug(f"{cut @ c}")
    return cut @ c

def grad_scipy_objective_function_parallelism(cut_lhs, cut_rhs, binvc, lp_solution, integral_indices):
    raise NotImplementedError

#def scipy_scip_standard_cut_score(cut_lhs, cut_rhs, binvc, lp_solution, integral_indices):
#    pass

def scipy_violation(cut_lhs, cut_rhs, binvc, lp_solution, integral_indices):
    scipy_interface_logger.debug(f"{np.array(cut_lhs)@ np.array(lp_solution) - cut_rhs}")
    return np.array(cut_lhs)@ np.array(lp_solution) - cut_rhs
    
def scipy_realitive_violation(cut_lhs, cut_rhs, binvc, lp_solution, integral_indices):
    if cut_rhs != 0:
        return scipy_violation(cut_lhs, cut_rhs, binvc, lp_solution, integral_indices)/abs(cut_rhs)
    else:
        return scipy_violation(cut_lhs, cut_rhs, binvc, lp_solution, integral_indices)

def scipy_cut_off_distance(cut_lhs, cut_rhs, binvc, lp_solution, integral_indices):
    return scipy_violation(cut_lhs, cut_rhs, binvc, lp_solution, integral_indices)/np.linalg.norm(np.array(cut_lhs))

def integer_support(cut_lhs, cut_rhs, binvc, lp_solution, integral_indices):
    support_count = 0
    for c in integral_indices:
        if abs(cut_lhs[c]) > 1e-9:
            support_count += 1
    return support_count
    
