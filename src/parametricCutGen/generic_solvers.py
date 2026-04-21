from cutgeneratingfunctionology.igp import *
from scipy.optimize import minimize, LinearConstraint, NonlinearConstraint
from cvxpy import Variable, Minimize, Problem
from .execptions import *
import logging

generic_solver_logger = logging.getLogger(__name__)

# TODO: Make base classes actually base classes rather than super classing things. This applies for cut_score module too.

class abstractCutGenProblemSolverInterface:
    r"""
    Interfaces types from ``cutgeratingfunctionolgy`` to a specified generic solver which solves the real optimization problem min f(x) s.t. g(x) <= 0; x\in R^n.
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

        constraints and objective should be of the solvers type.

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
        Return the correct solver type from a sage ring element.
        """
        raise NotImplementedError

# TODO: typing of BSA; or write a version which is typed. 

class scipyCutGenProbelmSolverInterface(abstractCutGenProblemSolverInterface):
    """
    Maps types and objects from ``cutgeneratingfunctiology`` to scipy.
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
            # mimic in a non-rigorous way <
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
        Given converted constraints and an objective function that is compatible with the solver,
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
        scipy supports inputs of floats, convert sage field element to its equivalent numerical (python) floating point value.
        """
        # be lazy and assume the_ring_element is something the converts to a rational number (or can be put into a floating point approximation).
        # this is a point where we lose the exactness of sage.
        return float(sage_ring_element)


class cvxpyCutGenProblemSolverInterface(abstractCutGenProblemSolverInterface):
    r"""
    Maps types and objects from ``cutgeneratingfunctiology`` to cvxpy.
    """
    @staticmethod
    def write_linear_constraints_from_bsa(bsa):
        r"""
        Given a BSA with only linear constraints, converts the bsa object into a format that the underlying solver can use.
        """
    def write_linear_constraints_from_bsa(bsa, epsilon=10**-9):
        r"""
        Given a BSA with only linear constraints, converts the bsa object into a format that the underlying solver can use.
        """
        # TODO:linear specialized code for linear BSAs?
        x = cp.Variable(bsa.ambient_dim())
        cons = []
        for polynomial in bsa.eq_poly():
            if polynomial.degree() != 1:
                raise ValueError(f"Constraint {polynomial} == 0 is not linear.")
            linear_coeffs = np.array([polynomial.coefficient(i) for i in polynomial.parent().gens()])
            cons.append(linear_coeffs @ x == float(-1*polynomial.constant_coefficient()))
        for polynomial in bsa.lt_poly():
            if polynomial.degree() != 1:
                raise ValueError(f"Constraint {polynomial} < 0 is not linear.")
            linear_coeffs = np.array([polynomial.coefficient(i) for i in polynomial.parent().gens()])
            cons.append(linear_coeffs @ x <= float(-1*polynomial.constant_coefficient() - epsilon))
        for polynomial in bsa.le_poly():
            if polynomial.degree() != 1:
                raise ValueError(f"Constraint {polynomial} <= 0 is not linear.")
            linear_coeffs = np.array([polynomial.coefficient(i) for i in polynomial.parent().gens()])
            cons.append(linear_coeffs @ x <= float(-1*polynomial.constant_coefficient()))
        if len(cons) != 0:
            return cons, x
        else:
            raise ValueError(f"No constraints have been written.")        

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

        Should return optimal objective value, optimal objective solution, solver success, and solver_output.
        """
        x = solver_options['x']# intnded to be form the constraints x
        prob = Problem(objective, constraints)
        result = prob.solve()        
        return prob.value, x.value, None, prob

    @staticmethod
    def nonlinear_solve(constraints, objective, **solver_options):
        r"""
        Interface to solver's min/max f(x) s.t. p_i(x) <= b_i, where p_i is a polynomial and at least 1 p_i has degree larger than 1.
        """
        raise NotImplementedError


    @staticmethod
    def sage_to_solver_type(sage_ring_element):
        r"""
        Return the correct solver type from a sage ring element.
        """
        raise NotImplementedError
