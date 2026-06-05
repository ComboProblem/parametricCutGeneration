from cutgeneratingfunctionology.igp import *

# specialized code for solving value poly and breakpoint as parameters cgps.
# This code is modified from cut generating functionology.

def to_sage_rationals(structure):
    return [QQ(i) for i in structure]

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


def generate_type_1_vertices_continuous_expr(fn, bkpt=None):
    r"""Output 6-tuples (x, y, z,xeps, yeps, zeps).
    """
    if bkpt is None:
        bkpt = fn.end_points()
    return ( delta_pi(fn,x,y) for x in bkpt for y in bkpt if x <= y ) # generator comprehension

def generate_type_2_vertices_continuous_expr(fn, bkpt=None):
    if bkpt is None:
        bkpt = fn.end_points()
    bkpt2 = bkpt[:-1] + [ x+1 for x in bkpt ]
    return ( delta_pi(fn, x, z-x) for x in bkpt for z in bkpt2 if x < z < 1+x ) # generator comprehension

def generate_assumed_symmetric_vertices_continuous_expr(fn, f, bkpt):
    """
    Silently assumes the symmetry condition holds for all vertices (x,y) in bkpt's breakpoints complex
    such that x+y equiv f.

    See function::``generate_symmetric_vertices_continuous``.
    """
    for i in range(len(bkpt)):
        x = bkpt[i]
        if x == f:
            continue
        if x < f:
            y = f - x
        else:
            y = 1 + f - x
        yield fn(x) + fn(y) 

def value_nnc_polyhedron_constraints(bkpt, f_index, val=None, coeff_type='rational', espilon=.00000000001, log_pw_functions=False):
    r"""
    Returns a list of constraints, possibiely evaluated, representing the value polyhedron for the given (proper) breakpoint sequence, f index, and values.  
    INPUTS::
    -`bkpt` - a proper breakpoint sequence, list like of sage rationals.
    -`val` - values, of the same legth as bkpt, list like.
    -`f_index` - pi(bkpt[f_index]) == 1 
    -`coeff_type` - python int, python float, or (sage) rational (default).
    -`log_pw_functions` - verboicity of logging peicewise function via the logging module.

    OUPUTS::
    - `cons` - a list of constraints repersenting the value polyhedron.
    EXAMPLES::
    >>> bkpt = to_sage_rationals([0, 4/5])
    >>> val = to_sage_rationals([0, 1])
    >>> all(value_nnc_polyhedron_constraint_like(bkpt, 1, val)) # equivlant to fininte minimality test on gmic
    True
    >>> import ppl
    >>> ppl_vars = [ppl.Variable(i) for i in range(2)]
    >>> cons = value_nnc_polyhedron_constraint_like(bkpt, ppl_vars, 1, 'int')
    >>> ppl_V =  ppl.Polyhedron(2)
    >>> for con in cons:
            ppl_V.add_constraint(con)
    >>> ppl_V
    """
    n = len(bkpt)
    if not log_pw_functions:
        pw_logging_level = logging.getLogger("cutgeneratingfunctionology.igp.functions").getEffectiveLevel()
        logging.getLogger("cutgeneratingfunctionology.igp.functions").setLevel(logging.ERROR)
    assert(n >= 2)
    assert(f_index >= 1)
    assert(f_index <= n - 1)
    if not isinstance(bkpt, list):
        bkpt = list(bkpt)
    coord_names = ['gamma'+str(i) for i in range(n)]
    if val is None:
        K = PolynomialRing(QQ, names=coord_names, order='lex')
        val = K.gens() #suitable rational ring
    h = pwl_with_value_parameters_and_bkpts_fixed(bkpt, f_index)
    cons = []
    cons.append(val[0] == 0)
    cons.append(val[f_index] == 1)
    for i in range(1,n):
        try:
            cons.append( val[i] > 0 )
        except Exception as e:
            cons.append( val[i] + espilon >= 0)
        cons.append( val[i] <= 1 )
    # Assumes minimality for the partially defined function.
    for expr in generate_type_1_vertices_continuous_expr(h):
        try:
            if len(expr.coefficients()) == 0:
                continue
        except Exception:
            continue
        try:
            lhs = [expr.coefficient(expr.parent().gens_dict()[name]) for name in coord_names]
            cst = expr.constant_coefficient()
            if coeff_type == "int":
                lcd = lcm(lcm([coeff.denominator() for coeff in lhs]), cst.denominator())
                cons.append(sum([int(lcd * lhs[i]) * val[i] for i in range(n)]) + lcd * cst >= 0)
            elif coeff_type == "float":
                cons.append(sum([float(lhs[i]) * val[i] for i in range(n)]) + float(cst) >= 0)
            else: #sage rational type assumed
                cons.append(sum([float(lhs[i]) * val[i] for i in range(n)]) + float(cst) >= 0)
        except Exception as e:
            if expr is 0:
                continue
            else:
                raise e
    for expr in generate_type_2_vertices_continuous_expr(h):
        try:
            if len(expr.coefficients()) == 0:
                continue
        except Exception:
            continue
        try:
            lhs = [expr.coefficient(expr.parent().gens_dict()[name]) for name in coord_names]
            cst = expr.constant_coefficient()
            if coeff_type == "int":
                lcd = lcm(lcm([coeff.denominator() for coeff in lhs]),  cst.denominator())
                cons.append(sum([int(lcd * lhs[i]) * val[i] for i in range(n)]) + lcd * cst >= 0)
            elif coeff_type == "float":
                cons.append(sum([float(lhs[i]) * val[i] for i in range(n)]) + float(cst) >= 0)
            else: #sage rational type assumed
                cons.append(sum( [float(lhs[i]) * val[i] for i in range(n)]) + float(cst) >= 0)
        except Exception as e:
            continue
    for expr in generate_assumed_symmetric_vertices_continuous_expr(h, bkpt[f_index], bkpt + [1]):
        try:
            if len(expr.coefficients()) == 0:
                continue
        except Exception:
            continue
        try:
            lhs = [expr.coefficient(expr.parent().gens_dict()[name]) for name in coord_names]
            cst = expr.constant_coefficient()
            if coeff_type == "int":
                lcd = lcm(lcm([coeff.denominator() for coeff in lhs]), cst.denominator())
                cons.append(sum([int(lcd * lhs[i]) * val[i] for i in range(n)]) + lcd * cst == 1*lcd)
            elif coeff_type == "float":
                cons.append(sum([float(lhs[i]) * val[i] for i in range(n)]) + float(cst) == 1)
            else: #sage rational type assumed
                cons.append(sum([lhs[i] * val[i] for i in range(n)]) + float(cst) == 1)
        except Exception as e:
            if expr is 1:
                continue
            else:
                raise e
    if not log_pw_functions:
        logging.getLogger("cutgeneratingfunctionology.igp.functions").setLevel(pw_logging_level)
    return cons


def expression_of_steepest_direction_score(binvc, bkpt, f_index, val):
    pi = pwl_with_value_parameters_and_bkpts_fixed(bkpt, f_index)
    cut_score_in_value_params = sum(pi(fractional(QQ(c))) for c in binvc)
    expr = sum([cut_score_in_value_params.coefficient(cut_score_in_value_params.parent().gens_dict()['gamma'+str(i)])*val[i] for i in range(len(bkpt))])
    return expr


