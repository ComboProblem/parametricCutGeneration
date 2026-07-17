"""
Specalized code based on cutgeneratingfunctionology for solving value poly lp and breakpoints as parameters cut generation problems. 
"""

from cutgeneratingfunctionology.igp import *
from pplite import Variable as pplite_Variable, Constraint as pplite_Con, Linear_Expression as pplite_Lin_Expr, Affine_Expression as pplite_Aff_expr, NNC_Polyhedron as pplite_NNC_Polyhedron, PPliteGenerator

cgf_special_logger = logging.getLogger(__name__)
cgf_special_logger.setLevel(logging.ERROR)


def to_sage_rationals(iterable):
    """Helper function to converte to sage rationals when needed."""
    return [QQ(i) for i in iterable]

def cgp_piecewise_function_from_breakpoints_slopes_and_values(bkpt, slopes, values, merge=True):
    r"""
    Create a continuous piecewise function from bkpt, slopes, and values.

    - bkpt and values are two parallel lists; it is assumed that bkpt is sorted in increasing order. 

    - slopes is one element shorter and represents the slopes of the interpolation.

    - The function is overdetermined by these data.  The consistency of the data is currently not checked.

    - The data are coerced into a common convenient field via ``nice_field_values``.

    - If merge is ``True`` (the default), adjacent pieces of equal slopes are merged into one.
    """
    # global symb_values
    if slopes is None: # make order assumptions in these functions, remove != when possible.
        slopes = [(values[i+1]-values[i])/(bkpt[i+1]-bkpt[i]) if bkpt[i+1] > bkpt[i] else 0 for i in range(len(bkpt)-1)] 
    else:
        symb_values = bkpt + slopes + values
        field_values = nice_field_values(symb_values, field)
        bkpt, slopes, values = field_values[0:len(bkpt)], field_values[len(bkpt):len(bkpt)+len(slopes)], field_values[-len(values):]
    intercepts = [ values[i] - slopes[i]*bkpt[i] for i in range(len(slopes)) ]
    # Make numbers nice
    ## slopes = [ canonicalize_number(slope) for slope in slopes ]
    ## intercepts = [ canonicalize_number(intercept) for intercept in intercepts ]
    #print slopes
    pieces = []
    for i in range(len(bkpt)-1):
        if bkpt[i] > bkpt[i+1]:
            raise ValueError("Breakpoints are not sorted in increasing order.")
        elif bkpt[i] == bkpt[i+1]:  #hasattr(field, '_big_cells') and field._big_cells, should be off-record
            functions_logger.warning("Degenerate interval occurs at breakpoint %s" % bkpt[i])
            if values[i] != values[i+1]:
                raise ValueError("Degeneration leads to a discontinuous function.")
        else:
            pieces.append( [(bkpt[i],bkpt[i+1]),
                            fast_linear_function(slopes[i], intercepts[i])] )
    return FastPiecewise(pieces, merge=merge)


def pwl_with_value_parameters_and_bkpts_fixed(bkpt, f_index, base_ring=QQ, poly_ring=None, log_pw_functions=False):
    """
    Piecewise linear function by breakpoints and values where the breakpoints are assumed to be fixed and the values of
    the function at each breakpoint are monomials gammai.

    Evaluating this function numerically returns a polynomial expressions (in fact linear) in the polyomial ring over the base ring
    with variables gammai for i in 0 ... len(bkpt) - 1.
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
    if poly_ring is None:
        coord_names = ['gamma'+str(i) for i in range(n)]
        K = PolynomialRing(base_ring, names=coord_names, order='lex')
    else:
        K = poly_ring
    one = K.one() # zero and one should be elements of the polynomial ring
    zero = K.zero()
    bkpt_K = [K(b) for b in bkpt]
    vals = [zero] + [K.gens()[i] if i != f_index else one for i in range(1,n)]
    if not log_pw_functions:
        logging.getLogger("cutgeneratingfunctionology.igp.functions").setLevel(pw_logging_level) 
    return cgp_piecewise_function_from_breakpoints_slopes_and_values(bkpt_K + [one], None,  vals + [zero], merge=False) # cgf bkpt val param has some unwanted coheasion.


def generate_type_1_vertices_continuous_expr(fn):
    r"""Output expressions for delta pi at type 1 vertices.
    """
    bkpt = fn.end_points()
    cgf_special_logger.debug(f"{type(bkpt[0])}")
    return ( delta_pi(fn,x,y) for x in bkpt for y in bkpt if x <= y ) # generator comprehension

def generate_type_2_vertices_continuous_expr(fn):
    """
    Output expressions for delta pi at type 1 vertices.
    """
    bkpt = fn.end_points()
    bkpt2 = bkpt[:-1] + [ x+1 for x in bkpt ]

    return ( delta_pi(fn, x, z-x) for x in bkpt for z in bkpt2 if x < z < 1+x ) # generator comprehension

def generate_symmetric_vertices_continuous_expr(fn, f):
    """
    Silently assumes the symmetry condition holds for all vertices (x,y) in bkpt's breakpoints complex
    such that x+y equiv f.

    See function::``generate_symmetric_vertices_continuous``.
    """
    bkpt = fn.end_points()
    for i in range(len(bkpt)):
        x = bkpt[i]
        if x == f:
            continue
        if x < f:
            y = f - x
        else:
            y = 1 + f - x
        yield fn(x) + fn(y) 

def value_nnc_polyhedron_constraints(bkpt, f_index, values=None, *, coeff_type='int', epsilon=.00000000001, log_pw_functions=False, backend=None, feasiblity_problem=False, ring=QQ):
    r"""
    Returns a list of constraints, possibiely evaluated, representing the value polyhedron for the given (proper) breakpoint sequence, f index, and values.

    Positional INPUTS::
    -`bkpt` - a proper breakpoint sequence, list like of sage rationals.
    -`f_index` - pi(bkpt[f_index]) == 1
    Keyword INPUTS::
    -`val` - values, of the same legth as bkpt, list like - default list of ppl.linear_algebra.Varible
    -`coeff_type` - string, 'int' (python.int),  'float' (python.float), or 'rational' (sage.rings.rational.Rational)
    -`log_pw_functions` - verboicity of logging peicewise function via the logging module.

    OUPUTS::
    - `cons` - a list of constraints repersenting the value polyhedron.
    
    Reproduce finite minimality test::
    >>> bkpt = to_sage_rationals([0, 4/5]) 
    >>> val = to_sage_rationals([0, 1])
    >>> all(value_nnc_polyhedron_constraints(bkpt, 1, val)) # equivlant to fininte minimality test on gmic with f=4/5 
    True
    
    Examples::
    >>> cons = value_nnc_polyhedron_constraints(bkpt, 1)
    >>> Val_Poly = NNC_Polyhedron(2)
    >>> for con in cons:
            Val_Poly.add_constraint(con)
    >>> Val_Poly
    A 0-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point
    >>> Val_Poly.generators() # The result is expected. The value for pi(4/5=bkpt[1]) == 1 as in gmic. 
    Generator_System {point(0/1, 1/1)}
    >>> from cvxpy import Maximize, Problem, cvxpy_Variable
    >>> cvxpy_val = cvxpy_Variable(4)
    >>> binvarow = [3.2, 4.1, 5.6, .2]
    >>> binvc = [1.2, 4.4, 5.6, -.1]
    >>> frac_bkpts = [fractional(i) for i in to_sage_rationals(binvarow)]
    >>> objective = Maximize(sum(binvc[i]*cvxpy_val[i] for i in range(4))) # paralleism cut score
    >>> val_poly_cons = value_nnc_polyhedron_constraints(frac_bkpts, 3, cvxpy_val, coeff_type='float') # it is reccomended to pass coeff_type='float' if using in conjuction with a solver.
    >>> prob = Problem(objective, val_poly_cons)
    >>> prob.solve()
    >>> cvxpy_val.values
    """
    n = len(bkpt)
    if not log_pw_functions:
        pw_logging_level = logging.getLogger("cutgeneratingfunctionology.igp.functions").getEffectiveLevel()
        logging.getLogger("cutgeneratingfunctionology.igp.functions").setLevel(logging.ERROR)
    assert(n >= 2)
    assert(f_index >= 1)
    assert(f_index <= n - 1)
    coord_names = ['gamma'+str(i) for i in range(n)]
    if values is None:
        if backend == 'pplite':
            Var = pplite_Variable
            coeff_type = 'int'
            values = [Var(i) for i in range(n)]
        elif backend == 'sage_mip':
            if backend_args["sage_mip"]["mip"] is not None:
                backend_args["sage_mip"]["mip"] = mip
                values = mip.default_variable()
                coeff_type = 'ring'
            else:
                raise ValueError("please provide a MixedIntegerProgram")
        else:
            Var = Variable
            coeff_type = 'int'
            values = [Var(i) for i in range(n)]
    if ring is not None:
        R=ring
    else:
        R=QQ
    h = pwl_with_value_parameters_and_bkpts_fixed(bkpt, f_index)
    cons = []
    cons.append( values[0] == 0 )
    cons.append( values[f_index] == 1 )
    for i in range(1,n):
        try:
            cons.append( values[i] > 0 )
        except Exception as e:
            if coeff_type == "int":
                d = QQ(epsilon).denominator()
                cons.append( int(d)*values[i] - int(d * epsilon) >= 0)
            elif coeff_type == "float":
                cons.append( values[i] - float(epsilon) >= 0 )
            else:
                cons.append( values[i] - float(epsilon) >= 0 )
        cons.append( values[i] <= 1 )
    # Assumes minimality for the partially defined function.
    if not feasiblity_problem:
        for linear_poly_expr in generate_type_1_vertices_continuous_expr(h):
            cgf_special_logger.debug(f"expr:{linear_poly_expr}, type:{type(linear_poly_expr)}")
            lhs = [linear_poly_expr.coefficient(linear_poly_expr.parent().gens_dict()[name]) for name in coord_names]
            cst = linear_poly_expr.constant_coefficient()
            if coeff_type == "int":
                lcd = lcm(lcm([coeff.denominator() for coeff in lhs]), cst.denominator())
                cons.append(sum([int(lcd * lhs[i]) * values[i] for i in range(n)]) + int(lcd * cst) >= 0)
            elif coeff_type == "float":
                cons.append(sum([float(lhs[i]) * values[i] for i in range(n)]) + float(cst) >= 0)
            else: #sage rational type assumed
                cons.append(sum([lhs[i] * values[i] for i in range(n)]) + float(cst) >= 0)
        for linear_poly_expr in generate_type_2_vertices_continuous_expr(h):
            lhs = [linear_poly_expr.coefficient(linear_poly_expr.parent().gens_dict()[name]) for name in coord_names]
            cst = linear_poly_expr.constant_coefficient()
            if coeff_type == "int":
                lcd = lcm(lcm([coeff.denominator() for coeff in lhs]),  cst.denominator())
                cons.append(sum([int(lcd * lhs[i]) * values[i] for i in range(n)]) + int(lcd * cst) >= 0)
            elif coeff_type == "float":
                cons.append(sum([float(lhs[i]) * values[i] for i in range(n)]) + float(cst) >= 0)
            else: #sage rational type assumed
                cons.append(sum( [lhs[i] * values[i] for i in range(n)]) + float(cst) >= 0)
        for linear_poly_expr in generate_symmetric_vertices_continuous_expr(h, bkpt[f_index]):
            lhs = [linear_poly_expr.coefficient(linear_poly_expr.parent().gens_dict()[name]) for name in coord_names]
            cst = linear_poly_expr.constant_coefficient()
            if coeff_type == "int":
                lcd = lcm(lcm([coeff.denominator() for coeff in lhs]), cst.denominator())
                cons.append(sum([int(lcd * lhs[i]) * values[i] for i in range(n)]) + int(lcd * cst) == 1*lcd)
            elif coeff_type == "float":
                cons.append(sum([float(lhs[i]) * values[i] for i in range(n)]) + float(cst) == 1)
            else: #sage rational type assumed
                cons.append(sum([lhs[i] * values[i] for i in range(n)]) + float(cst) == 1)
    else:
        # assert len(values) == n+1
        cons.append(values[n] >= 0)
        for linear_poly_expr in generate_type_1_vertices_continuous_expr(h):
            cgf_special_logger.debug(f"expr:{linear_poly_expr}, type:{type(linear_poly_expr)}")
            lhs = [linear_poly_expr.coefficient(linear_poly_expr.parent().gens_dict()[name]) for name in coord_names]
            cst = linear_poly_expr.constant_coefficient()
            if coeff_type == "int":
                lcd = lcm(lcm([coeff.denominator() for coeff in lhs]), cst.denominator())
                cons.append(sum([int(lcd * lhs[i]) * values[i] for i in range(n)]) + int(lcd * cst) + values[n] >= 0)
            elif coeff_type == "float":
                cons.append(sum([float(lhs[i]) * values[i] for i in range(n)]) + float(cst) + values[n] >= 0)
            else: #sage rational type assumed
                cons.append(sum([lhs[i] * values[i] for i in range(n)]) + float(cst) + values[n] >= 0)
        for linear_poly_expr in generate_type_2_vertices_continuous_expr(h):
            lhs = [linear_poly_expr.coefficient(linear_poly_expr.parent().gens_dict()[name]) for name in coord_names]
            cst = linear_poly_expr.constant_coefficient()
            if coeff_type == "int":
                lcd = lcm(lcm([coeff.denominator() for coeff in lhs]),  cst.denominator())
                cons.append(sum([int(lcd * lhs[i]) * values[i] for i in range(n)]) + int(lcd * cst) + values[n] >= 0)
            elif coeff_type == "float":
                cons.append(sum([float(lhs[i]) * values[i] for i in range(n)]) + float(cst) + values[n] >= 0)
            else: #sage rational type assumed
                cons.append(sum( [lhs[i] * values[i] for i in range(n)]) + float(cst) + values[n] >= 0)
        for linear_poly_expr in generate_symmetric_vertices_continuous_expr(h, bkpt[f_index]):
            lhs = [linear_poly_expr.coefficient(linear_poly_expr.parent().gens_dict()[name]) for name in coord_names]
            cst = linear_poly_expr.constant_coefficient()
            if coeff_type == "int":
                lcd = lcm(lcm([coeff.denominator() for coeff in lhs]), cst.denominator())
                cons.append(sum([int(lcd * lhs[i]) * values[i] for i in range(n)]) + int(lcd * cst) + values[n] == 1*lcd)
            elif coeff_type == "float":
                cons.append(sum([float(lhs[i]) * values[i] for i in range(n)]) + float(cst) + values[n] == 1)
            else: #sage rational type assumed
                cons.append(sum([lhs[i] * values[i] for i in range(n)]) + float(cst) + values[n] == 1)


    if not log_pw_functions:
        logging.getLogger("cutgeneratingfunctionology.igp.functions").setLevel(pw_logging_level)
    return cons

