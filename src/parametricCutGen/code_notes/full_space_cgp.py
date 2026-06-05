def lipschitz_constant_linear_constraints(bkpt, f_index, val=None, M=1e8, *,coeff_type='int', backend=None):
    n = len(bkpt)
    if backend == 'pplite':
        Var = pplite_Variable
    else:
        Var = Variable
    if val is None:
        val = [Var(i) for i in range(n)]
    val_ext = val+[0]
    cons = []
    ext_bkpt = bkpt+[1]
    for i in range(n):
        if coeff_type == "int":
            lcd = lcm(ext_bkpt[i].denominator(), ext_bkpt[i+1].denominator()) 
            cons.append(lcd*(val_ext[i+1]-val_ext[i]) <= M*(lcd*(ext_bkpt[i+1]-ext_bkpt[i]))  )
        elif coeff_type == "float":
            cons.append(lcd*(val_ext[i+1]-val_ext[i]) <= M*(ext_bkpt[i+1]-ext_bkpt[i])  )
        else: #sage rational type assumed
            cons.append(val_ext[i+1]-val_ext[i] <= M*(ext_bkpt[i+1]-ext_bkpt[i])  )
    return cons


# Outline for full space algorithm.
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
    if self._objective_sense == "maximize": # cgp written as max cut_score(cut) st. cut \in cutSpace
        self._cut_score.set_objective_sense("minimize") # solvers are by default in the of min f(x) s.t. x\in S =  max -f(x) s.t. x\in S 
        best_result = -1*np.inf                         # f(x) = -1*cut_score(point) with objective sense min. 
    elif self._objective_sense == "minimize":           # f(x) = cut_score(point) with objective sense max.
        self._cut_score.set_objective_sense("maximize")
        best_result = np.inf
    best_point = None
    if self._cut_space is None: # load the semi algebraic cell descriptions of pi min.
         self._cut_space = PiMinContContainer(self._max_num_of_bkpts, backend=self._backend)
    # start the clock when the actual portion of the solving processes starts.
    problem_timer = cgpTimer(self._max_cgp_solver_time)
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
                cut_generation_problem_logger.debug(f"Non empty bsa found.")
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
                    continue_solving = True
                    self._solver.nonlinear_solve(cut_score, point, subdomain_solver_constraints)
                except ModelViolation: # _cut_score.get_feasible_point() for last feasible point.
                    cut_generation_problem_logger.debug(f"Cell {subdomain_with_f_constraint} solved with exit {ModelViolation}")
                except SolverTimeOut:
                    continue_solving = False
                    cut_generation_problem_logger.error(f"Cell {subdomain_with_f_constraint} not solved fully with exit {SolverTimeOut}")
                except SolverRelTolReached:
                    cut_generation_problem_logger.debug(f"Cell {subdomain_with_f_constraint} solved with exit {SolverRelTolReached}")
                # Get the most recent feasible point and value according to cutScore
                point = self._cut_score.get_feasible_point()  # this should always be defined.
                result = self._cut_score.get_prev_result() # this should always be defined.
                self._cut_score.set_prev_result(None)
                self._cut_score.set_feasible_point(None)
                cut_generation_problem_logger.debug(f"Cell {subdomain_with_f_constraint} result:{result}, best_result:{best_result}")
                if self._objective_sense == "maximize": # multiply resutl by -1 to get objective value in original problem phrasing.
                    if best_result < -1*result:
                        best_result = -1*result
                        best_point = point
                elif self._objective_sense == "minimize":  # f(x) = cut_score(point) with objective sense max.
                    if best_result > result:
                        best_result = result
                        best_point = point
                if not continue_solving:
                    break
        except EmptyBSA:
            cut_generation_problem_logger.debug(f"BSA EMPTY.")
            continue
    # If best_result is +/- np.inf, the solver has failed to find any meaningful result.
    # The statement below should always be false.
    if best_result == np.inf or best_result == -1*np.inf:
        cut_generation_problem_logger.error(f"The solver has failed.  best_result is unbouded {best_result}")
        raise SolverError("The solver has failed.  best_result is unbouded {best_result}")
    sol_vals = [QQ(gamma_i) for gamma_i in best_point[self._max_num_of_bkpts:]]
    sol_bkpts = [QQ(lambda_i) for lambda_i in best_point[:self._max_num_of_bkpts]]
    pi_p = piecewise_function_from_breakpoints_and_values(sol_bkpts+[1],sol_vals+[0])
    log_problem_result(sol_bkpts, sol_vals, binvarow, binvc, f)
    if self._prove_seperator:
        res = minimality_test(pi_p) # add someway to log certificates.
        if not res:
            cut_generation_problem_logger.error(f"minimality of  {pi_p}: {res}")
    return pi_p, best_result
