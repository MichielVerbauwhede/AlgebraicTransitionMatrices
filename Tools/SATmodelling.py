from pysat.solvers import Glucose4 as Solver

def enum_projected_models(clauses, projected_vars, assumptions=[]):
    res = []
    with Solver(bootstrap_with=clauses) as solver:
        while solver.solve(assumptions=assumptions):
            model = solver.get_model()
            solver.add_clause(tuple(-model[v-1] for v in projected_vars))
            yield tuple(model[v-1] for v in projected_vars)