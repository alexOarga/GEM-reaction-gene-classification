import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd

import shutil
import sys
import os.path

from pyomo.environ import *
from pyomo.gdp import *
import pyomo.environ as pyo
import pandas as pd


CONST_EPSILON = 1e-08


'''
    param: cobrapy model
    type: cobra.core.model
    return: dictionary of flux distribution
    type: dict[cobra.reaction.id] = flux(float)
'''
def pyomo_flux_balance_analysis(model, BIOMASS, verbose = True):
    
    m = ConcreteModel()
    S = {}

    REACTIONS = []
    for r in model.reactions:
        REACTIONS.append(r.id)

    METABOLITES = [metboj.id for metboj in model.metabolites]

    # Initiate stoichiometry matrix
    for met in model.metabolites:
        for r in model.reactions:
            if r in met.reactions:
                S[met.id, r.id] = r.get_coefficient(met)
            else:
                S[met.id, r.id] = float(0.0)

    def s_init(model, metabolite, reaction):
        return S[metabolite, reaction]    

    m.S = Param(METABOLITES, REACTIONS, initialize=s_init, mutable=True)

    def flux_bounds(m, reaction):
        for rr_m in model.reactions: # search reaction by id
            if rr_m.id == reaction: 
                rr = rr_m
                break

        return (rr.lower_bound, rr.upper_bound)

    # decision variables
    #    flux of each reaction
    m.v = Var(REACTIONS, domain = Reals, bounds=flux_bounds)

    # Structural constraints
    #     steady state constraint
    m.steady_state_constraint = ConstraintList()
    r = None
    for met in METABOLITES:
        m.steady_state_constraint.add(sum([ m.v[r] * m.S[met, r] for r in REACTIONS]) == 0.0)

    # objective function
    #    maximize biomass flux
    m.objective = Objective(expr = m.v[BIOMASS],
                            sense = maximize)

    print("starting ...")
    results = SolverFactory('cplex').solve(m)
    if verbose:
        results.write()

    # generate flux result list [(cobra.reaction, flux(float))]
    solution = {}
    for r in model.reactions:
        solution[r.id] = m.v[r.id]()
            
    return solution


# This is a function to compare cobrapy and pyomo flux balance analysiss
def test_pyomo_flux_balance_analysis(cobra_model, objective_id, tolerance = 1e-6):
    pyomo_fba = pyomo_flux_balance_analysis(cobra_model, BIOMASS, verbose = False)
    cobra_fba = cobra_model.optimize()
    
    if abs(pyomo_fba[objective_id] - cobra_fba[objective_id]) > tolerance:
        raise AssertionError("Objective flux does not match: pyomo: {}, cobra: {}".format(pyomo_fba[objective_id], cobra_fba[objective_id] ))
        
    count = 0
    for r in model.reactions:
        if abs(pyomo_fba[r.id] - cobra_fba[r.id]) > tolerance:
            warn("{}.Reaction {} flux does not match: pyomo: {}, cobra: {}".format(count, r.id, pyomo_fba[r.id], cobra_fba[r.id]))
            count = count + 1
            
            

def pyomo_min(model, BIOMASS, max_growth, flex, verbose = True):
    """
    Solves a linear optimization problem using Pyomo for minimizing a set of flexible reactions 
    subject to constraints on growth and fluxes.

    Parameters:
    - model: A metabolic model object.
    - BIOMASS: Identifier for the biomass reaction.
    - max_growth: Maximum growth rate constraint.
    - flex: List of flexible reactions.
    - verbose: (optional) If True, prints solver results (default is True).

    Returns:
    - solution: Dictionary containing flux values for non-zero flux reactions.
    - solution_g: Dictionary containing binary values for flexible reactions.
    - solution_w: Dictionary containing flux bounds for all reactions.
    """
    m = ConcreteModel()
    S = {}

    REACTIONS = []
    for r in model.reactions:
        REACTIONS.append(r.id)

    METABOLITES = [metboj.id for metboj in model.metabolites]
    FLEXIBLE = [r.id for r in flex]
    
    # Initiate stoichiometry matrix
    for met in model.metabolites:
        for r in model.reactions:
            if r in met.reactions:
                S[met.id, r.id] = r.get_coefficient(met)
            else:
                S[met.id, r.id] = float(0.0)

    def s_init(model, metabolite, reaction):
        return S[metabolite, reaction]    

    m.S = Param(METABOLITES, REACTIONS, initialize=s_init, mutable=True)

    def flux_bounds(m, reaction):
        for rr_m in model.reactions: # search reaction by id
            if rr_m.id == reaction: 
                rr = rr_m
                break

        if rr.id == BIOMASS:
            #print(f"setting bounds for growth: ({rr.lower_bound}, {rr.upper_bound})")
            return (rr.lower_bound, rr.upper_bound)
                
        if rr.lower_bound > rr.upper_bound:
            raise RuntimeError(f"Lower flux greater than upper in reaction: {rr.id}")
                
        return (rr.lower_bound, rr.upper_bound)

    # decision variables
    #    flux of each reaction
    m.w = Var(REACTIONS, domain = Reals, bounds=flux_bounds)
    
    # decision variables
    #    final flux
    m.v = Var(REACTIONS, domain = Reals, bounds=flux_bounds)

    # boolean variables to minimise
    #m.g = Var(FLEXIBLE, domain = Boolean)
    m.g = pyo.Var(FLEXIBLE, within=pyo.Binary)
    
    # Structural constraints
    #     steady state constraint
    m.steady_state_constraint = ConstraintList()
    r = None
    for met in METABOLITES:
        m.steady_state_constraint.add(sum([ m.v[r] * m.S[met, r] for r in REACTIONS]) == 0.0)

    # constrain
    #  keep maximum flux
    #print(f"set constraint: v[{BIOMASS}]={max_growth}")
    m.max_growth = Constraint(expr = m.v[BIOMASS] >= max_growth)
        
    # constraints product v = w*g
    #     flux_boolean_constraint
    m.flux_boolean_constraint = ConstraintList()
    for r in FLEXIBLE:
        lower, upper = flux_bounds(m, r)
        #m.flux_boolean_constraint.add( m.v[r] = m.w[r] * m.g[r] )
        m.flux_boolean_constraint.add( m.v[r] <= upper * m.g[r] )
        m.flux_boolean_constraint.add( m.v[r] >= lower * m.g[r] ) 
        m.flux_boolean_constraint.add( m.v[r] <= m.w[r] - (lower * (1.0 - m.g[r])) )
        m.flux_boolean_constraint.add( m.v[r] >= m.w[r] - (upper * (1.0 - m.g[r])) )
        
    # objective function
    m.objective = Objective(expr = sum([m.g[r] for r in FLEXIBLE]),
                            sense = minimize)

    #print("starting ...")
    solver = SolverFactory('gurobi')
    #solver.options['MIPFocus'] = 3
    solver.options['Heuristics'] = 0.1
    solver.options['TimeLimit'] = 600 # 10 minutes
    #solver.options['MIPGap'] = 5e-3
    #solver.options['IntFeasTol'] = 1e-7
    #solver.options['mip limits solutions'] = 1
    #results = solver.solve(m, tee=True)
    results = solver.solve(m)
    if verbose:
        results.write()

    '''
    if (results.solver.status != SolverStatus.ok) or (results.solver.termination_condition != TerminationCondition.optimal):
        return None
    '''
        
    # generate flux result list [(cobra.reaction, flux(float))]
    solution = {}
    solution_g = {}
    solution_w = {}
    for r in model.reactions:
        if m.v[r.id]() > CONST_EPSILON or m.v[r.id]() < -CONST_EPSILON:
            solution[r.id] = m.v[r.id]()
        else:
            solution[r.id] = 0.0
            
        if r.id in FLEXIBLE:
            solution_g[r.id] = m.g[r.id]()
            solution_w[r.id] = m.w[r.id]()
            
            
    return solution, solution_g, solution_w
