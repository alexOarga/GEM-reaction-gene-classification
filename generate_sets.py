#%matplotlib inline
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd

import shutil
import sys
import os, os.path
            
assert(shutil.which("cbc") or os.path.isfile("cbc"))
from pyomo.environ import *
from pyomo.gdp import *
import pandas as pd

import numpy as np
import sys
import xlwt
import warnings
import time
import math
import networkx as nx
import matplotlib.pyplot as plt
import cobra
from cobra.flux_analysis import *
from cobra.flux_analysis.deletion import single_reaction_deletion
from enum import Enum
from math import isnan

from contrabass import CobraMetabolicModel
from pyomo_utils import pyomo_flux_balance_analysis, pyomo_min



CONST_EPSILON = 1e-08

PATH = "<YOUR-MODELS-PATH HERE>"

COMPUTE_MINIMAL = True
PRINT_MINIMAL_ERROR = True
OUTPUT = "output.csv"

growth_rate = 1.0 # 1.0 for optimal growth, <1.0 for suboptimal

MODELS = [
    "MODEL1507180052",
    "MODEL1507180006",
    "MODEL1106080000",
    "MODEL1507180007",
    "MODEL1507180030",
    "MODEL1507180048",
    "MODEL1507180070",
    "MODEL1710040000",
    "MODEL1507180024",
    "MODEL1507180036",
    "MODEL1507180021",
    "MODEL1507180044",
    "MODEL1507180049",
    "MODEL1507180068",
    "MODEL1507180060",
    "MODEL1507180059",
    "MODEL1507180013",
    #"MODEL1507180058",
    "MODEL1507180022",
    "MODEL1507180012",
    "MODEL1507180027",
    "MODEL1507180011",
    "MODEL1507180033",
    "MODEL1507180015",
    "MODEL1507180064",
    "MODEL1212060001",
    "MODEL1507180054",
    "MODEL1105030000",
    "MODEL1507180010",
    "MODEL1507180017"
]



def current_milli_time():
    return round(time.time() * 1000)

class Direction(Enum):
    FORWARD = 0
    BACKWARD = 1
    REVERSIBLE = 2
    
def reaction_direction(reaction):
    if - CONST_EPSILON < reaction.upper_bound < CONST_EPSILON:
        upper = 0
    else:
        upper = reaction.upper_bound
    if - CONST_EPSILON < reaction.lower_bound < CONST_EPSILON:
        lower = 0
    else:
        lower = reaction.lower_bound

    if lower == 0 and upper == 0:
        # Note that reactions with upper and lower bound equal to 0 are considered FORWARD by default
        return Direction.FORWARD
    elif lower >= 0:
        return Direction.FORWARD
    elif upper > 0:
        return Direction.REVERSIBLE
    else:
        return Direction.BACKWARD

def is_dead_reaction(reaction):
    return abs(reaction.upper_bound) < CONST_EPSILON and abs(reaction.lower_bound) < CONST_EPSILON
    
def metabolite_unique_producer_unique_consumer(findCPcore_model, metabolite):
    
    model = findCPcore_model
    
    unique_producer = None
    unique_consumer = None
    
    for r_cp, m_cp in model.chokepoints():
        if m_cp.id == metabolite.id:
            if metabolite.id in [aux.id for aux in reaction_reactants(r_cp)]:
                unique_producer = copy(r_cp.id)
            elif metabolite.id in [aux.id for aux in reaction_products(r_cp)]:
                unique_consumer = copy(r_cp.id)
            else:
                raise RuntimeError("Error not reactions...")
                
    if unique_producer is not None and unique_consumer is not None:
        unique_producer = model.model().reactions.get_by_id(unique_producer)
        unique_consumer = model.model().reactions.get_by_id(unique_consumer)
        return unique_producer, unique_consumer
    else:
        return None, None
    
    
def flux_dependent_reactants_products(reaction):
    if reaction_direction(reaction) == Direction.FORWARD:
        reactants1 = reaction.reactants
        products1  = reaction.products
    elif reaction_direction(reaction) == Direction.REVERSIBLE:
        reactants1 = reaction.reactants
        products1  = reaction.products
    elif reaction_direction(reaction) == Direction.BACKWARD:
        reactants1 = reaction.products
        products1  = reaction.reactants
    else:
        raise RuntimeError("Error flux_dependent_reactants_products()")
        
    return reactants1, products1


def flux_dependent_reactants_products_doubling_reversibles(reaction):   
    if reaction_direction(reaction) != Direction.REVERSIBLE:
        return flux_dependent_reactants_products(reaction)
    else:
        reactants = reaction.reactants
        products  = reaction.products
        combined = reactants + products
        return combined, combined

    
def flux_dependent_producer_consumers(m):
    producers = set()
    consumers = set()
    for r in m.reactions:
        if not is_dead_reaction(r):
            if reaction_direction(r) == Direction.REVERSIBLE:
                producers.add(r)
                consumers.add(r)
            else:
                reactants, products = flux_dependent_reactants_products(r)

                if m in reactants:
                    consumers.add(r)
                elif m in products:
                    producers.add(r)
                else:
                    raise RuntimeError("code error!")
    return producers, consumers

def reactions_optimal_growth(solution, model):
    ROG_id = set()
    for r in model.reactions:
        if solution[r.id] > CONST_EPSILON or solution[r.id] < -CONST_EPSILON:
            ROG_id.add(r.id)
    return ROG_id

def _get_biomass_reaction(cobra_model):
    # NOTE: We dont expect to have more than one objective reaction
    return list(cobra.util.solver.linear_reaction_coefficients(cobra_model).keys())[0]

def is_open_flux_reaction(r):
    return r.upper_bound > -CONST_EPSILON and r.lower_bound < CONST_EPSILON

def is_forced_flux_reaction(r):
    return r.lower_bound > CONST_EPSILON or r.upper_bound < -CONST_EPSILON

def is_loose_flux_reaction(r):
    return abs(r.upper_bound - r.lower_bound) > CONST_EPSILON

def is_tight_flux_reaction(r):
    return abs(r.upper_bound - r.lower_bound) < CONST_EPSILON


def growth_dependent_essential_reactions(cobra_model, growth_rate=1.0):
    biomass = _get_biomass_reaction(cobra_model)
    max_growth = cobra_model.optimize()[biomass.id]
    max_growth_rate = max_growth * growth_rate
    essential_reactions = set()
    deletions = single_reaction_deletion(cobra_model, method='fba')
    reactions_knockout_growth = deletions.loc[:, :]['growth']
    for r, g in reactions_knockout_growth.items():
        reaction = cobra_model.reactions.get_by_id(list(r)[0])
        if isnan(g) or g < CONST_EPSILON or (g + CONST_EPSILON < max_growth_rate):
            essential_reactions.add(reaction)
    return essential_reactions



def count_reactions(model):
    DR_1 = [r for r in model.reactions() if is_dead_reaction(r)]
    FR_1 = [r for r in model.reactions() if is_forced_flux_reaction(r)]
    RR_1 = [r for r in model.reactions() if not is_dead_reaction(r) and reaction_direction(r) == Direction.REVERSIBLE]
    NR_1 = [r for r in model.reactions() if not is_dead_reaction(r) 
            and (reaction_direction(r) == Direction.FORWARD 
            or reaction_direction(r) == Direction.BACKWARD)]
    
    OR_1 = [r for r in model.reactions() if not is_dead_reaction(r) and is_open_flux_reaction(r)]
    TR_1 = [r for r in model.reactions() if not is_dead_reaction(r) and is_tight_flux_reaction(r)]
    LR_1 = [r for r in model.reactions() if not is_dead_reaction(r) and is_loose_flux_reaction(r)]
    assert len(OR_1) < len(LR_1)
    for r in RR_1: # all RR in OR and LR
        assert r in OR_1
        assert r in LR_1
    for r in OR_1: # all OR in LR and all OR not in TR
        assert r in LR_1
        assert r not in TR_1
    for r in TR_1: # TR_1 and LR_1 are disjoint
        assert r not in LR_1
    for r in LR_1:
        assert r not in TR_1
    # disjoint sets
    assert len(DR_1) + len(RR_1) + len(NR_1) == len(model.reactions())
    for r in DR_1:
        assert  r not in RR_1 and r not in NR_1 \
            and r not in OR_1 and r not in LR_1 \
            and r not in TR_1
    print("NR_1 - OR_1: ", len(set(NR_1).difference(set(OR_1))))
    return DR_1, RR_1, NR_1, OR_1, TR_1, LR_1, FR_1
    

def build_EROG_model(EROG, biomass):
    # generate new model for bottom-up construction
    from cobra import Model
    step = Model('base-model')
    step.add_reactions(list(EROG))
    step.objective = biomass.id
    return step
        
        
def compute_reactions_sets(MODEL, growth_rate=1.0):
    cobra.core.model.configuration.solver = 'glpk'
    cbmodel_fva = CobraMetabolicModel(MODEL)
    cbmodel_fva.set_epsilon(CONST_EPSILON)
    #cbmodel_fva.set_objective(RG_ID)
    cbmodel_fva.model().optimize()
    
    forced_ids = [r.id for r in cbmodel_fva.reactions() if is_forced_flux_reaction(r)]
    
    cbmodel_fva.fva(update_flux = True, threshold=growth_rate)
    model = cbmodel_fva.model()

    #cbmodel_fva.print_model_info()

    dead_reactions_ids = set()
    for r in model.reactions:
        if is_dead_reaction(r):
            dead_reactions_ids.add(r.id)

    cobra.core.model.configuration.solver = 'glpk'
    cbmodel = CobraMetabolicModel(MODEL)
    cbmodel.set_epsilon(CONST_EPSILON)
    #cbmodel.set_objective(RG_ID)
    #cbmodel.fva(update_flux = True, threshold=0.0)
    model = cbmodel.model()
    biomass = _get_biomass_reaction(cbmodel.model())
    max_growth = cbmodel.model().optimize()[biomass.id]

    #cbmodel.knockout_reactions_growth()
    #cbmodel.find_optimal_growth_essential_reactions()
    #essential_reactions = cbmodel.optimal_growth_essential_reactions().union(set([biomass]))
    essential_reactions = growth_dependent_essential_reactions(model, growth_rate=growth_rate)
    essential_reactions = essential_reactions.union(set([biomass]))
    essential_reactions_ids = [r.id for r in essential_reactions]

    # removing blocked reactions from the model
    cbmodel.model().remove_reactions(dead_reactions_ids)
    # removing metabolites with no connection
    metabolites_to_delete = []
    for m in cbmodel.model().metabolites:
        if len(m.reactions) == 0:
            metabolites_to_delete.append(m)
    cbmodel.model().remove_metabolites(metabolites_to_delete)

    # add redundant reactions
    redundant_reactions_ids = set([r.id for r in cbmodel.model().reactions]) - set(essential_reactions_ids)

    if growth_rate < 1.0:
        # we set a limit to the biomass upper bound to simulate a suboptimal growth
        biomass.upper_bound = max_growth * growth_rate 
        #print("New biomass upper bound: ", biomass.upper_bound)
    solution = cbmodel.model().optimize()
    RG_FBA_id = reactions_optimal_growth(solution, cbmodel.model())
    #pfba_sol = cobra.flux_analysis.pfba(cbmodel.model())
    #ROG_pFBA_id = reactions_optimal_growth(pfba_sol, cbmodel.model())
    #print("ROG (pFBA):", len(ROG_pFBA_id))
    
    return cbmodel, cbmodel_fva, biomass, \
        essential_reactions_ids, redundant_reactions_ids, dead_reactions_ids, forced_ids, \
        RG_FBA_id
        

def computation_minimum_reduced_sets(PATH, M, growth_rate=1.0):
    start = current_milli_time()
    MODEL = PATH + M + ".xml" 
    print("    Computing set of reactions for growth in model: ", M)
    cbmodel, cbmodel_fva, biomass, \
        essential_reactions_ids, redundant_reactions_ids, dead_reactions_ids, forced_ids, \
        RG_FBA_id = compute_reactions_sets(MODEL, growth_rate=growth_rate)
    
    # Compute max growth
    max_growth = cbmodel.model().optimize()[biomass.id]
    #pyomo_fba = pyomo_flux_balance_analysis(cbmodel_fva.model(), BIOMASS, verbose = False)
    #max_growth = pyomo_fba[BIOMASS]
    assert max_growth >= 0, "growth is : " + str(max_growth)
    #print("    Reactions: ", len(cbmodel.model().reactions))
    #print("    Metabolites: " ,len(cbmodel.model().metabolites))
    #print("    Growth: ", max_growth)

    print("    Computing minimum set of reactions for model: ", M)
    redundant_reactions = [cbmodel.model().reactions.get_by_id(rid) for rid in redundant_reactions_ids]
    res_min = pyomo_min(cbmodel.model(), biomass.id, max_growth, redundant_reactions, verbose = False)
    if res_min is None:
        print("    #### INFEASIBLE ####    ")
        return
    
    sol_v, sol_g, sol_w = res_min    
    min_rg_ids = [k for k, v in sol_g.items() if v > 0.5]
    #print("    Length min variables: ", len(min_rg_ids))
    #print("    Length ER: ", len(essential_reactions_ids))
    #print("    Length MRG: ", len(min_rg_ids) + len(essential_reactions_ids))

    minimum_set_growth_ids = essential_reactions_ids.copy() + min_rg_ids.copy()

    return essential_reactions_ids, redundant_reactions_ids, dead_reactions_ids, forced_ids, \
        RG_FBA_id, minimum_set_growth_ids
        
        

if __name__ == "__main__":

    print("ModelID\tER\tRR\tDR\tFR\tRG_FBA\tMRG")

    for i in range(0, len(MODELS)):
        last_time = time.time()
        ER, RR, DR, FR, RG_FBA, MRG = computation_minimum_reduced_sets(PATH, MODELS[i], growth_rate=growth_rate)
        now = time.time() - last_time
        print(f"{MODELS[i]}\t{len(ER)}\t{len(RR)}\t{len(DR)}\t{len(FR)}\t{len(RG_FBA)}\t{len(MRG)}")
        print(now % 60, " seconds ")
