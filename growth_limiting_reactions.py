import csv
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

PATH = "<YOUR-MODELS-PATH-HERE>"
TOLERANCE=1E-9

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
    "MODEL1507180058",
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


def compare_bounds(reactions_list, L, U):
  
  reactions_equal_U = []  # reactions whose flux bounds are lb_1 = ub_1 = U
  reactions_equal_L = []  # reactions whose flux bounds are lb_1 = ub_1 = L

  for r in reactions_list:

    # if flux bounds are not zero 
    if not math.isclose(r.upper_bound, 0.0, abs_tol=TOLERANCE): 
    
      # if (lb_1 = ub_1 = L) OR (lb_1 = ub_1 = U) 
      if math.isclose(r.upper_bound, r.lower_bound, abs_tol=TOLERANCE) \
      and math.isclose(r.upper_bound, U[r], abs_tol=TOLERANCE):

        #print_reaction_bounds(r, "with lb_1=ub_1=U")
        reactions_equal_U.append(r)

      elif math.isclose(r.upper_bound, r.lower_bound, rel_tol=TOLERANCE) \
      and math.isclose(r.lower_bound, L[r], abs_tol=TOLERANCE):

        #print_reaction_bounds(r, "with lb_1=ub_1=L")
        reactions_equal_L.append(r)
        
  return reactions_equal_U, reactions_equal_L


with open('fva.csv', mode='w') as csv_file:
    fieldnames = ['Biomass model ID', 'upper', 'lower']
    writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
    writer.writeheader()

    for i, MODEL in enumerate(MODELS):
        print(i, MODEL)
        MODEL_FILENAME = PATH + MODEL + ".xml"
        model = CobraMetabolicModel(MODEL_FILENAME)

        U = dict([(r, r.upper_bound) for r in model.reactions()])
        L = dict([(r, r.lower_bound) for r in model.reactions()])

        # run FVA and constrain model flux
        model.fva(update_flux=True, threshold=1.0)

        reactions_equal_U, reactions_equal_L = compare_bounds(model.reactions(), L, U)

        writer.writerow({'Biomass model ID': MODEL, 'upper': str(len(reactions_equal_U)), 'lower': str(len(reactions_equal_L))})


