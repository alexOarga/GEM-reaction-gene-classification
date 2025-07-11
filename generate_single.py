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
from generate_sets import *


CONST_EPSILON = 1e-08
PATH = "<YOUR-MODELS-PATH-HERE>"       
MODELS = [
    #"iML1515",
    "MODEL1507180007"
]

COMPUTE_MINIMAL = True
PRINT_MINIMAL_ERROR = True
OUTPUT = "output.csv"


print("ModelID\tER\tRR\tDR\tFR\tRG_FBA\tMRG")
for i in np.arange(0.0, 1.05, 0.05):
    ER, RR, DR, FR, RG_FBA, MRG = computation_minimum_reduced_sets(PATH, MODELS[0], growth_rate=i)
    print(f"{MODELS[0]}\t{len(ER)}\t{len(RR)}\t{len(DR)}\t{len(FR)}\t{len(RG_FBA)}\t{len(MRG)}")
