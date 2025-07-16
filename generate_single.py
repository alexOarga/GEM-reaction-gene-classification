import argparse
import numpy as np
import shutil
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import warnings
import math
import networkx as nx
import cobra
from cobra.flux_analysis import *
from cobra.flux_analysis.deletion import single_reaction_deletion
from enum import Enum
from math import isnan

from contrabass import CobraMetabolicModel
from pyomo_utils import pyomo_flux_balance_analysis, pyomo_min
from generate_sets import *
from pyomo.environ import *
from pyomo.gdp import *


def main():
    CONST_EPSILON = 1e-08 # Default epsilon value

    parser = argparse.ArgumentParser(description="Run functional classification of reactions for a single GEM")
    parser.add_argument("path", help="Path to the model files")
    #parser.add_argument("--compute-minimal", action="store_true", help="Whether to compute minimal sets")
    #parser.add_argument("--print-minimal-error", action="store_true", help="Whether to print minimal error")
    #parser.add_argument("--output", default="output.csv", help="Output CSV filename")
    parser.add_argument("--epsilon", type=float, default=CONST_EPSILON, help="Epsilon constant used for numerical tolerance")
    
    args = parser.parse_args()
    CONST_EPSILON = args.epsilon  # Set it based on CLI input

    # Check for CBC solver
    assert(shutil.which("cbc") or os.path.isfile("cbc")), "CBC solver not found"

    path = args.path
    if not os.path.exists(path):
        raise FileNotFoundError(f"Path {path} does not exist.")
    if not path.endswith(".xml"):
        raise ValueError("Only XML model files are supported.")
    # separate directory and file name
    model_dir = os.path.dirname(path)
    model_name = os.path.basename(path)
    # get model_name without extension
    model_name = os.path.splitext(model_name)[0]

    print("ModelID\tER\tRR\tDR\tFR\tRG_FBA\tMRG")
    for i in np.arange(0.0, 1.05, 0.05):
        ER, RR, DR, FR, RG_FBA, MRG = computation_minimum_reduced_sets(
            model_dir, model_name, growth_rate=i)
        print(f"{model_name}\t{len(ER)}\t{len(RR)}\t{len(DR)}\t{len(FR)}\t{len(RG_FBA)}\t{len(MRG)}")

if __name__ == "__main__":
    main()