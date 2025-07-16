import argparse
import csv
import numpy as np
import sys, os
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

def compare_bounds(reactions_list, L, U, tolerance):
    reactions_equal_U = []
    reactions_equal_L = []

    for r in reactions_list:
        if not math.isclose(r.upper_bound, 0.0, abs_tol=tolerance):
            if math.isclose(r.upper_bound, r.lower_bound, abs_tol=tolerance) \
            and math.isclose(r.upper_bound, U[r], abs_tol=tolerance):
                reactions_equal_U.append(r)
            elif math.isclose(r.upper_bound, r.lower_bound, rel_tol=tolerance) \
            and math.isclose(r.lower_bound, L[r], abs_tol=tolerance):
                reactions_equal_L.append(r)
    
    return reactions_equal_U, reactions_equal_L


def main():
    parser = argparse.ArgumentParser(description="Run FVA to compute the number of growth limitting reactions (on upper or lower bound).")
    parser.add_argument("--models", nargs='*', help="List of model path to analyze")
    parser.add_argument("--tolerance", type=float, default=1e-9, help="Tolerance for bound comparison")
    parser.add_argument("--output", default="fva.csv", help="CSV file to write the results to")

    args = parser.parse_args()

    default_models = [
        "MODEL1507180052.xml", "MODEL1507180006.xml", "MODEL1106080000.xml", "MODEL1507180007.xml",
        "MODEL1507180030.xml", "MODEL1507180048.xml", "MODEL1507180070.xml", "MODEL1710040000.xml",
        "MODEL1507180024.xml", "MODEL1507180036.xml", "MODEL1507180021.xml", "MODEL1507180044.xml",
        "MODEL1507180049.xml", "MODEL1507180068.xml", "MODEL1507180060.xml", "MODEL1507180059.xml",
        "MODEL1507180013.xml", "MODEL1507180058.xml", "MODEL1507180022.xml", "MODEL1507180012.xml",
        "MODEL1507180027.xml", "MODEL1507180011.xml", "MODEL1507180033.xml", "MODEL1507180015.xml",
        "MODEL1507180064.xml", "MODEL1212060001.xml", "MODEL1507180054.xml", "MODEL1105030000.xml",
        "MODEL1507180010.xml", "MODEL1507180017.xml"
    ]

    models_to_run = args.models if args.models else default_models
    for path in models_to_run:
        if not os.path.exists(path):
            raise FileNotFoundError(f"Path {path} does not exist.")
        if not path.endswith(".xml"):
            raise ValueError("Only XML model files are supported.")
        # separate directory and file name

    with open(args.output, mode='w', newline='') as csv_file:
        fieldnames = ['Biomass model ID', 'upper', 'lower']
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()

        for i, model_path in enumerate(models_to_run):
            print(i, model_path)
            model = CobraMetabolicModel(model_path)

            U = {r: r.upper_bound for r in model.reactions()}
            L = {r: r.lower_bound for r in model.reactions()}

            model.fva(update_flux=True, threshold=1.0)

            reactions_equal_U, reactions_equal_L = compare_bounds(model.reactions(), L, U, args.tolerance)

            writer.writerow({
                'Biomass model ID': model_path,
                'upper': str(len(reactions_equal_U)),
                'lower': str(len(reactions_equal_L))
            })


if __name__ == "__main__":
    main()
