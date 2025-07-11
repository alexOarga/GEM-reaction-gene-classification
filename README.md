# Linking Cellular Functions and Gene Essentiality through Growth Phenotype in Genome-Scale Models


This repository contains the implementation of the methods presented in our paper:

> ** Linking Cellular Functions and Gene Essentiality through Growth Phenotype in Genome-Scale Models**

The scripts reproduce the main experiments and results presented in the study, which introduces a novel method to classify metabolic reactions based on their functional role in the production of growth. These tools leverage Flux Balance Analysis (FBA) and Flux Variability Analysis (FVA) to:

- Identify essential and redundant reactions and dead reactions at optimal or suboptimal growth.
- Compute bottleneck reactions limiting optimal growth.
- Explain gene essentiality in terms of reaction redundancy.

## 📁 Repository Contents

- ```growth_limiting_reactions.py```
Coomputes the number of bottleneck reactions limiting optimal growth 
using Flux Variability Analysis (FVA) as described in the paper.

- ```generate_sets.py```
Processes a list of Genome-Scale Models (GEMs) and computes the functional classification of reactions for each model (essential, redundant, unused).

- ```generate_single.py```
Computes the functional classification of reactions for a single GEM
for different ranges of growth rates, identifying essential, redundant, and unused reactions.

## ⚙️ Requirements

This codebase is written in Python and uses the following main dependencies:

- [COBRApy](https://opencobra.github.io/cobrapy/)
- [CONTRABASS](https://github.com/openCONTRABASS/CONTRABASS)
- [Pyomo](https://www.pyomo.org/)
- [Gurobi Optimizer](https://www.gurobi.com/) (version ≥ 9.5.2)

The essential dependencies can be installed with:
```bash
pip install -r requirements.txt
```

## 🧪 Reproducing the Results

1. Place your genome-scale models (in SBML or JSON format) in a directory.

2. For reactions limiting growth,
change the model path in `growth_limiting_reactions.py` to point to your model file and run:
 run:
    ```bash
    python growth_limiting_reactions.py
    ```
3. For a batch analysis of models,
change the model paths in `generate_sets.py` to point to your models directory and run:
    ```bash
    python generate_sets.py
    ```

4. For a single model analysis,
change the model path in `generate_single.py` to point to your model file and run:
    ```bash
    python generate_single.py
    ```

## 📄 License

This project is licensed under the GPL-3.0 License. See the [LICENSE](LICENSE) file for details.

## 📫 Contact

For questions please contact [alexOarga](https://github.com/alexOarga)