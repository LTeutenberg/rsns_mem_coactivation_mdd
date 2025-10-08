🧠 Project Description

This repository contains the full analysis pipeline accompanying the manuscript
“Synergistic Co-Activation Probabilities of Large-Scale Resting-State Networks in Major Depressive Disorder.”

The code implements a pairwise maximum entropy model (MaxEnt) to estimate co-activation probabilities across seven large-scale resting-state networks (RSNs) and examines their association with clinical features of major depressive disorder (MDD). Analyses include both univariate (GLM) and multivariate (CCA) approaches to relate network co-activation states to diagnostic and symptom measures.

The workflow is modular and organized into eight parts, ranging from probability estimation to predictive validation. Each part can be run independently or executed as a complete pipeline via main.m.

🧩Key features

Maximum entropy modeling of RSN co-activation states (128 binary patterns).

General linear models (GLMs) testing group and remission effects (HC vs. MDD).

Canonical correlation analysis (CCA) linking brain-state probabilities to depression severity.

Stability and generalizability analyses for CCA using repeated cross-validation and subgroup matching.

Predictive evaluation of canonical variates in held-out data.

Fully configurable through a single configuration file (config.m).

Synthetic validation of the maximum entropy model and canonical correlation analysis.

Reproducibility

All analyses are implemented in MATLAB. Data-dependent components (paths, column indices) can be modified in config.m, ensuring transparency and adaptability to other datasets. Data used in the original paper cannot be shared publicly.

🧩 Requirements

To run the analyses, the following software and dependencies are required:

MATLAB (R2021b or later)
Required for all analyses and scripts. The code has been tested with MATLAB versions R2021b–R2024a.

Statistics and Machine Learning Toolbox
Used for GLM fitting (fitglm), canonical correlation analysis (canoncorr), and cross-validation (cvpartition).

MaxEnt Toolbox
Required for estimating pairwise maximum entropy (Ising) models of network co-activation states.
→ Available at: https://orimaoz.github.io/maxent_toolbox/
