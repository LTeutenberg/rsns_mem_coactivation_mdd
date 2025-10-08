🧠 Project Description

This repository contains the full analysis pipeline accompanying the manuscript
“Synergistic Co-Activation Probabilities of Large-Scale Resting-State Networks in Major Depressive Disorder.”

The code implements a pairwise maximum entropy model (MaxEnt) to estimate co-activation probabilities across seven large-scale resting-state networks (RSNs) and examines their association with clinical features of major depressive disorder (MDD). Analyses include both univariate (GLM) and multivariate (CCA) approaches to relate network co-activation states to diagnostic and symptom measures.

The workflow is modular and organized into eight parts, ranging from probability estimation to predictive validation. Each part can be run independently or executed as a complete pipeline via main.m.

Key features

Maximum entropy modeling of RSN co-activation states (128 binary patterns).

General linear models (GLMs) testing group and remission effects (HC vs. MDD).

Canonical correlation analysis (CCA) linking brain-state probabilities to HAMD-17 symptom items.

Stability and generalizability analyses using repeated cross-validation and subgroup matching.

Predictive evaluation of canonical variates in held-out data.

Fully configurable through a single configuration file (config.m).

Reproducibility

All analyses are implemented in MATLAB and are deterministic given the provided random seed. Data-dependent components (paths, column indices) can be modified in config.m, ensuring transparency and adaptability to other datasets.
