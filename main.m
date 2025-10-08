function main()
% Main wrapper script for all analysis parts (1–8)
% Author: [Your Name]
% Project: BP:CNNI manuscript
% Paper: "Synergistic Co-Activation Probabilities of Large-Scale Resting State Networks in Major Depressive Disorder"

addpath(genpath(fileparts(mfilename('fullpath'))));

% Load configuration
if exist('config.m','file')
    C = config();
else
    C = config_example();
end

% Create output directory if not present
if ~exist(C.paths.outdir,'dir')
    mkdir(C.paths.outdir);
end

% === Load data ===
S1 = load(C.paths.rsn_data_mat, 'RSN_data');             
RSN_data = S1.RSN_data;

S2 = load(C.paths.behavioral_mat, 'behavioral_data');    
behavioral_data = S2.behavioral_data;

% === Set random seed for reproducibility ===
seed_rng(C.seed);

% === Part 1 ===
fprintf('\n[1/8] Maximum Entropy Model...\n');
[Energy_states, Prob_states] = part1_maxent_probabilities(RSN_data, C);

% === Part 2 ===
fprintf('\n[2/8] GLM: HC vs MDD...\n');
T_group = part2_glm_group(Prob_states, behavioral_data, C);

% === Part 3 ===
fprintf('\n[3/8] GLM: HC vs MDD (acute/remitted)...\n');
T_status = part3_glm_status(Prob_states, behavioral_data, C);

% === Part 4 ===
fprintf('\n[4/8] Partial correlations with symptom severity...\n');
T_partial = part4_partialcorr(Prob_states, behavioral_data, C);

% === Part 5 ===
fprintf('\n[5/8] Canonical Correlation Analysis (full sample)...\n');
[A_full, B_full, r_full, stats_full, X, Y, validIdx] = part5_cca_full(Prob_states, behavioral_data, C);

% === Part 6 ===
fprintf('\n[6/8] CCA Stability Analysis...\n');
part6_cca_stability(X, Y, A_full, C);

% === Part 7 ===
fprintf('\n[7/8] CCA Generalizability (subgroups)...\n');
results_subgroups = part7_cca_generalizability(Prob_states, behavioral_data, C);

% === Part 8 ===
fprintf('\n[8/8] Predictive Utility of CCA...\n');
part8_cca_predictive(X, Y, A_full, C);

fprintf('\nAll analyses completed successfully.\nResults saved to: %s\n', C.paths.outdir);
end
