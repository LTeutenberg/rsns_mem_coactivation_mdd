function C = config_example()
% Copy to src/config.m and adjust.

C.paths.rsn_data_mat   = fullfile('..','data','rsn_data.mat');
C.paths.behavioral_mat = fullfile('..','data','behavioral_data.mat');

C.paths.outdir   = fullfile('..','results');
C.paths.out_prob = fullfile(C.paths.outdir,'probabilities');
C.paths.out_glm  = fullfile(C.paths.outdir,'glm');
C.paths.out_corr = fullfile(C.paths.outdir,'correlations');
C.paths.out_cca  = fullfile(C.paths.outdir,'cca');

% --- MaxEnt / binarization ---
C.maxent.threshold_binarize = 0.1;
C.maxent.ncells             = 7;
C.maxent.train_threshold    = 1.3;
C.maxent.binsize            = 0.1;
C.maxent.depth              = 25;

% --- Column indices in behavioral_data ---
C.cols.site_vars    = 7:8;
C.cols.age          = 29;
C.cols.gender       = 30;
C.cols.edu_years    = 31;
C.cols.group_HC_MDD = 44;
C.cols.group_status = 228;
C.cols.severity     = 73;

% --- CCA CV settings ---
C.cca.kfolds_stability   = 5;
C.cca.repeats_stability  = 50;
C.cca.kfolds_predictive  = 30;
C.cca.repeats_predictive = 50;
C.cca.k_use              = 2;

C.seed = 42;
end
