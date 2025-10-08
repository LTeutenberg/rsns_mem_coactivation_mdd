%%%%%%   SYNTHETIC MAXENT + DEMOGRAPHICS + CCA PIPELINE
% Author: Lea Teutenberg, Hamidreza Jamalabadi - University of Marburg - October 2025
% Purpose:
% - Generate synthetic binary symptom datasets across multiple subjects under a ground-truth pairwise MaxEnt (Ising) model.
% - Create demographic covariates (age, gender, etc.).
% - Fit a MaxEnt model (k-pairwise/Ising) on a train/test split (if 'maxent' toolbox is available).
% - Evaluate model fit (DKL, marginal comparisons).
% - Run CCA between demographics and symptom features across subjects.
% - the code uses the Maxnet Toolbox to estimate maximum entropy models https://github.com/orimaoz/maxent_toolbox

%%
clear
close all
clc
rng(7);  % reproducible

%% ------------------------- Synthetic Data Generation Params ---------------------------------------
ncells         = 15;     % number of fMRI regions
nSubjects      = 120;    % number of synthetic subjects
samplesPerSubj = 300;    % timepoints per subject
burnin         = 200;    % Gibbs burn-in
thin           = 3;      % keep every 'thin'-th sample

% Ground-truth coupling graph density and scales
edgeDensity    = 0.15;   % sparsity of couplings
J_scale        = 0.8;    % coupling magnitude
h_scale        = 0.6;    % base field magnitude

% Demographics
age_mean = 40;  
age_sd = 12;
gender_p = 0.5;
edu_mean = 16;  
edu_sd = 2;

% Effect size knob (how strongly demographics modulate symptoms in Case A)
effect_scale = 1.0;  % try 0.5, 1.0, 1.5, 2.0 to boost the CCA difference

% Inter-subject noise on fields
demog_noise_sd = 0.15;


% How many subjects to visualize as time series examples
nShow = 3;           % how many subjects to show per case
nfMRIToShow = 6; % how many fMRI rows to display


%% ---------------------- DEFINE EFFECT VECTORS ----------------------------
% Per-symptom sensitivities (these are scaled by 'effect_scale' in Case A)
a_age = 0.35*randn(ncells,1);   % age effect
b_g   = 0.50*randn(ncells,1);   % gender effect
c_edu = 0.25*randn(ncells,1);   % education effect

%% ---------------------- CASE A: WITH DEMO EFFECTS -----------------------
cohortA = simulate_cohort(ncells, nSubjects, samplesPerSubj, ...
    edgeDensity, J_scale, h_scale, ...
    age_mean, age_sd, gender_p, edu_mean, edu_sd, ...
    a_age, b_g, c_edu, effect_scale, demog_noise_sd, ...
    burnin, thin);

resA = run_cca_and_report(cohortA.X_demo, cohortA.Y_means, cohortA.Y_sync);

%% ---------------------- CASE B: NO DEMO EFFECTS -------------------------
% Same everything, but effect_scale=0 (demographics do not influence symptoms)
cohortB = simulate_cohort(ncells, nSubjects, samplesPerSubj, ...
    edgeDensity, J_scale, h_scale, ...
    age_mean, age_sd, gender_p, edu_mean, edu_sd, ...
    a_age, b_g, c_edu, 0.0, demog_noise_sd, ...
    burnin, thin);

resB = run_cca_and_report(cohortB.X_demo, cohortB.Y_means, cohortB.Y_sync);

%% ---------------------- PRINT SUMMARIES ---------------------------------
fprintf('\n=== Canonical correlations (Case A: WITH effects, effect_scale=%.2f) ===\n', effect_scale);
disp(resA.r(1:min(6,end))');
if isfield(resA.stats,'pF')
    fprintf('Wilks’ Lambda p-values (sequential):\n');
    disp(resA.stats.pF(1:min(6,end))');
end

fprintf('\n=== Canonical correlations (Case B: NO effects) ===\n');
disp(resB.r(1:min(6,end))');
if isfield(resB.stats,'pF')
    fprintf('Wilks’ Lambda p-values (sequential):\n');
    disp(resB.stats.pF(1:min(6,end))');
end

%% ---------------------- MAXENT FIT & DIAGNOSTICS (A vs B) ----------------

    fprintf('\n[MaxEnt] Fitting/diagnosing Case A (WITH effects)\n');
    diagA = maxent_fit_eval(cohortA.TS, ncells, 'Case A (with effects)');

    fprintf('\n[MaxEnt] Fitting/diagnosing Case B (NO effects)\n');
    diagB = maxent_fit_eval(cohortB.TS, ncells, 'Case B (no effects)');


%% ---------------------- VISUAL COMPARISON PLOTS -------------------------
% Compare canonical correlations
kA = numel(resA.r); kB = numel(resB.r);
k  = max(kA,kB);
figure; hold on;
stem(1:kA, resA.r, 'filled'); 
stem((1:kB)+0.15, resB.r, 'filled'); 
xlabel('Canonical dimension'); ylabel('Correlation');
title(sprintf('CCA comparison: WITH effects (A) vs NO effects (B); effect\\_scale=%.2f', effect_scale));
legend('Case A: with effects','Case B: no effects','Location','best'); grid on;

% Show loadings on the first canonical variate
figure;
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

nexttile;
bar(resA.Xload(:,1)); 
set(gca,'XTick',1:numel(resA.Xnames),'XTickLabel',resA.Xnames,'XTickLabelRotation',30);
ylabel('corr(X var, U_1)'); title('Case A: X loadings on CV1');

nexttile;
bar(resB.Xload(:,1)); 
set(gca,'XTick',1:numel(resB.Xnames),'XTickLabel',resB.Xnames,'XTickLabelRotation',30);
ylabel('corr(X var, U_1)'); title('Case B: X loadings on CV1');

% Top Y loadings on CV1
figure;
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% Case A
nexttile;
[~,idxA] = sort(abs(resA.Yload(:,1)),'descend');
kshow = min(12, numel(idxA));
bar(resA.Yload(idxA(1:kshow),1));
set(gca,'XTick',1:kshow,'XTickLabel',resA.Ynames(idxA(1:kshow)),'XTickLabelRotation',30);
ylabel('corr(Y var, V_1)'); title('Case A: top Y loadings on CV1');

% Case B
nexttile;
[~,idxB] = sort(abs(resB.Yload(:,1)),'descend');
kshow = min(12, numel(idxB));
bar(resB.Yload(idxB(1:kshow),1));
set(gca,'XTick',1:kshow,'XTickLabel',resB.Ynames(idxB(1:kshow)),'XTickLabelRotation',30);
ylabel('corr(Y var, V_1)'); title('Case B: top Y loadings on CV1');

%% ---------------------- EXAMPLE TIME SERIES (RASTERS) -------------------
% Pick nShow random subjects from each cohort and display first nSymptomsToShow
timepts = 1:samplesPerSubj;

% Case A rasters
exA = randperm(nSubjects, nShow);
figure; 
tiledlayout(nShow,1,'Padding','compact','TileSpacing','compact');
for i=1:nShow
    Xs = cohortA.TS{exA(i)};
    nexttile;
    imagesc(timepts, 1:nfMRIToShow, Xs(1:nfMRIToShow,:));
    colormap(gray);
    ylabel('Symptom');
    xlabel('Time');
    title(sprintf('Case A (with effects) - Subject %d', exA(i)));
end
sgtitle(sprintf('Binary symptom time series (first %d symptoms): Case A', nfMRIToShow));

% Case B rasters
exB = randperm(nSubjects, nShow);
figure; 
tiledlayout(nShow,1,'Padding','compact','TileSpacing','compact');
for i=1:nShow
    Xs = cohortB.TS{exB(i)};
    nexttile;
    imagesc(timepts, 1:nfMRIToShow, Xs(1:nfMRIToShow,:));
    colormap(gray);
    ylabel('Symptom');
    xlabel('Time');
    title(sprintf('Case B (no effects) - Subject %d', exB(i)));
end
sgtitle(sprintf('Binary symptom time series (first %d symptoms): Case B', nfMRIToShow));


%% ---------------------- HELPERS -----------------------------------------
function samples01 = sampleIsingGibbs(n, T, h, J, burnin, thin)
    % Gibbs sampler for Ising model with spins in {-1,+1}, returns {0,1}
    totalIters = burnin + T*thin;
    s = sign(randn(n,1)); s(s==0)=1;
    kept = zeros(n, T);
    k = 0;
    for t = 1:totalIters
        for i=1:n
            localField = h(i) + J(i,:)*s;
            p_plus = 1/(1+exp(-2*localField));
            s(i) = (rand < p_plus)*2 - 1;
        end
        if t > burnin && mod(t-burnin, thin)==0
            k = k + 1;
            kept(:,k) = s;
            if k == T, break; end
        end
    end
    samples01 = (kept+1)/2;
end

function S = upperTriuVector(X)
    n = size(X,1);
    mask = triu(true(n),1);
    S = X(mask);
end

function [means, pairCorr] = empiricalStats01(samples01)
    n = size(samples01,1);
    T = size(samples01,2);
    means = mean(samples01, 2);
    Pij = (samples01 * samples01.')/T;  % P11
    pairCorr = upperTriuVector(Pij);
end

function J = makeSymmetricCouplings(n, density, scale)
    vals = randn(n);
    mask = triu(rand(n),1) < density;
    J = triu(scale * (vals .* mask),1);
    J = J + J.';
    J(1:n+1:end) = 0;
end

function cohort = simulate_cohort(ncells, nSubjects, samplesPerSubj, ...
        edgeDensity, J_scale, h_scale, ...
        age_mean, age_sd, gender_p, edu_mean, edu_sd, ...
        a_age, b_g, c_edu, effect_scale, demog_noise_sd, ...
        burnin, thin)

    % constant ground-truth base fields & couplings
    h0 = h_scale * randn(ncells,1);
    J0 = makeSymmetricCouplings(ncells, edgeDensity, J_scale);

    % storage
    X_demo = zeros(nSubjects,3);              % [age, gender, edu]
    Y_means = zeros(nSubjects, ncells);       % per-subject symptom means
    Y_sync  = zeros(nSubjects, 1);            % excess synchrony feature
    TS      = cell(nSubjects,1);              % time series per subject (ncells x T)

    for s = 1:nSubjects
        % demographics
        age    = age_mean + age_sd*randn;
        gender = double(rand < gender_p);
        edu    = max(6, edu_mean + edu_sd*randn);

        % z-scores for linear modulation
        age_z = (age - age_mean)/age_sd;
        edu_z = (edu - edu_mean)/edu_sd;

        % fields with/without demographic effect (scaled)
        h_subj = h0 + effect_scale*(a_age*age_z + b_g*gender + c_edu*edu_z) ...
                 + demog_noise_sd*randn(ncells,1);

        % sample
        Xs = sampleIsingGibbs(ncells, samplesPerSubj, h_subj, J0, burnin, thin);

        % features
        [m_i, ~] = empiricalStats01(Xs);
        Y_means(s,:) = m_i.';
        X_demo(s,:) = [age, gender, edu];
        TS{s} = Xs;
    end

    % pack
    cohort.X_demo = X_demo;
    cohort.Y_means = Y_means;
    cohort.Y_sync  = Y_sync;
    cohort.TS      = TS;
end

function out = run_cca_and_report(X_demo, Y_means, Y_sync)
    % Prepare X (z-score continuous; keep gender binary)
    age    = X_demo(:,1);
    gender = X_demo(:,2);
    edu    = X_demo(:,3);
    X = [zscore(age), gender, zscore(edu)];


        Y = Y_means;


    ok = all(isfinite(X),2) & all(isfinite(Y),2);
    X = X(ok,:); Y = Y(ok,:);

    [A,B,r,U,V,stats] = canoncorr(X,Y); %#ok<ASGLU>
    Xload = corr(X,U);
    Yload = corr(Y,V);

    out.r = r;
    out.stats = stats;
    out.Xload = Xload;
    out.Yload = Yload;
    out.Xnames = {'age_z','gender','edu_z'};
    featNames = arrayfun(@(i) sprintf('symptom_%02d_mean', i), 1:size(Y_means,2), 'UniformOutput', false);

    out.Ynames = featNames;
end

function diag = maxent_fit_eval(TS, ncells, plotLabel)
    % Pool all subjects' time series into one sample matrix (ncells x N)
    samples_all = [];
    for i = 1:numel(TS)
        Xi = TS{i};
        if size(Xi,1) ~= ncells
            error('TS{%d} has %d rows, expected %d', i, size(Xi,1), ncells);
        end
        samples_all = [samples_all, Xi]; %#ok<AGROW>
    end
    nsamples = size(samples_all,2);

    % Train/test split
    idx_train = randperm(nsamples, ceil(nsamples/2));
    idx_test  = setdiff(1:nsamples, idx_train);
    samples_train = samples_all(:, idx_train);
    samples_test  = samples_all(:, idx_test);

    % Build & train k-pairwise (Ising) model
    model = maxent.createModel(ncells,'kising');
    model = maxent.trainModel(model, samples_train, 'threshold', 1);

    % --- (1) KL divergence on held-out patterns
    empirical_distribution = maxent.getEmpiricalModel(samples_test);
    model_logprobs = maxent.getLogProbability(model, empirical_distribution.words);
    test_dkl = maxent.dkl(empirical_distribution.logprobs, model_logprobs);
    fprintf('Kullback-Leibler divergence (test): %f\n', test_dkl);

    % --- (2) Marginals: empirical vs predicted (log–log)
    marginals_data  = maxent.getEmpiricalMarginals(samples_test, model);
    marginals_model = maxent.getMarginals(model);

    figure;
    loglog(marginals_data, marginals_model, 'b*'); hold on;
    minval = min(marginals_data(marginals_data>0));
    plot([minval 1],[minval 1],'-r','LineWidth',1.25); % identity
    xlabel('empirical marginal'); ylabel('predicted marginal');
    title(sprintf('%s: marginals (%d cells)', plotLabel, ncells));
    grid on;

    % --- (3) Pattern probability scatter (log–log)
    emp_probs = exp(empirical_distribution.logprobs);
    mdl_probs = exp(model_logprobs);

    % keep only patterns observed in test (avoid zeros)
    keep = emp_probs > 0;
    emp_probs = emp_probs(keep);
    mdl_probs = mdl_probs(keep);

    figure;
    loglog(emp_probs, mdl_probs, 'k.'); hold on;
    mv = min(emp_probs(emp_probs>0));
    M  = max([emp_probs; mdl_probs]);
    plot([mv M],[mv M],'-r','LineWidth',1.25);
    xlabel('empirical pattern probability');
    ylabel('model pattern probability');
    title(sprintf('%s: pattern probabilities (test set)', plotLabel));
    grid on;

    % return diagnostics
    diag.dkl     = test_dkl;
    diag.H_emp   = empirical_distribution.entropy;
    diag.n_test_patterns = numel(emp_probs);
end
