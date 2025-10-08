function results = part7_cca_generalizability(Prob_states, behavioral_data, C)
% 3.2 Generalizability — subgroup CCAs and greedy matching (no covariates)
% Changes vs previous draft: remove covariates entirely; keep original matching logic.

ta = @(col) table2array(behavioral_data(:,col));

group = ta(C.cols.group_status);
HC     = (group == 3);
cMDD   = (group == 0);
rMDD   = (group == 1) | (group == 2);
all_MDD = cMDD | rMDD;
all_subjects = HC | all_MDD;

groups = {all_subjects, HC, all_MDD, cMDD, rMDD};
names  = {'all_subjects','HC','all_MDD','cMDD','rMDD'};

% Precompute HAMD items once (full sample), then subset by G
HAMD5_item   = nanmean(ta(44:45),2);
HAMD17_items_all = [ta([39:40 42:43 50:52 54 56:57 59:64]) HAMD5_item];

for g = 1:numel(groups)
    G = groups{g};

    Y = fillmissing(HAMD17_items_all(G,:), "nearest", 1);
    X = Prob_states(G,:);

    valid = ~isnan(sum(X,2)) & ~isnan(sum(Y,2));
    [A, B, r, ~, ~, stats] = canoncorr(X(valid,:), Y(valid,:));
    A_std = zscore(A); B_std = zscore(B);

    results.(names{g}).r = r;
    results.(names{g}).p = stats.p;
    results.(names{g}).sigCV = find(stats.p < 0.05);
    results.(names{g}).A_loadings = arrayfun(@(k) A_std(1:128,k)', 1:length(r), 'uni', 0); % first 128 = brain states
    results.(names{g}).B_loadings = arrayfun(@(k) B_std(:,k)',       1:length(r), 'uni', 0);

    fprintf('Significant canonical variates for %s: %s\n', names{g}, mat2str(results.(names{g}).sigCV));
    for idx = results.(names{g}).sigCV(:)'
        fprintf('  CV %d: r = %.2f, p = %.3f\n', idx, r(idx), stats.p(idx));
    end
end

%% === Matching of subgroup CVs to full-sample CVs ===
if ~isfield(results, 'all_subjects') || isempty(results.all_subjects.sigCV)
    error('No significant CVs found for all_subjects. Cannot proceed with matching.');
end

gold_standard_CVs = results.all_subjects.sigCV;
nGold = length(gold_standard_CVs);

% Similarities among gold CVs
gold_similarities = NaN(nGold, nGold);
for i = 1:nGold
    gi = [results.all_subjects.A_loadings{gold_standard_CVs(i)} results.all_subjects.B_loadings{gold_standard_CVs(i)}];
    for j = 1:nGold
        gj = [results.all_subjects.A_loadings{gold_standard_CVs(j)} results.all_subjects.B_loadings{gold_standard_CVs(j)}];
        sim = abs(dot(gi, gj) / (norm(gi) * norm(gj)));
        gold_similarities(i, j) = sim;
    end
end

for g = 2:length(groups)  % Skip all_subjects
    group_name = names{g};
    if ~isfield(results, group_name)
        fprintf('Skipping %s: No CVs computed.\n', group_name);
        continue;
    end

    nGroup = length(results.(group_name).r);
    if nGroup == 0
        fprintf('Skipping %s: No CVs available.\n', group_name);
        continue;
    end

    % Compute similarity matrix using raw loadings
    similarity_matrix = NaN(nGold, nGroup);
    for i = 1:nGold
        gold_idx = gold_standard_CVs(i);
        gold_loadings = [results.all_subjects.A_loadings{gold_idx} results.all_subjects.B_loadings{gold_idx}];
        for j = 1:nGroup
            grp_loadings = [results.(group_name).A_loadings{j} results.(group_name).B_loadings{j}];
            similarity_matrix(i, j) = abs(dot(gold_loadings, grp_loadings) / (norm(gold_loadings) * norm(grp_loadings)));
        end
    end

    % Greedy matching: pick highest similarity for each gold CV
    results.(group_name).matched_CVs = NaN(nGold, 1);
    results.(group_name).matched_similarities = NaN(nGold, nGold);
    results.(group_name).matched_r = NaN(nGold, 1);
    results.(group_name).matched_p = NaN(nGold, 1);

    available_CVs = 1:nGroup;  % Track available CVs
    for i = 1:nGold
        if isempty(available_CVs)
            break;
        end
        [~, local_idx] = max(similarity_matrix(i, available_CVs));
        matched_idx = available_CVs(local_idx);
        results.(group_name).matched_CVs(i) = matched_idx;
        results.(group_name).matched_r(i) = results.(group_name).r(matched_idx);
        results.(group_name).matched_p(i) = results.(group_name).p(matched_idx);

        % Remove matched from availability
        available_CVs(local_idx) = [];

        % Similarities of this matched CV to all gold CVs
        grp_loadings = [results.(group_name).A_loadings{matched_idx} results.(group_name).B_loadings{matched_idx}];
        for k = 1:nGold
            gidx = gold_standard_CVs(k);
            gold_load = [results.all_subjects.A_loadings{gidx} results.all_subjects.B_loadings{gidx}];
            results.(group_name).matched_similarities(i, k) = abs(dot(gold_load, grp_loadings) / (norm(gold_load) * norm(grp_loadings)));
        end
    end

    % Similarities for significant CVs in the group
    sig_CVs = results.(group_name).sigCV;
    nSig = length(sig_CVs);
    results.(group_name).sig_similarities = NaN(nSig, nGold);
    for i = 1:nSig
        sig_idx = sig_CVs(i);
        sig_load = [results.(group_name).A_loadings{sig_idx} results.(group_name).B_loadings{sig_idx}];
        for k = 1:nGold
            gidx = gold_standard_CVs(k);
            gold_load = [results.all_subjects.A_loadings{gidx} results.all_subjects.B_loadings{gidx}];
            results.(group_name).sig_similarities(i, k) = abs(dot(gold_load, sig_load) / (norm(gold_load) * norm(sig_load)));
        end
    end

    % Display matching results
    fprintf('\nMatching Results for %s vs all_subjects (Greedy, Raw Loadings):\n', group_name);
    fprintf('Gold CV (all_subjects) | Matched CV (%s) | r-value | p-value | Similarities to all Gold CVs\n', group_name);
    for i = 1:nGold
        if ~isnan(results.(group_name).matched_CVs(i))
            fprintf('CV %d               | CV %d            | %.2f    | %.3f    |', ...
                gold_standard_CVs(i), results.(group_name).matched_CVs(i), ...
                results.(group_name).matched_r(i), results.(group_name).matched_p(i));
            for k = 1:nGold
                fprintf(' %.3f (CV %d)', results.(group_name).matched_similarities(i, k), gold_standard_CVs(k));
            end
            fprintf('\n');
        else
            fprintf('CV %d               | No match         | NaN     | NaN     | NaN\n', gold_standard_CVs(i));
        end
    end

    % Display similarities of significant CVs
    if ~isempty(sig_CVs)
        fprintf('\nSimilarities of Significant CVs in %s to all_subjects Significant CVs (Raw Loadings):\n', group_name);
        fprintf('Sig CV (%s) | r-value | p-value | Similarities to all Gold CVs\n', group_name);
        for i = 1:nSig
            fprintf('CV %d         | %.2f    | %.3f    |', sig_CVs(i), ...
                results.(group_name).r(sig_CVs(i)), results.(group_name).p(sig_CVs(i)));
            for k = 1:nGold
                fprintf(' %.3f (CV %d)', results.(group_name).sig_similarities(i, k), gold_standard_CVs(k));
            end
            fprintf('\n');
        end
    end
end

% Display gold standard similarities
fprintf('\nSimilarities Between Significant CVs in all_subjects (Raw Loadings):\n');
fprintf('CV | ');
for k = 1:nGold
    fprintf('CV %d ', gold_standard_CVs(k));
end
fprintf('\n');
for i = 1:nGold
    fprintf('CV %d | ', gold_standard_CVs(i));
    for k = 1:nGold
        fprintf('%.3f ', gold_similarities(i, k));
    end
    fprintf('\n');
end

fprintf('\nGold Standard (all_subjects) Significant CVs:\n');
fprintf('CV | r-value | p-value\n');
for k = 1:nGold
    fprintf('%d  | %.2f    | %.3f\n', gold_standard_CVs(k), results.all_subjects.r(gold_standard_CVs(k)), ...
        results.all_subjects.p(gold_standard_CVs(k)));
end
end
