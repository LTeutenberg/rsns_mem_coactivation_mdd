function T = part4_partialcorr(Prob_states, behavioral_data, C)
% 2.3 Partial Pearson correlations — severity

if ~exist(C.paths.out_corr,'dir'); mkdir(C.paths.out_corr); end
ta = @(col) table2array(behavioral_data(:,col));

site_variables = ta(C.cols.site_vars);
age    = ta(C.cols.age);
gender = ta(C.cols.gender);
Edu    = ta(C.cols.edu_years);
covariates = [age gender Edu site_variables];

severity = ta(C.cols.severity);

nstates = size(Prob_states,2);
pvals = NaN(nstates,1);
rvals = NaN(nstates,1);

valid_rows = ~isnan(severity) & ~isnan(sum(covariates,2));
Xcov = covariates(valid_rows,:);

for s = 1:nstates
    y = Prob_states(valid_rows, s);
    [r, p] = partialcorr(y, severity(valid_rows), Xcov);
    rvals(s) = r; pvals(s) = p;
end

if exist('fdr_bh','file')
    [h, crit_p, ~, adj_p] = fdr_bh(pvals, 0.05);
else
    [h, crit_p, ~, adj_p] = fdr_bh_local(pvals, 0.05);
end

T = table((1:nstates)', rvals, pvals, adj_p, h, crit_p, ...
    'VariableNames', {'state','r','p','p_fdr','h_sig','crit_p'});
writetable(T, fullfile(C.paths.out_corr,'partialcorr_results.csv'));
end
