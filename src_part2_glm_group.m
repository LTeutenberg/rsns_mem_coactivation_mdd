function T = part2_glm_group(Prob_states, behavioral_data, C)
% 2.1 General linear model — HC vs MDD

if ~exist(C.paths.out_glm,'dir'); mkdir(C.paths.out_glm); end
ta = @(col) table2array(behavioral_data(:,col));

site_variables = ta(C.cols.site_vars);
age    = ta(C.cols.age);
gender = ta(C.cols.gender);
Edu    = ta(C.cols.edu_years);
covariates = [age gender site_variables Edu];

group = ta(C.cols.group_HC_MDD);

nstates = size(Prob_states,2);
p_group = NaN(nstates,1);
t_group = NaN(nstates,1);

for s = 1:nstates
    y = Prob_states(:,s);
    mdl = fitglm([group, covariates], y);
    p_group(s) = mdl.Coefficients.pValue('x1');
    t_group(s) = mdl.Coefficients.tStat('x1');
end

if exist('fdr_bh','file')
    [h, crit_p, ~, adj_p] = fdr_bh(p_group, 0.05);
else
    [h, crit_p, ~, adj_p] = fdr_bh_local(p_group, 0.05);
end

T = table((1:nstates)', p_group, adj_p, t_group, h, crit_p, ...
    'VariableNames', {'state','p','p_fdr','t','h_sig','crit_p'});
writetable(T, fullfile(C.paths.out_glm,'HC_vs_MDD_results.csv'));
end
