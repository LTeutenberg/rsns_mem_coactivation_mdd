function T = part3_glm_status(Prob_states, behavioral_data, C)
% 2.2 General linear model — HC vs MDD acute/remitted (collapsed)

if ~exist(C.paths.out_glm,'dir'); mkdir(C.paths.out_glm); end
ta = @(col) table2array(behavioral_data(:,col));

site_variables = ta(C.cols.site_vars);
age    = ta(C.cols.age);
gender = ta(C.cols.gender);
Edu    = ta(C.cols.edu_years);
covariates = [age, gender, Edu, site_variables];

group = ta(C.cols.group_status);
group_new = zeros(size(group));
group_new(group==1 | group==2) = 1; % any MDD
group_new(group==3) = 2;            % keep mapping as in original

nstates = size(Prob_states,2);
p_status = NaN(nstates,1);

for s = 1:nstates
    y = Prob_states(:, s);
    mdl = fitglm([group_new, covariates], y);
    if ismember('x1', mdl.Coefficients.Properties.RowNames)
        p_status(s) = mdl.Coefficients.pValue('x1');
    end
end

if exist('fdr_bh','file')
    [h, crit_p, ~, adj_p] = fdr_bh(p_status, 0.05);
else
    [h, crit_p, ~, adj_p] = fdr_bh_local(p_status, 0.05);
end

T = table((1:nstates)', p_status, adj_p, h, crit_p, ...
    'VariableNames', {'state','p','p_fdr','h_sig','crit_p'});
writetable(T, fullfile(C.paths.out_glm,'Status_results.csv'));
end
