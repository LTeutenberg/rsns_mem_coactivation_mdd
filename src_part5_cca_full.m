function [A_full, B_full, r_full, stats_full, X, Y, validIdx] = part5_cca_full(Prob_states, behavioral_data, C)
% 3. Canonical Correlation Analysis (CCA) on the whole sample

if ~exist(C.paths.out_cca,'dir'); mkdir(C.paths.out_cca); end
ta = @(col) table2array(behavioral_data(:,col));

HAMD5_item   = nanmean(ta(44:45), 2);
HAMD17_items = [ta([39:40 42:43 50:52 54 56:57 59:64]) HAMD5_item];
  
Y = HAMD17_items;
Y = fillmissing(Y, "nearest", 1);
X = [Prob_states];

validIdx = ~isnan(sum(X,2)) & ~isnan(sum(Y,2));
X = X(validIdx,:);  Y = Y(validIdx,:);

[A_full, B_full, r_full, ~, ~, stats_full] = canoncorr(X, Y);

save(fullfile(C.paths.out_cca,'cca_full.mat'), 'A_full','B_full','r_full','stats_full');
fprintf('Full-sample canonical correlations (first two): %.3f, %.3f | p=%.3g, %.3g\n', ...
    r_full(1), r_full(2), stats_full.p(1), stats_full.p(2));
end
