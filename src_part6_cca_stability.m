function out = part6_cca_stability(X, Y, A_full, C)
% 3.1 CCA Stability analysis — 5-fold CV repeated 50 times (defaults)


kFolds   = C.cca.kfolds_stability;
nRepeats = C.cca.repeats_stability;
nSub = size(X,1);

r_matched = zeros(nRepeats * kFolds, 2);
cos_sim_matched = zeros(nRepeats * kFolds, 2);

cvCounter = 1;
for rep = 1:nRepeats
    % rng(rep); % optional reproducibility
    cv = cvpartition(nSub, 'KFold', kFolds);
    for fold = 1:kFolds
        tr = training(cv, fold);
        te = test(cv, fold);

        X_train = X(tr,:);  Y_train = Y(tr,:); 
        X_test  = X(te,:);  Y_test  = Y(te,:);

        [A_fold, ~, r_fold] = canoncorr(X_test, Y_test);

        % Cosine similarity vs A_full (first two)
        cs = zeros(size(A_fold,2), 2);
        for i = 1:size(A_fold,2)
            for j = 1:2
                cs(i,j) = abs(dot(A_fold(:,i), A_full(:,j)) / (norm(A_fold(:,i))*norm(A_full(:,j))));
            end
        end

        % Greedy match
        matchedIdx = zeros(2,1);
        avail = 1:size(A_fold,2);
        for j = 1:2
            [~, idx] = max(cs(avail,j));
            matchedIdx(j) = avail(idx);
            avail(idx) = [];
            if isempty(avail), break; end
        end

        r_matched(cvCounter,:) = r_fold(matchedIdx(1:2));
        cos_sim_matched(cvCounter,:) = [cs(matchedIdx(1),1), cs(matchedIdx(2),2)];
        cvCounter = cvCounter + 1;
    end
end

mean_cos = mean(cos_sim_matched, 1, 'omitnan');
std_cos  = std(cos_sim_matched, [], 1, 'omitnan');
fprintf('Stability — cosine similarity (mean±sd): CV1 %.3f±%.3f, CV2 %.3f±%.3f\n', ...
    mean_cos(1), std_cos(1), mean_cos(2), std_cos(2));

out.r_matched = r_matched;
out.cos_sim_matched = cos_sim_matched;
save(fullfile(C.paths.out_cca,'cca_stability.mat'), '-struct','out');
end
