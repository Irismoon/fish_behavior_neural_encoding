%align good bout idx to calcium
idx = ismember(good_bout_idx(:,1),find(align_with_fluo_low==1));
good_bout_idx = good_bout_idx(idx,:);
fluo_good_bout_idx = good_bout_idx;
idx = find(align_with_fluo_low);
idx = idx(1:5:end);
[~,fluo_good_bout_idx(:,1)] = min(abs(good_bout_idx(:,1) - col2row(idx,1)),[],2);
[~,fluo_good_bout_idx(:,2)] = min(abs(good_bout_idx(:,2) - col2row(idx,1)),[],2);

