function fluo_idx = align_behavior2neuro(bout_idx,align_with_fluo_low)
%function fluo_idx = align_behavior2neuro(bout_idx,align_with_fluo_low)
bout_label = arrayfun(@(i) all(align_with_fluo_low(bout_idx(i,1)):align_with_fluo_low(bout_idx(i,2))==1),1:size(bout_idx,1));
bout_idx = bout_idx(bout_label,:);
low_idx = find(align_with_fluo_low);
new_bout_idx = arrayfun(@(i) [find(low_idx==bout_idx(i,1)) find(low_idx==bout_idx(i,2))],1:size(bout_idx,1),'un',0);
new_bout_idx = cat(1,new_bout_idx{:});
bout_idx(:,1:2) = new_bout_idx;
%delete very short bout
bout_len = bout_idx(:,2) - bout_idx(:,1) + 1;
bout_idx(bout_len<4,:) = [];
%align to fluo index
mask = 1:5:nnz(align_with_fluo_low);
fluo_idx = bout_idx;
[~,fluo_idx(:,1)] = min(abs(row2col(bout_idx(:,1),1) - col2row(mask)),[],2);
[~,fluo_idx(:,2)] = min(abs(row2col(bout_idx(:,2),1) - col2row(mask)),[],2);
end

