function fluo_idx = align_behavior2neuro(behavior_idx,behavior_len)
mask = 1:5:behavior_len;
[~,fluo_idx] = min(abs(row2col(behavior_idx,1) - col2row(mask)),[],2);
end

