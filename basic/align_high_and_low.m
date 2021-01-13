%align low and high
if abs(length(conv_or_not)-length(sum_curv))>5
    load(fullfile(getpath('behavior',sessionID,fishID),'align_with_fluo'));
    conv_or_not = conv_or_not(align_with_fluo_high==1);
    [sum_curv,param_head_angle_all,param_head_dist_all] = samfnmultvar(@(x) x(align_with_fluo_low==1),sum_curv,param_head_angle_all,param_head_dist_all);
end