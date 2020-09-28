function visuo2motor(sessionID,fishID)
%%
load(fullfile(getpath('behavior',sessionID,fishID),'tail_swing'),'sum_curv','left_tail_swing','right_tail_swing');
load(fullfile(getpath('behavior',sessionID,fishID),'low_video_analysis_result'),'param_head_angle_all','no_param','param_head_dist_all','param_pos_all');
load(fullfile(getpath('behavior',sessionID,fishID),'leftright_eyes_to_param_angle_dist'));
load(fullfile(getpath('behavior',sessionID,fishID),'align_with_fluo'));
%%
%check if when param is at different location, the induced behavior is
%different
startFrame = detectMoveBlockStart(sessionID,fishID);
[param_head_angle_move,param_head_dist_move,left_eye_to_param_angle_move,right_eye_to_param_angle_move,...
    param_lefteye_dist_move,param_righteye_dist_move] = samfnmultvar(@(x) arrayfun(@(i) mean(x(startFrame(i)-20:startFrame(i)-1)),1:length(startFrame)),...
    param_head_angle_all,param_head_dist_all,left_eye_to_param_angle,right_eye_to_param_angle,...
    param_lefteye_dist,param_righteye_dist);
sum_curv_move = arrayfun(@(i) max(sum_curv(startFrame(i):startFrame(i)+5)),1:length(startFrame));
%%
mask = ~isoutlier(param_head_dist_move);
[param_head_angle_move,param_head_dist_move,sum_curv_move] = samfnmultvar(@(x) x(mask),param_head_angle_move,param_head_dist_move,sum_curv_move);
%%
figure,
param_head_angle_move_normal = angle_midline2normal(param_head_angle_move);
polarscatter(param_head_angle_move_normal,(param_head_dist_move),10,sum_curv_move,'filled');
colormap('jet');cbar=colorbar;title(cbar,'sum_curv','Interpreter','none');
caxis([-max(sum_curv_move) max(sum_curv_move)]);
title('param to fish position just before move (radius log2 scale)');
%%
figure,
param_head_angle_whole = angle_midline2normal(param_head_angle_all);
polarscatter(param_head_angle_whole,log2(param_head_dist_all),10,sum_curv,'filled');
colormap('cool');cbar=colorbar;title(cbar,'sum_curv','Interpreter','none');
title('param to fish position over all time (radius log2 scale)');
%%
figure,
histogram2Polar(param_head_angle_whole(mask),param_head_dist_all(mask),2,'ThetaZeroLocation','right','Normalization','probability');
end