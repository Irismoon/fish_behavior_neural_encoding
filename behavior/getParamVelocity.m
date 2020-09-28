function [outputArg1,outputArg2] = getParamVelocity()
load(fullfile(getpath('behavior',sessionID,fishID),'low_video_analysis_result'),'param_head_angle_all','param_head_dist_all');
param_head_angle_all = deg2rad(param_head_angle_all);
move_dir = param_head_dist_all.*(cos(param_head_angle_all)+1i*sin(param_head_angle_all));
move_dir = diff(move_dir);
move_dir = move_dir(align_with_fluo_high(2:end)==1);
move_dir = [0;move_dir(1:end-1)];%such that the second element of move_dir reflect the movement from 1 to 2, while the 2nd frame of calcium imaging 
%reflect the reaction to param movement from 1 to 2
param_head_velocity = move_dir(1:5:length(move_dir));
save(fullfile(getpath('behavior',sessionID,fishID),'aligned_behavior'),'param_head_velocity','-append');
end

