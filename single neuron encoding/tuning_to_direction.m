function tuning_to_direction()
%1.group the data segment when the param is located at the left and right
%of fish, (next step is finer angle tuning), 2.for each region, compare its
%activity in these two different conditions
%%
load(fullfile(getpath('behavior',sessionID,fishID),'aligned_behavior'),'param_head_dist_all','param_head_angle_all','conv_or_not','left_tail_swing','right_tail_swing',...
    'left_eye_to_param_angle','right_eye_to_param_angle','param_lefteye_dist','param_righteye_dist','no_param');
load(fullfile(getpath('neural activity',sessionID,fishID),'Firing_Rate'));
%%
%remove the initial invisible part
mask_have_param = no_param==1;
Firing_Rate = Firing_Rate(:,mask_have_param);
[param_head_angle_all,left_eye_to_param_angle,right_eye_to_param_angle,conv_or_not,left_tail_swing,right_tail_swing,param_head_dist_all,...
    move_dir_angle] = samfnmultvar(@(x) x(mask_have_param),param_head_angle_all,left_eye_to_param_angle,right_eye_to_param_angle,conv_or_not,...
    left_tail_swing,right_tail_swing,param_head_dist_all,move_dir_angle);
%%
%remove the data segment when fish is moving eye or tail
move_or_not = conv_or_not+abs(left_tail_swing)+abs(right_tail_swing) ~=1;
Firing_Rate = samfnmultvar(@(x) x(:,move_or_not),Firing_Rate);
[param_head_angle_all,left_eye_to_param_angle,right_eye_to_param_angle,param_head_dist_all,move_dir_angle] = ...
    samfnmultvar(@(x) x(move_or_not),param_head_angle_all,left_eye_to_param_angle,right_eye_to_param_angle,param_head_dist_all,...
    move_dir_angle);
%%
numRegion = size(Firing_Rate,1);
%%
%to head 
mask_head_left = param_head_angle_all>0;
mask_head_right = ~mask_head_left;
fr_head_left = Firing_Rate(:,mask_head_left);%
fr_head_right = Firing_Rate(:,mask_head_right);%region x time
%1.directly subtract from each other
fr_head_ratio = (mean(fr_head_left,2) - mean(fr_head_right,2))./(mean(fr_head_left,2) + mean(fr_head_right,2));
%%
%2.permutation test
p_head_ratio = zeros(numRegion,1);h_head_ratio = p_head_ratio;
for iR = 1:numRegion
    [p_head_ratio(iR),h_head_ratio(iR)] = shuffleTest(fr_head_left(iR,:),fr_head_right(iR,:));
end
%%
figure,
imshow(regionID2MIP(1/p_head_ratio,sessionID,fishID));
%%
%moving direction, approaching or away from

%%
%use PCs rather than raw trace
end

