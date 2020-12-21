csessionID = '200108';cfishID = '';
numT = 11;
%%
load(fullfile(getpath('behavior',csessionID,cfishID),'tail_swing'),'good_bout_idx','angledata','bout_idx','sum_curv');
load(fullfile(getpath('behavior',csessionID,cfishID),'align_with_fluo'));
load(fullfile(getpath('behavior',csessionID,cfishID),'high_analysis'),'conv_or_not');
load(fullfile(getpath('behavior',csessionID,cfishID),'low_video_analysis_result'),'no_param','param_head_angle_all','param_head_dist_all','param_speed_all');
load(fullfile(getpath('behavior',csessionID,cfishID),'leftright_eyes_to_param_angle_dist'),'left_eye_to_param_angle','right_eye_to_param_angle','left_eye_to_param_angle_visible', 'right_eye_to_param_angle_visible');
load(fullfile(getpath('neural activity',csessionID,cfishID),'spikes_new'),'spikes_OASIS');
load(fullfile(getpath('neural activity',csessionID,cfishID),'Coherence3_correct'),'center');
spikes_OASIS(:,1651:1779) = [];
Spike_X_EstTrace = smoothdata(spikes_OASIS,2,"gaussian",5);
%0108 good bout idx was already aligned to fluo data because its high and
%low fluo is very different from each other
angledata = angledata(align_with_fluo_low==1,:);
angledata_fluo = angledata(1:5:length(angledata),:);
[sum_curv,param_head_angle_all,param_head_dist_all,param_speed_all,no_param] = samfnmultvar(@(x) x(align_with_fluo_low==1),sum_curv,param_head_angle_all,param_head_dist_all,param_speed_all,no_param);
[left_eye_to_param_angle,right_eye_to_param_angle,left_eye_to_param_angle_visible,right_eye_to_param_angle_visible] = ...
    samfnmultvar(@(x) x(align_with_fluo_low==1),left_eye_to_param_angle,right_eye_to_param_angle,left_eye_to_param_angle_visible,right_eye_to_param_angle_visible);
[sum_curv_fluo,param_head_angle_fluo,param_head_dist_fluo,param_speed_fluo,no_param_fluo] = samfnmultvar(@(x) x(1:5:length(param_head_angle_all)),sum_curv,param_head_angle_all,param_head_dist_all,param_speed_all,no_param);
[lefteye_angle_fluo,righteye_angle_fluo,lefteye_angle_vis_fluo,righteye_angle_vis_fluo] = samfnmultvar(@(x) x(1:5:length(param_head_angle_all)),...
    left_eye_to_param_angle,right_eye_to_param_angle,left_eye_to_param_angle_visible,right_eye_to_param_angle_visible);
param_head_angle_hunting = mean(reshape(param_head_angle_all(reshape((good_bout_idx(:,1)-[1:20])',[],1)),20,[]),1);
mask = 1:5:nnz(align_with_fluo_low);
[~,hunting_idx_start] = min(abs(row2col(good_bout_idx(:,1),1) - col2row(mask)),[],2);
[~,hunting_idx_end] = min(abs(row2col(good_bout_idx(:,2),1) - col2row(mask)),[],2);
% fluo_idx = randperm(size(Spike_X_EstTrace,2)-numT,nnz(lia));
leftIdx = param_head_angle_hunting>0;
rightIdx = param_head_angle_hunting<0;
[numRegion,numWholeTime] = size(Spike_X_EstTrace);
numTrial = max(nnz(leftIdx),nnz(rightIdx));
if nnz(leftIdx)>=nnz(rightIdx)
    moreIdx = leftIdx;
    lessIdx = rightIdx;
    disp('left right');
    direction = {'left','right'};
else
    moreIdx = rightIdx;
    lessIdx = leftIdx;
    disp('right left');
    direction = {'right','left'};
end
conv_or_not = conv_or_not(align_with_fluo_high==1);
conv_or_not_fluo = conv_or_not(mask);
%%
%hunting bouts and nonhunting bouts
hunting_idx_wholecourse = arrayfun(@(i) (hunting_idx_start(i):hunting_idx_end(i))',1:length(hunting_idx_start),'un',0);
hunting_idx_wholecourse = cat(1,hunting_idx_wholecourse{:});
[~,bout_idx_aligned_start] = min(abs(row2col(bout_idx(:,1),1) - mask),[],2);
[~,bout_idx_aligned_end] = min(abs(row2col(bout_idx(:,2),1) - mask),[],2);
bout_idx_aligned_wholecourse = arrayfun(@(i) (bout_idx_aligned_start(i):bout_idx_aligned_end(i))',1:length(bout_idx),'un',0);
bout_idx_aligned_wholecourse = cat(1,bout_idx_aligned_wholecourse{:});
nonbout_idx = setdiff(1:numWholeTime,bout_idx_aligned_wholecourse);
simul_bout_idx_wholecourse = setdiff(bout_idx_aligned_wholecourse,hunting_idx_wholecourse);
rest_idx_wholecourse = setdiff(1:numWholeTime,bout_idx_aligned_wholecourse);
%%
visuo_idx = find(removeStrangeVisuoData(no_param_fluo,param_head_dist_fluo,param_speed_fluo));
hunting_idx_visuo = ismember(hunting_idx_wholecourse,visuo_idx);
hunting_idx_visuo = hunting_idx_wholecourse(hunting_idx_visuo);
nonhunting_idx = setdiff(1:numWholeTime,hunting_idx_wholecourse);
nonhunting_idx_visuo = ismember(nonhunting_idx,visuo_idx);
nonhunting_idx_visuo = nonhunting_idx(nonhunting_idx_visuo);
simul_idx_visuo = ismember(simul_bout_idx_wholecourse,visuo_idx);
simul_idx_visuo = simul_bout_idx_wholecourse(simul_idx_visuo);
rest_idx_visuo = ismember(rest_idx_wholecourse,visuo_idx);
rest_idx_visuo = rest_idx_wholecourse(rest_idx_visuo);
