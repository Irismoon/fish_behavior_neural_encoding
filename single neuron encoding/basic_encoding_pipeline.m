function [outputArg1,outputArg2] = basic_encoding_pipeline(csessionID,cfishID)
%%
%loading
load(fullfile(getpath('behavior',csessionID,cfishID),'tail_swing'),'good_bout_idx');
load(fullfile(getpath('behavior',csessionID,cfishID),'align_with_fluo'),'align_with_fluo_low');
load(fullfile(getpath('behavior',csessionID,cfishID),'low_video_analysis_result'),'param_head_angle_all','param_head_dist_all','no_param');
load(fullfile(getpath('behavior',csessionID,cfishID),'leftright_eyes_to_param_angle_dist'),'left_eye_to_param_angle','right_eye_to_param_angle');
if strcmp(csessionID,'200108')
    load(fullfile(getpath('neural activity',csessionID,cfishID),'spikes_new'),'spikes_OASIS');
    spikes_OASIS(:,1651:1779) = [];
elseif strcmp(csessionID,'200703')
    load(fullfile(getpath('neural activity',cSession,cFish),'spikes.mat'), 'spikes_OASIS');
end 
spike_raw = smoothdata(spikes_OASIS,2,"gaussian",5);
%%
%data preprocessing
%no param
idx_low_keep = no_param==1
%too far
%not moving
%then align with fluo
%%
%encoding for visual stimulus

%encoding motor direction
end

