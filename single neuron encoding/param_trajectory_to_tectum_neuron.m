load(fullfile(getpath('behavior',sessionID,fishID),'aligned_behavior'),'param_head_dist_all','param_head_angle_all','conv_or_not','left_tail_swing','right_tail_swing',...
    'left_eye_to_param_angle','right_eye_to_param_angle','param_lefteye_dist','param_righteye_dist','no_param');
load(fullfile(getpath('neural activity',sessionID,fishID),'Firing_Rate'));
%when the param is moving in front of the fish, how the neuron in tectum
%changes
%we observe that around tectum there is a burst activity at around 600
%frames, located in (362,262,144). Let's look at how it changes with param
%moving
%find the region it belongs to

%%
load(fullfile(getpath('neural activity',sessionID,fishID),'Firing_Rate'));
load(fullfile(getpath('neural activity',sessionID,fishID),'DFoF'));
load(fullfile(getpath('neural activity',sessionID,fishID),'Coherence3'),'A3');
%%
[startFrame,len] = detectMoveBlockStart(sessionID,fishID);
orig_idx = 1:len;
orig_idx = orig_idx(align_with_fluo_high==1);
startFrame = startFrame(ismember(startFrame,find(align_with_fluo_high))');
% [~,startFrame] = min(abs(startFrame - orig_idx(1:5:length(orig_idx))),[],2);
startFrame = arrayfun(@(i) find((startFrame(i)-orig_idx(1:5:length(orig_idx)))>0,1,'last'),1:length(startFrame)');
%%
tidx = 602;
dfof = reshape(A3*DFoF(:,tidx),600,600,280);
[~,page] = max(squeeze(dfof(396,279,:)));
idx = L_temp3(396,279,page);
%%
figure,hold on;
plot(zscore(Firing_Rate(idx,:)));
axis tight;
plot(zscore(param_head_angle_all));
plot(zscore(param_head_dist_all));
ylim = get(gca,'YLim');
arrayfun(@(i) line([startFrame(i) startFrame(i)],ylim,'Color','k','LineWidth',1,'LineStyle','--'),1:length(startFrame),'un',0);
legend({'firing rate','param head angle','param head dist'});
title(['region ' num2str(idx)]);
