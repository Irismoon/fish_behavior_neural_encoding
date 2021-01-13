%%
sessionID = '201117';fishID = '1';
%this script is to find out the effect of visual feedback on subsequent
%behavior. Specifically, how the positive and negative visual feedback
%affects the subsequent behavior
%%
%1.only when eye converged and classify the effect of first bout as
%positive or negative
load(fullfile(getpath('behavior',sessionID,fishID),'tail_swing'),'good_bout_idx','bout_idx','sum_curv');
load(fullfile(getpath('behavior',sessionID,fishID),'low_video_analysis_result'),'param_head_angle_all','param_head_dist_all');
load(fullfile(getpath('behavior',sessionID,fishID),'high_analysis'),'conv_or_not');
align_high_and_low;
idx = row2col(find(diff(conv_or_not)),1);
idx = idx(1:floor(length(idx)/2)*2);
conv_idx = reshape(idx,2,[]);
conv_duration = diff(conv_idx,1,1);
%the first conv bout->its conv phase->other bouts belonging to the conv
%phase
first_conv_bout = find(good_bout_idx(:,4)==1);
[~,index] = min(abs(row2col(good_bout_idx(first_conv_bout,1),1) - conv_idx(1,:)),[],2);%the conv idx corresponding to 1st conv bout
otherbouts = arrayfun(@(i) row2col(find(bout_idx(:,1)>conv_idx(1,index(i))+5 & bout_idx(:,1)<conv_idx(2,index(i))),1),1:length(first_conv_bout),'un',0);%other bouts belonging to this conv phase
%2.identify if the 1st conv bout is a positive or negative feedback,
%comparing the angle before and after the bout
angle = arrayfun(@(i) [mean(abs(param_head_angle_all(good_bout_idx(i,1)+[-3:0]))) mean(abs(param_head_angle_all(good_bout_idx(i,2)+[0:3])))],first_conv_bout,'un',0);
angle = cat(1,angle{:});
positive_or_not = diff(angle,1,2)<0;
figure,%compare the convergence duration of positive and negative bouts
boxplot(conv_duration(index),positive_or_not);set(gca,'XTickLabel',{'negative','positive'});
%2.the subsequent behavior coupling with visual stimuli
stim_pos_good = arrayfun(@(i) max(param_head_angle_all(good_bout_idx(i,1)+[-5:0])),1:size(good_bout_idx,1));
motor_amp_good = arrayfun(@(i) max(sum_curv(good_bout_idx(i,1):good_bout_idx(i,2))),1:size(good_bout_idx,1));
figure,
scatter(stim_pos_good,motor_amp_good);
hold on;
scatter(stim_pos_good(first_conv_bout),motor_amp_good(first_conv_bout));
title(sessionID);
%3.compare the bouts with positive or negative feedback
stim_pos = arrayfun(@(i) max(param_head_angle_all(max([bout_idx(i,1)+[-5:0];ones(1,6)]))),1:size(bout_idx,1));
motor_amp = arrayfun(@(i) mean(sum_curv(bout_idx(i,1):bout_idx(i,2))),1:size(bout_idx,1));
pos_otherbouts = otherbouts(positive_or_not);
neg_otherbouts = otherbouts(~positive_or_not);
pos_otherbouts = cat(1,pos_otherbouts{:});
neg_otherbouts = cat(1,neg_otherbouts{:});
figure,
scatter(stim_pos(pos_otherbouts),motor_amp(pos_otherbouts));
hold on;
scatter(stim_pos(neg_otherbouts),motor_amp(neg_otherbouts));
ax = gca;
xline(0,'--k');yline(0,'--k');
xlabel('param head angle');ylabel('mean(sum curv)');
legend({'positive feedback','negative feedback'});
title(sessionID);
%4.compare the convergence duration between +/- feedback