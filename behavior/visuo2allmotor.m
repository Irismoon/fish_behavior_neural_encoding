function [outputArg1,outputArg2] = visuo2allmotor()
%comprehensive analysis of all behavior variables including tail, param,
%and eye
%%
load(fullfile(getpath('behavior',sessionID,fishID),'tail_swing'),'curvdata','sum_curv','left_tail_swing','right_tail_swing');
load(fullfile(getpath('behavior',sessionID,fishID),'low_video_analysis_result'),'param_head_angle_all','no_param','param_head_dist_all','param_pos_all');
load(fullfile(getpath('behavior',sessionID,fishID),'leftright_eyes_to_param_angle_dist'));
load(fullfile(getpath('behavior',sessionID,fishID),'align_with_fluo'));
load(fullfile(getpath('behavior',sessionID,fishID),'high_analysis'),'converge_angle');
%%
%lower dimension of curvdata, each PC is one pattern of tail swing, we
%could project one swing onto many different patterns and see which one
%it's mostly similar to.
[coeff,score,latent] = pca(curvdata);
%%
figure,plot(cumsum(latent/sum(latent)));axis tight
xlabel('# PC');ylabel('explained variance');
%%
mask_outlier = isoutlier(param_head_dist_all);
%%
%compare the similarity between the behavior data space and visual data
%space
behavior_space = [score(:,1:5) converge_angle];
[coeff2,score2,latent2] = pca(behavior_space);
visuospace = [param_head_dist_all.*cos(deg2rad(param_head_angle_all)) param_head_dist_all.*sin(deg2rad(param_head_angle_all))];
%%
figure,subplot(1,2,1),hold on;colormap('jet');
caxis = [1 length(score2)];
set(gca,'XLim',[min(score2(:,1)) max(score2(:,1))],'YLim',[min(score2(:,2)) max(score2(:,2))],'ZLim',[min(score2(:,3)) max(score2(:,3))]);
subplot(1,2,2),hold on;colormap('jet');
caxis = [1 length(score2)];grid on;
set(gca,'XLim',[min(visuospace(:,1)) max(visuospace(:,1))],'YLim',[min(visuospace(:,2)) max(visuospace(:,2))]);
subplot(1,2,1),
scatter3(score2(:,1),score2(:,2),score2(:,3),10,1:length(score2),'filled');
view(45,45);
subplot(1,2,2),
scatter(visuospace(:,1),visuospace(:,2),10,1:length(score2),'filled');
%%
figure,
subplot(2,2,1)
scatter3(score2(:,1),score2(:,2),score2(:,3),10,visuospace(:,1),'filled');colormap('jet');
title('x');
subplot(2,2,2)
scatter3(score2(:,1),score2(:,2),score2(:,3),10,visuospace(:,2),'filled');colormap('jet');
title('y');
subplot(2,2,3)
scatter3(score2(:,1),score2(:,2),score2(:,3),10,param_head_dist_all,'filled');colormap('jet');
title('dist');
subplot(2,2,4)
scatter3(score2(:,1),score2(:,2),score2(:,3),10,param_head_angle_all,'filled');colormap('jet');
title('angle');
%%
%cluster visuomotor space, this is to see if there exist clusters in the
%whole visuomotor space, which indicates that there exists different
%stim-response modes
visuomotor_space = [behavior_space visuospace];
idx = kmeans(normalize(visuomotor_space),4);
%%

end

