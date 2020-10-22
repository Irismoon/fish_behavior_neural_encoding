%%
load(fullfile(getpath('behavior',sessionID,fishID),'tail_swing_detect'),'bout_idx','tail_swing');
load(fullfile(getpath('behavior',sessionID,fishID),'low_video_analysis_result'),'param_head_angle_all_11210_chazhi','no_param','param_head_dist_all_11210','param_pos_all_11210','tail_amp_all');
% load(fullfile(getpath('behavior',sessionID,fishID),'leftright_eyes_to_param_angle_dist'));
%%
param_head_angle_all = param_head_angle_all_11210_chazhi;%this session record right head as negative and vice versa for left
param_head_dist_all = param_head_dist_all_11210;
%%
%check if when param is at different location, the induced behavior is
%different
% startFrame = detectMoveBlockStart(sessionID,fishID);
bout_idx = bout_idx(2:end-1);
[param_head_angle_move,param_head_dist_move] = samfnmultvar(@(x) arrayfun(@(i) mean(x(bout_idx(i)+(-20:-1))),1:length(bout_idx)),...
    param_head_angle_all,param_head_dist_all);
sum_curv_move = zeros(length(bout_idx),1);
for ibout = 1:length(bout_idx)
    tmp = tail_amp_all(bout_idx(ibout)+(1:5));
    [~,I] = max(abs(tmp));
    sum_curv_move(ibout) = tmp(I);
end
%%
numbout = length(bout_idx);
figure('Position',[2055 406 1087 362]),
subplot(1,2,1)
param_head_angle_move_normal = angle_midline2normal(param_head_angle_move);
polarscatter(param_head_angle_move_normal,log2(param_head_dist_move),10,sum_curv_move,'filled');
colormap('jet');cbar=colorbar;caxis([-max(abs(sum_curv_move)),max(abs(sum_curv_move))]);
title(cbar,'sum_curv','Interpreter','none');
title('param to fish position just before move (radius log2 scale)');
%%
subplot(1,2,2)
boxplot(sum_curv_move,param_head_angle_move>0,'notch','on');
p = ranksum(sum_curv_move(param_head_angle_move>0),sum_curv_move(param_head_angle_move<0));
hold on;
sigstar([1 2],p);
line(get(gca,'XLim'),[0 0],'LineStyle','--','Color','k');
ax = gca;
ax.XTickLabel = {'at left','at right'};
ylabel('tail swing direction');
title('')

savefig(gcf,fullfile(getpath('behavior',sessionID,fishID),'tail_swing_direction_versus_prey_location'));
