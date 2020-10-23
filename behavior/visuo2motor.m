function visuo2motor(sessionID,fishID)
%%
load(fullfile(getpath('behavior',sessionID,fishID),'tail_swing'),'sum_curv','bout_idx');
load(fullfile(getpath('behavior',sessionID,fishID),'low_video_analysis_result'),'param_head_angle_all','no_param','param_head_dist_all','param_pos_all');
load(fullfile(getpath('behavior',sessionID,fishID),'leftright_eyes_to_param_angle'));
%%
%check if when param is at different location, the induced behavior is
%different
startFrame = bout_idx(:,1);

for iframe=1:length(startFrame)
    if startFrame(iframe)>20
        [param_head_angle_move(iframe),param_head_dist_move(iframe)] = samfnmultvar(@(x) mean(x(startFrame(iframe)-20:startFrame(iframe)-1)),...
            param_head_angle_all,param_head_dist_all);
    else
        disp('first frame detected!');
        [param_head_angle_move(iframe),param_head_dist_move(iframe)] = samfnmultvar(@(x) x(startFrame(iframe)),...
            param_head_angle_all,param_head_dist_all);
    end
end
[sum_curv_move,Imaxabs] = arrayfun(@(i) maxabs(sum_curv(startFrame(i):startFrame(i)+5)),1:length(startFrame));
%%
mask = ~isoutlier(param_head_dist_move);
[startFrame,param_head_angle_move,param_head_dist_move,sum_curv_move,Imaxabs] = samfnmultvar(@(x) x(mask),startFrame,param_head_angle_move,param_head_dist_move,sum_curv_move,Imaxabs);
%%
figure('Position',[1927 430 1765 420]),
subplot(1,2,1)
param_head_angle_move_normal = angle_midline2normal(param_head_angle_move);
polarscatter(param_head_angle_move_normal,(param_head_dist_move),10,sum_curv_move,'filled');
hold on;
[~,idx] = max(abs(sum_curv_move));
polarscatter(param_head_angle_move_normal(idx),(param_head_dist_move(idx)),30,sum_curv_move(idx),'filled','MarkerEdgeColor','k');
colormap('jet');cbar=colorbar;title(cbar,'sum_curv','Interpreter','none');
cbar.Ticks = [-1,1],cbar.TickLabels={'right swing','left swing'};
caxis([-max(abs(sum_curv_move)) max(abs(sum_curv_move))]);
title('param to fish position just before move (radius log2 scale)');
%%
% subplot(1,3,2)
% reader = VideoReader(fullfile(getpath('behavior',sessionID,fishID),'low_unalign.avi'));
% frame = read(reader,startFrame(idx)+Imaxabs(idx)-1);
% imshow(frame);
% title('video frame at the largest point in left scatter');
%%
subplot(1,2,2)
boxplot(sum_curv_move,param_head_angle_move>0+1,'notch','on');
p = ranksum(sum_curv_move(param_head_angle_move>0),sum_curv_move(param_head_angle_move<0));
hold on;
sigstar([1 2],p);
line(get(gca,'XLim'),[0 0],'LineStyle','--','Color','k');
ax = gca;
ax.XTickLabel = {'at right','at left'};
ax.YTick = [-1,1];
ax.YTickLabel = {'right swing','left swing'};ax.YTickLabelRotation=90;
ylabel('tail swing direction');xlabel('prey location');
title('left swing:>0, right swing:<0');
% param_head_angle_whole = angle_midline2normal(param_head_angle_all);
% polarscatter(param_head_angle_whole,log2(param_head_dist_all),10,sum_curv,'filled');
% colormap('cool');cbar=colorbar;title(cbar,'sum_curv','Interpreter','none');
% title('param to fish position over all time (radius log2 scale)');
sgtitle(sessionID);
savefig(gcf,fullfile(getpath('behavior',sessionID,fishID),'tail_swing_direction_versus_prey_location'));
end
function [y,I] = maxabs(x)
[~,I] = max(abs(x));
y = x(I);
end