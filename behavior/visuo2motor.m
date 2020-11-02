function visuo2motor(sessionID,fishID)
%%
load(fullfile(getpath('behavior',sessionID,fishID),'tail_swing'),'sum_curv','bout_idx');
load(fullfile(getpath('behavior',sessionID,fishID),'high_analysis'),'conv_or_not');
load(fullfile(getpath('behavior',sessionID,fishID),'low_video_analysis_result'),'param_head_angle_all','no_param','param_head_dist_all','param_pos_all');
% load(fullfile(getpath('behavior',sessionID,fishID),'leftright_eyes_to_param_angle'));
if (length(conv_or_not)-length(sum_curv))>-5 && (length(conv_or_not)-length(sum_curv))<=0
    conv_or_not = [conv_or_not;zeros(length(sum_curv)-length(conv_or_not),1)];
elseif (length(conv_or_not)-length(sum_curv))<5 && (length(conv_or_not)-length(sum_curv))>0
    conv_or_not(end:end-(length(conv_or_not) - length(sum_curv))+1) = [];
else
    load(fullfile(getpath('behavior',sessionID,fishID),'align_with_fluo'));
    conv_or_not = conv_or_not(align_with_fluo_high==1);
    [sum_curv,param_head_angle_all,param_head_dist_all,param_pos_all] = samfnmultvar(@(x) x(align_with_fluo_low==1,:,:),sum_curv,param_head_angle_all,param_head_dist_all,param_pos_all);
end
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
conv_or_not_move = arrayfun(@(i) mean(conv_or_not(startFrame(i):startFrame(i)+5)),1:length(startFrame))>0;
%%
% mask = ~isoutlier(param_head_dist_move);
% [startFrame,param_head_angle_move,param_head_dist_move,sum_curv_move,Imaxabs] = samfnmultvar(@(x) x(mask),startFrame,param_head_angle_move,param_head_dist_move,sum_curv_move,Imaxabs);
%%
%compare bouts with to without eye convergence
figure,
[sum_curv_move_plt,param_head_angle_move_plt] = samfnmultvar(@(x) x(conv_or_not_move==true),sum_curv_move,param_head_angle_move);
subplot(1,2,1)
try
p = ranksum(sum_curv_move(param_head_angle_move_plt>0),sum_curv_move(param_head_angle_move_plt<=0));
catch ME
    return;
end
boxplot(sum_curv_move_plt,param_head_angle_move_plt>0);
sigstar([1 2],p);
title(['converged ' num2str(nnz((sum_curv_move_plt.*param_head_angle_move_plt)>0)/length(sum_curv_move_plt))]);
[sum_curv_move_plt,param_head_angle_move_plt] = samfnmultvar(@(x) x(conv_or_not_move==false),sum_curv_move,param_head_angle_move);
subplot(1,2,2)
try
p = ranksum(sum_curv_move(param_head_angle_move_plt>0),sum_curv_move(param_head_angle_move_plt<=0));
catch ME
    return;
end
boxplot(sum_curv_move_plt,param_head_angle_move_plt>0);
sigstar([1 2],p);
title(['unconverged ' num2str(nnz((sum_curv_move_plt.*param_head_angle_move_plt)>0)/length(sum_curv_move_plt))]);
sgtitle([sessionID ' fish ' fishID]);
%%
% label_prey_angle = param_head_angle_move>0;
% label_tail_direction = (sum_curv_move>0)*2;
% label = label_prey_angle+label_tail_direction;
% transition = [row2col(label(1:end-1),1) row2col(label(2:end),1)];
% trans_prob = zeros(4,4);
% for i=1:4
%     for j=1:4
%         trans_prob(i,j) = nnz(transition(:,1)==i-1&transition(:,2)==j-1)/length(transition);
%     end
% end
% figure('Position',[1994 517 1694 420]),
% subplot(1,4,1)
% imagesc(trans_prob);colorbar;
% set(gca,'XTick',1:4,'XTickLabel',{'RR','LR','RL','LL'},'YTick',1:4,'YTickLabel',{'RR','LR','RL','LL'});
% tail_trans = [row2col(label_tail_direction(1:end-1),1) row2col(label_tail_direction(2:end),1)]==0;%1 is right, 0 is left
% tail_trans_prob = zeros(2,2);
% tail_trans_prob(1,1) = nnz(tail_trans(:,1)==0&tail_trans(:,2)==0)/length(tail_trans);
% tail_trans_prob(1,2) = nnz(tail_trans(:,1)==0&tail_trans(:,2)==1)/length(tail_trans);
% tail_trans_prob(2,1) = nnz(tail_trans(:,1)==1&tail_trans(:,2)==0)/length(tail_trans);
% tail_trans_prob(2,2) = nnz(tail_trans(:,1)==1&tail_trans(:,2)==1)/length(tail_trans);
% subplot(1,4,2)
% imagesc(tail_trans_prob);title('tail transition');
% set(gca,'XTick',[1 2],'XTickLabel',{'L','R'},'YTick',[1 2],'YTickLabel',{'L','R'});
% prey_trans = [row2col(label_prey_angle(1:end-1),1) row2col(label_prey_angle(2:end),1)];%1 is left, 0 is right
% prey_trans_prob = zeros(2,2);
% for i=1:2
%     for j=1:2
%         prey_trans_prob(i,j) = nnz(prey_trans(:,1)==i-1&prey_trans(:,2)==j-1)/length(prey_trans);
%     end
% end
% prey_trans_prob = rot90(prey_trans_prob,2);
% subplot(1,4,3)
% imagesc(prey_trans_prob);
% set(gca,'XTick',[1 2],'XTickLabel',{'L','R'},'YTick',[1 2],'YTickLabel',{'L','R'});
% title('prey location transition');
% % trans_prob(i,i) = transition(:,1)==i-1&transition(:,2)==1;
% subplot(1,4,4)
% prop = tail_trans_prob-prey_trans_prob;
% imagesc(prop);colorbar;caxis([-1 1]);
% set(gca,'XTick',[1 2],'XTickLabel',{'L','R'},'YTick',[1 2],'YTickLabel',{'L','R'});
% title('tail/prey');
% sgtitle([sessionID ' ' fishID]);

%%
% figure('Position',[1927 430 1765 420]),
% subplot(1,2,1)
% param_head_angle_move_normal = angle_midline2normal(param_head_angle_move);
% polarscatter(param_head_angle_move_normal,(param_head_dist_move),10,sum_curv_move,'filled');
% hold on;
% [~,idx] = max(abs(sum_curv_move));
% polarscatter(param_head_angle_move_normal(idx),(param_head_dist_move(idx)),30,sum_curv_move(idx),'filled','MarkerEdgeColor','k');
% colormap('jet');cbar=colorbar;title(cbar,'sum_curv','Interpreter','none');
% cbar.Ticks = [-1,1],cbar.TickLabels={'right swing','left swing'};
% caxis([-max(abs(sum_curv_move)) max(abs(sum_curv_move))]);
% title('param to fish position just before move (radius log2 scale)');
%%
% subplot(1,3,2)
% reader = VideoReader(fullfile(getpath('behavior',sessionID,fishID),'low_unalign.avi'));
% frame = read(reader,startFrame(idx)+Imaxabs(idx)-1);
% imshow(frame);
% title('video frame at the largest point in left scatter');
%%
% subplot(1,2,2)
% boxplot(sum_curv_move,param_head_angle_move>0+1,'notch','on');
% p = ranksum(sum_curv_move(param_head_angle_move>0),sum_curv_move(param_head_angle_move<0));
% hold on;
% sigstar([1 2],p);
% line(get(gca,'XLim'),[0 0],'LineStyle','--','Color','k');
% ax = gca;
% ax.XTickLabel = {'at right','at left'};
% ax.YTick = [-1,1];
% ax.YTickLabel = {'right swing','left swing'};ax.YTickLabelRotation=90;
% ylabel('tail swing direction');xlabel('prey location');
% title('left swing:>0, right swing:<0');
% % param_head_angle_whole = angle_midline2normal(param_head_angle_all);
% % polarscatter(param_head_angle_whole,log2(param_head_dist_all),10,sum_curv,'filled');
% % colormap('cool');cbar=colorbar;title(cbar,'sum_curv','Interpreter','none');
% % title('param to fish position over all time (radius log2 scale)');
% sgtitle(sessionID);
% savefig(gcf,fullfile(getpath('behavior',sessionID,fishID),'tail_swing_direction_versus_prey_location'));
end
function [y,I] = maxabs(x)
[~,I] = max(abs(x));
y = x(I);
end