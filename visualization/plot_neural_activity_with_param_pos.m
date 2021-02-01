function plot_neural_activity_with_param_pos(spike_focus,center,param_head_angle_fluo,param_head_dist_fluo,angledata_fluo,name)
param_head_dist_fluo = row2col(mapminmax(col2row(param_head_dist_fluo),1,4),1);
writer = VideoWriter(fullfile(getpath('result'),name));
open(writer);
[N_x,bin_x] = discretize(center(:,1),[min(center(:,1)):5:(max(center(:,1))+5)]);
[N_y,bin_y] = discretize(center(:,2),[min(center(:,2)):5:(max(center(:,2))+5)]);
tmp = arrayfun(@(i) arrayfun(@(j) N_x==i&N_y==j,unique(N_y),'un',0),unique(N_x),'un',0);
tmp = cellfun(@(x) cat(2,x{:}),tmp,'un',0);
tmp = permute(cat(3,tmp{:}),[3 2 1]);%region x y x x
[row,col] = find(any(tmp,3));
bin_center = zeros(size(row,1),2);
bin_region = cell(size(bin_center,1),1);
for i=1:size(bin_center,1)
    bin_center(i,:) = [bin_x(row(i)) bin_y(col(i))];
    bin_region{i} = squeeze(tmp(row(i),col(i),:));
end

fig = figure('Position',[2044 169 606 673]);
cmap = [min(spike_focus,[],'all') quantile(spike_focus,.6,'all')];
for iT=1:size(spike_focus,2)
%     ax = plot_brain_contour(center);
%     set(gca,ax(1));
    clf;
    sgtitle(['frame: ' num2str(iT)]),
    ax2 = axes('Position',[[0.2000 0.3000 0.6000 0.5000]]);
    spike_max_bin = arrayfun(@(i) max(spike_focus(bin_region{i},iT)),1:size(bin_region,1));
    
    scatter(bin_center(:,1),bin_center(:,2),20,spike_max_bin,'filled');
    caxis(cmap);grid on;axis tight;view(90,90);
    ax1 = polaraxes('Position',[0.1 0.1 0.8 0.7],'Color','none');
    hold on;
    polarscatter(pi/2-deg2rad(param_head_angle_fluo(iT)),param_head_dist_fluo(iT),20,'r','o','filled');
    %tail
    polarplot(pi+angledata_fluo(iT,:),1+linspace(0,1,21),'LineWidth',2);
%     %two eyes
%     polarscatter(pi/2-[deg2rad(lefteye_angle_fluo(focus_idx_conv_start(idx)+i-6)) deg2rad(righteye_angle_fluo(focus_idx_conv_start(idx)+i-6))],[3 3],20,rgb({'green','blue'}),'o','filled');
%     polarscatter(pi/2-[deg2rad(left_eye_angle_fluo(focus_idx_conv_start(idx)+i-6)) deg2rad(right_eye_angle_fluo(focus_idx_conv_start(idx)+i-6))],[2.5 2.5],20,rgb({'green','blue'}),'o','filled');
%     title(['conv or not:' num2str(conv_or_not_fluo(focus_idx_conv_start(idx)+i-6))]);
    set(gca,'Color','none');
%     subplot(1,3,3),
%     scatter(center(:,2),center(:,3),20,spike_focus(:,i,idx),'filled');
%     caxis(cmap);grid on;
%     multiview;
%     pause(.5);
%     pause;
    frame = getframe(fig);
    writeVideo(writer,frame);
end
close(writer);
disp('run done！');
end

