function wholebrainmapping(spike_focus,param_head_angle_fluo,csessionID,cfishID,center,angledata,postfix)
% function wholebrainmapping(spike_focus,param_head_angle_fluo,csessionID,cfishID,center,angledata,postfix)
%angledata: tail x duration
%spike_focus, region x duration
duration = length(param_head_angle_fluo);
writer = VideoWriter(fullfile(getpath('result'),csessionID,['functional_mapping_frame_' postfix]));
open(writer);
fig = figure('Position',[2044 169 606 673]);
cmap = [min(spike_focus(:,:),[],'all') max(spike_focus(:,:),[],'all')];
for i=1:duration
%     ax = plot_brain_contour(center);
%     set(gca,ax(1));
    clf;
    sgtitle(['frame: ' num2str(i)]),
    ax2 = axes('Position',[[0.2000 0.3000 0.6000 0.5000]]);
    scatter(center(:,1),center(:,2),20,spike_focus(:,i),'filled');
    caxis(cmap);grid on;axis tight;colorbar;
    ax1 = polaraxes('Position',[0.1 0.1 0.8 0.7],'Color','none');
    %prey relative to head
    polarplot(pi/2-deg2rad(param_head_angle_fluo),4*ones(duration,1),'-k*');
    hold on;
    polarscatter(pi/2-deg2rad(param_head_angle_fluo(i)),4,20,'r','o','filled');
    %tail
    polarplot(pi+angledata(:,i),1+linspace(0,1,size(angledata,1)),'LineWidth',2);
    %two eyes
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
end

