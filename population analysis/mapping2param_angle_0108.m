idx = hunting_idx_visuo;
% idx = nonhunting_idx_visuo;
% idx = simul_idx_visuo;
% idx = rest_idx_visuo;
param_head_angle_visuo = param_head_angle_fluo(idx);
spikes_visuo = Spike_X_EstTrace(:,idx);
tuned2angle(param_head_angle_visuo,spikes_visuo,center,'single');
%%
%investigate 20 degree regional activity during hunting
focus_idx_conv_start = hunting_idx_visuo();
spike_focus = zeros(numRegion,15,length(focus_idx_conv_start));
for i=1:length(focus_idx_conv_start)
    spike_focus(:,:,i) = Spike_X_EstTrace(:,focus_idx_conv_start(i)+[-5:9]);
end
idx = 4;
figure('Position',[2044 169 1278 420]),ax(1) = subplot(1,3,1,polaraxes),ax(2) = subplot(1,3,2),ax(3) = subplot(1,3,3),
cmap = [min(spike_focus(:,:,idx),[],'all') max(spike_focus(:,:,idx),[],'all')];
for i=1:15
%     ax = plot_brain_contour(center);
%     set(gca,ax(1));
    sgtitle(num2str(i)),
    subplot(1,3,1),cla,
    polarplot(deg2rad(param_head_angle_fluo(focus_idx_conv_start(idx)+[-5:5])),ones(11,1),'-k*');
    polarplot(angledata_fluo(focus_idx_conv_start(idx)+i-6,:)+pi/2,linspace(0,1,21));
    hold on;
    polarscatter(deg2rad(param_head_angle_fluo(focus_idx_conv_start(idx)+i-6)),1,20,'r','o','filled');
    polarscatter([deg2rad(lefteye_angle_fluo(focus_idx_conv_start(idx)+i-6)) deg2rad(righteye_angle_fluo(focus_idx_conv_start(idx)+i-6))],[.5 .5],20,rgb({'green','blue'}),'o','filled');
    polarscatter([deg2rad(left_eye_angle_fluo(focus_idx_conv_start(idx)+i-6)) deg2rad(right_eye_angle_fluo(focus_idx_conv_start(idx)+i-6))],[.25 .25],20,rgb({'green','blue'}),'o','filled');
    title(num2str(conv_or_not_fluo(focus_idx_conv_start(idx)+i-6)));
    
    subplot(1,3,2),
    scatter(center(:,1),center(:,2),20,spike_focus(:,i,idx),'filled');
    caxis(cmap);grid on;

    subplot(1,3,3),
    scatter(center(:,2),center(:,3),20,spike_focus(:,i,idx),'filled');
    caxis(cmap);grid on;
%     multiview;
    pause(.5);
end
%%
%investigate 20 degree regional activity during hunting
disp('running this section....');
duration = 15;
focus_idx_conv_start = hunting_idx_visuo(60);
% focus_idx = simul_idx_visuo(267);
% focus_idx = rest_idx_visuo(140)
spike_focus = zeros(numRegion,duration,length(focus_idx_conv_start));
for i=1:length(focus_idx_conv_start)
    spike_focus(:,:,i) = Spike_X_EstTrace(:,focus_idx_conv_start(i)+[-5:duration-6]);
end
idx = 1;
writer = VideoWriter(fullfile(getpath('result'),csessionID,['functional_mapping_frame_' num2str(focus_idx_conv_start(idx))]));
open(writer);
fig = figure('Position',[2044 169 606 673]);
cmap = [min(spike_focus(:,:,idx),[],'all') max(spike_focus(:,:,idx),[],'all')];
for i=1:duration
%     ax = plot_brain_contour(center);
%     set(gca,ax(1));
    clf;
    sgtitle(['frame: ' num2str(i)]),
    ax2 = axes('Position',[[0.2000 0.3000 0.6000 0.5000]]);
    scatter(center(:,1),center(:,2),20,spike_focus(:,i,idx),'filled');
    caxis(cmap);grid on;axis tight;
    ax1 = polaraxes('Position',[0.1 0.1 0.8 0.7],'Color','none');
    %prey relative to head
    polarplot(pi/2-deg2rad(param_head_angle_fluo(focus_idx_conv_start(idx)+[-5:5])),4*ones(11,1),'-k*');
    hold on;
    polarscatter(pi/2-deg2rad(param_head_angle_fluo(focus_idx_conv_start(idx)+i-6)),4,20,'r','o','filled');
    %tail
    polarplot(pi+angledata_fluo(focus_idx_conv_start(idx)+i-6,:),1+linspace(0,1,21),'LineWidth',2);
    %two eyes
    polarscatter(pi/2-[deg2rad(lefteye_angle_fluo(focus_idx_conv_start(idx)+i-6)) deg2rad(righteye_angle_fluo(focus_idx_conv_start(idx)+i-6))],[3 3],20,rgb({'green','blue'}),'o','filled');
    polarscatter(pi/2-[deg2rad(left_eye_angle_fluo(focus_idx_conv_start(idx)+i-6)) deg2rad(right_eye_angle_fluo(focus_idx_conv_start(idx)+i-6))],[2.5 2.5],20,rgb({'green','blue'}),'o','filled');
    title(['conv or not:' num2str(conv_or_not_fluo(focus_idx_conv_start(idx)+i-6))]);
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
%%
%average across hunting
disp('running this section....');
duration = 15;
focus_idx_conv_start = hunting_idx_visuo(ismember(hunting_idx_visuo,hunting_idx_start));
left_param_mask = arrayfun(@(i) mean(param_head_angle_fluo(focus_idx_conv_start(i)+[-5:1])),1:length(focus_idx_conv_start))>0;
focus_idx_conv_start = focus_idx_conv_start(~left_param_mask);
spike_focus = zeros(numRegion,duration,length(focus_idx_conv_start));
for i=1:length(focus_idx_conv_start)
    spike_focus(:,:,i) = Spike_X_EstTrace(:,focus_idx_conv_start(i)+[-5:duration-6]);
end
spike_focus = mean(spike_focus,3);
idx = 1;
writer = VideoWriter(fullfile(getpath('result'),csessionID,['functional_mapping_frame_' num2str(focus_idx_conv_start(idx))]));
open(writer);
fig = figure('Position',[2044 169 606 673]);
cmap = [min(spike_focus(:,:,idx),[],'all') max(spike_focus(:,:,idx),[],'all')];
for i=1:duration
%     ax = plot_brain_contour(center);
%     set(gca,ax(1));
    clf;
    sgtitle(['frame: ' num2str(i)]),
    ax2 = axes('Position',[[0.2000 0.3000 0.6000 0.5000]]);
    scatter(center(:,1),center(:,2),20,spike_focus(:,i,idx),'filled');
    caxis(cmap);grid on;axis tight;
    ax1 = polaraxes('Position',[0.1 0.1 0.8 0.7],'Color','none');
    %prey relative to head
    polarplot(pi/2-deg2rad(param_head_angle_fluo(focus_idx_conv_start(idx)+[-5:5])),4*ones(11,1),'-k*');
    hold on;
    polarscatter(pi/2-deg2rad(param_head_angle_fluo(focus_idx_conv_start(idx)+i-6)),4,20,'r','o','filled');
    %tail
    polarplot(pi+angledata_fluo(focus_idx_conv_start(idx)+i-6,:),1+linspace(0,1,21),'LineWidth',2);
    %two eyes
    polarscatter(pi/2-[deg2rad(lefteye_angle_fluo(focus_idx_conv_start(idx)+i-6)) deg2rad(righteye_angle_fluo(focus_idx_conv_start(idx)+i-6))],[3 3],20,rgb({'green','blue'}),'o','filled');
    polarscatter(pi/2-[deg2rad(left_eye_angle_fluo(focus_idx_conv_start(idx)+i-6)) deg2rad(right_eye_angle_fluo(focus_idx_conv_start(idx)+i-6))],[2.5 2.5],20,rgb({'green','blue'}),'o','filled');
    title(['conv or not:' num2str(conv_or_not_fluo(focus_idx_conv_start(idx)+i-6))]);
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
%%
idx = idx & ~(lefteye_angle_vis_fluo==0) & ~(righteye_angle_vis_fluo==0);
lefteye_angle_fluo = lefteye_angle_fluo(idx);
spikes_visuo = Spike_X_EstTrace(:,idx);
tuned2angle(lefteye_angle_fluo,spikes_visuo,center,'single');
%%
duration = 15;
focus_idx_conv_start = hunting_idx_visuo(ismember(hunting_idx_visuo,hunting_idx_start));
left_param_mask = arrayfun(@(i) mean(param_head_angle_fluo(focus_idx_conv_start(i)+[-5:1])),1:length(focus_idx_conv_start))>0;
focus_idx_conv_start = focus_idx_conv_start(~left_param_mask);
spike_focus = zeros(numRegion,duration,length(focus_idx_conv_start));
for i=1:length(focus_idx_conv_start)
    spike_focus(:,:,i) = Spike_X_EstTrace(:,focus_idx_conv_start(i)+[-5:duration-6]);
end%ignore the actual difference of visual stimulus, like different angle, distance, direction, speed, etc.
%clustering the areas with very similar response 
%cross-correlation between regions to find out the info flow route
del_region = sum(all(spike_focus==0,2),3)>0.1*length(focus_idx_conv_start);
keep_region = find(~del_region);
spike_focus(del_region,:,:) = [];
spike_focus = reshape(spike_focus,size(spike_focus,1),[]);

orig_r = xcorr(spike_focus',3,'coeff');%(2 x duration + 1)xregion^2, earlier and later
[r,lag] = max(orig_r,[],1);%if lag<6, x is earlier than y; else, later
%use the r to find out the cluster
r = reshape(r,sqrt(length(r)),sqrt(length(r)));
r = tril(r) + tril(r,-1)';
r(index_diag(r)) = 1;
% idx = spectralcluster(r,20,'Distance','precomputed');
Z = linkage(squareform(1-r),'ward');
idx = cluster(Z,'maxclust',2);
thld = median([Z(end-5,3) Z(end-4,3)]);
dendrogram(Z,'ColorThreshold',thld);
cophenet(Z,squareform(1-r));
% color = brewermap(length(unique(idx)),'Set3');
N = histcounts(idx,[1:length(unique(idx))+1]-0.5);
color = colormap(parula(length(unique(idx))));
[N,I] = sort(N,'ascend');
color = color(I,:);
figure,hold on,
for i=unique(idx)'
    scatter3(center(keep_region(idx==i),1),center(keep_region(idx==i),2),center(keep_region(idx==i),3),10,color(i,:),'filled');
end
%%
%how the variability of neural activity around movement onset change with
%time
figure,
plot(mean(var(spike_focus,[],3)));ylabel('variance')
yyaxis right;
plot(mean(mean(spike_focus,3),1));ylabel('average firing rate');
figure,plot(mean(var(spike_focus,[],3))./mean(mean(spike_focus,3),1));
hold on;
plot(mean(var(spike_focus,[],3))./mean(mean(spike_focus,3),1));
focus_idx_conv_start = randperm(size(Spike_X_EstTrace,2)-duration,length(focus_idx_conv_start));
spike_focus = zeros(length(:),duration,length(focus_idx_conv_start));
for i=1:length(focus_idx_conv_start)
    spike_focus(:,:,i) = Spike_X_EstTrace(:,focus_idx_conv_start(i)+[-5:duration-6]);
end%ignore the actual difference of visual stimulus, like different angle, distance, direction, speed, etc.
%%
spike_focus = mean(spike_focus,3);
thld = quantile(spike_focus(:),.6);
detect_region = median(spike_focus(:,4:13),2)>thld;
[activation,act_time] = max(spike_focus(detect_region,4:13),[],2);
plot_brain_contour(center,'3d');
scatter3(center(keep_region(detect_region),1),center(keep_region(detect_region),2),center(keep_region(detect_region),3),10,act_time,'filled');
colormap('jet');
%%
%the single trials trace of focus region
roi_pos = [362 291 104;306 331 164;284 166 99;312 75 128;231 298 121;228 159 100];%
chnname = {'ot1','fb1','mb1','hb1','contra-ot','hb2'};
%ot1, fb1,mb1, hb1, contra-ot,hb2
[~,roi_idx] = arrayfun(@(i) min(sum((center -roi_pos(i,:)).^2,2)),1:size(roi_pos,1));
% multchnPlot(permute(mean(spike_focus(roi_idx,:,:),3),[2 3 1]),[-0.5,0.1*(duration-5)],[0],[],chnname,'plotstyle','plot');
for i=1:size(spike_focus,3)
    multchnPlot(permute(spike_focus(roi_idx,:,i),[2 3 1]),[-0.5,0.1*(duration-5)],[0],[],chnname,'plotstyle','plot');
end
%%
%plot an interesting segment
x_idx = 7;
idx = focus_idx_conv_start(x_idx)+[-5:30];
t = [-5:30]*0.1;
fig = multchnPlot(reshape(Spike_X_EstTrace(roi_idx,idx)',length(idx),1,[]),t,t(find(diff(conv_or_not_fluo(idx)))+1),[],chnname,'plotstyle','plot');
load(fullfile(getpath('behavior',csessionID,cfishID),'tail_swing'),'sum_curv');
sum_curv_fluo = sum_curv(align_with_fluo_low==1);
ax = flipud(fig.findobj('type','axes'));
for i=1:length(ax)
    set(fig,'currentaxes',ax(i));
    yyaxis right;
%     plot(linspace(t(1),t(end),length(mask(idx(1)):mask(idx(end)))),...
%         sum_curv_fluo(mask(idx(1)):mask(idx(end))),'k');
    plot(linspace(t(1),t(end),length(idx)),sum_curv_fluo(mask(idx)));
end
sgtitle(focus_idx_conv_start(x_idx));
%%
figure,
subplot(2,1,1),
mspk = mean(spike_focus(roi_idx,:,:),3)';%t x region
sspk = std(spike_focus(roi_idx,:,:),0,3)'/sqrt(size(spike_focus,3));
plot(-0.5:0.1:0.9,mspk);axis tight;
ylabel('mean firing rate');
xlabel('t/s');
legend(chnname);
% boundedline(-0.5:0.1:0.9,mspk,reshape(sspk,duration,1,length(roi_idx)),'alpha');
subplot(2,1,2),
plot(-0.5:0.1:0.9,sspk);axis tight;
ylabel('std(firing rate)');
xlabel('t/s');
% plot(mean(diff(spike_focus(roi_idx,:,:),[],2),3)'./mean(spike_focus(roi_idx,1:end-1,:),3)')
%%
disp('this section running...');
start = 3;
duration = 7;
[label,I] = ismember(hunting_idx_visuo,hunting_idx_start);
%separate convergence start and continuing bouts
focus_idx_conv = hunting_idx_start(I(label));
focus_idx_conv_start = focus_idx_conv(good_bout_idx(I(label),4)==1);%convergence start
focus_idx_conv_cont = focus_idx_conv(good_bout_idx(I(label),4)==0);%convergence
% continuing
left_param_mask = arrayfun(@(i) mean(param_head_angle_fluo(focus_idx_conv_start(i)+[-5:1])),1:length(focus_idx_conv_start))>0;
focus_idx_conv_start = focus_idx_conv_start(~left_param_mask);
left_param_mask = arrayfun(@(i) mean(param_head_angle_fluo(focus_idx_conv_cont(i)+[-5:1])),1:length(focus_idx_conv_cont))>0;
focus_idx_conv_cont = focus_idx_conv_cont(~left_param_mask);
spike_focus_conv_start = zeros(numRegion,duration,length(focus_idx_conv_start));
sum_curv_conv_start = zeros(duration,length(focus_idx_conv_start));
param_head_dist_conv_start = zeros(duration,length(focus_idx_conv_start));
param_head_angle_conv_start = param_head_dist_conv_start;
for i=1:length(focus_idx_conv_start)
    spike_focus_conv_start(:,:,i) = Spike_X_EstTrace(:,focus_idx_conv_start(i)+[-start:duration-start-1]) - mean(Spike_X_EstTrace(:,focus_idx_conv_start(i)-start-[5:-1;0]),2);
    sum_curv_conv_start(:,i) = sum_curv_fluo(focus_idx_conv_start(i)+[-start:duration-start-1]);
    param_head_dist_conv_start(:,i) = param_head_dist_fluo(focus_idx_conv_start(i)+[-start:duration-start-1]);
    param_head_angle_conv_start(:,i) = param_head_angle_fluo(focus_idx_conv_start(i)+[-start:duration-start-1]);
end
spike_focus_conv_cont = zeros(numRegion,duration,length(focus_idx_conv_cont));
sum_curv_conv_cont = zeros(duration,length(focus_idx_conv_cont));
param_head_dist_conv_cont = zeros(duration,length(focus_idx_conv_cont));
param_head_angle_conv_cont = param_head_dist_conv_cont;
for i=1:length(focus_idx_conv_cont)
    spike_focus_conv_cont(:,:,i) = Spike_X_EstTrace(:,focus_idx_conv_cont(i)+[-start:duration-start-1]) - Spike_X_EstTrace(:,focus_idx_conv_cont(i)-start);
    sum_curv_conv_cont(:,i) = sum_curv_fluo(focus_idx_conv_cont(i)+[-start:duration-start-1]);
    param_head_dist_conv_cont(:,i) = param_head_dist_fluo(focus_idx_conv_cont(i)+[-start:duration-start-1]);
    param_head_angle_conv_cont(:,i) = param_head_angle_fluo(focus_idx_conv_cont(i)+[-start:duration-start-1]);
end
disp('run done!');
%%
disp('this section running!');
%the single trials trace of focus region
roi_pos = [362 291 104;306 331 164;284 166 99;312 75 128;231 298 121;228 159 100];%
chnname = {'ot1','fb1','mb1','hb1','contra-ot','hb2'};
%ot1, fb1,mb1, hb1, contra-ot,hb2
[~,roi_idx] = arrayfun(@(i) min(sum((center -roi_pos(i,:)).^2,2)),1:size(roi_pos,1));
% multchnPlot(permute(mean(spike_focus(roi_idx,:,:),3),[2 3 1]),[-0.5,0.1*(duration-5)],[0],[],chnname,'plotstyle','plot');
% for i=1:size(spike_focus,3)
m_spike_focus_conv_diff = permute(cat(3,median(spike_focus_conv_start,3),median(spike_focus_conv_cont,3)),[2 3 1]);
std_spike_focus_conv_diff = permute(cat(3,std(spike_focus_conv_start,[],3)/sqrt(length(focus_idx_conv_start)),std(spike_focus_conv_cont,[],3)/sqrt(length(focus_idx_conv_cont))),[2 3 1]);
multchnPlot(m_spike_focus_conv_diff(:,:,roi_idx),[-0.1*start,0.1*(duration-start-1)],[0],std_spike_focus_conv_diff(:,:,roi_idx),chnname,'plotstyle','boundplot');
% end
ax = findobj(gcf,'type','axes');
for i=1:length(roi_idx)
    set(gcf,'currentaxes',ax(i));hold on;
    yyaxis right;
    plot([-0.1*start:0.1:0.1*(duration-start-1)],median(abs(sum_curv_conv_start),2),'Color',defaultColor(1),'LineStyle','--');
    plot([-0.1*start:0.1:0.1*(duration-start-1)],median(abs(sum_curv_conv_cont),2),'Color',defaultColor(2),'LineStyle','--');
end

disp('run done!');
%%
disp('running this section...');
wholebrainmapping(squeeze(m_spike_focus_conv_diff(:,1,:))',median(param_head_angle_conv_start,2),csessionID,cfishID,center,cat(1,zeros(duration,1)',median(abs(sum_curv_conv_start),2)'),'conv_start');
wholebrainmapping(squeeze(m_spike_focus_conv_diff(:,2,:))',median(param_head_angle_conv_cont,2),csessionID,cfishID,center,cat(1,zeros(duration,1)',median(abs(sum_curv_conv_cont),2)'),'conv_cont');
disp('this section done!');
%%
figure,
boxplot2(m_spike_focus_conv_diff(:,:,:));
xlabel('t');ylabel('activity amplitude');
%%
disp('this section running...');
%granger causality
X = spike_focus_conv_start(roi_idx,:,:);%roi x time x trial
[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,5,'LWR');
morder = moAIC;
[A,SIG] = tsdata_to_var(X,morder,'LWR');
assert(~isbad(A),'VAR estimation failed');
[G,info] = var_to_autocov(A,SIG,4);
var_info(info,true);
F = autocov_to_pwcgc(G);
assert(~isbad(F,false),'GC calculation failed');
pval = mvgc_pval(F,morder,size(X,2),size(X,3),1,1,size(X,1)-2, tstat);
sig = significance(pval,0.05,'FDR');
figure,
subplot(1,3,1),
plot_pw(F);
title('Pairwise-conditional GC');
subplot(1,3,2),
plot_pw(pval);
title('p-values');
subplot(1,3,3),
plot_pw(sig);
title('significant at p=0.05' );
disp('this section done!');