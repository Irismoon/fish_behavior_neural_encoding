% from the selected good bouts, decode the prey location in the visual field in a time-varying way.
load(fullfile(getpath('behavior','200705',''),'tail_swing'),'good_bout_idx');
load(fullfile(getpath('behavior','200705',''),'align_with_fluo'),'align_with_fluo_low');
load(fullfile(getpath('behavior','200705',''),'low_video_analysis_result'),'param_head_angle_all','param_head_dist_all','param_speed_all','no_param');
load(fullfile(getpath('neural activity','200705',''),'Spike_X_EstTrace'));
param_head_angle = mean(reshape(param_head_angle_all(reshape((good_bout_idx(:,1)-[1:20])',[],1)),20,[]),1);
%0705 high and low was not aligned to fluo because good_bout_idx was based
%on high and low video (they have the same length) 
[lia,idx] = ismember(good_bout_idx(:,1),find(align_with_fluo_low));
mask = 1:5:nnz(align_with_fluo_low);
[~,fluo_idx] = min(abs(row2col(idx(lia)) - col2row(mask)),[],2);
param_head_angle = param_head_angle(lia);
leftIdx = param_head_angle>0;
rightIdx = param_head_angle<0;
numRegion = size(Spike_X_EstTrace,1);
numTrial = max(nnz(leftIdx),nnz(rightIdx));

numT = 9;
spk = zeros(length(fluo_idx),size(Spike_X_EstTrace,1),numT);
for i = 1:length(fluo_idx)
    spk(i,:,:) = Spike_X_EstTrace(:,fluo_idx(i)+(-floor(numT/2):floor(numT/2)));%region x 5
end
m_spk = squeeze(mean(spk(param_head_angle>0,:,:),1));
%map them on MIP
generate_MIP(m_spk+0.01*max(m_spk(:)),'200705','');

%fill mean value as missing data
missing_data = mean(spk(rightIdx,:,:));%1 x region x time
missing_data = repmat(missing_data,nnz(leftIdx)-nnz(rightIdx),1,1);
spk_balanced = cat(1,spk(leftIdx,:,:),spk(rightIdx,:,:),missing_data);%
spk_balanced = permute(reshape(spk_balanced,nnz(leftIdx),2,numRegion,numT),[3,2,4,1]);%region x direction x time x trial
% [W,V,whichMarg] = drtoolbox.techniques.dpca(mean(spk_balanced,4),5,'combinedParams',{{1},{2},{1,2}});
spk_balanced_ave = mean(spk_balanced,4);%neuron x direction x time
X = spk_balanced_ave(:,:);
X = X - mean(X,2);
[W,explVar,V] = svd(X,'econ');
W = W(:,1:18);
drtoolbox.techniques.dpca_plot(permute(spk_balanced,[1 4 2 3]), W, W, @drtoolbox.techniques.dpca_plot_default);
figure,
plot(explVar/sum(explVar));axis tight;xlabel('comp');ylabel('explained variance');
stdx = permute(std(pagemtimes(W',spk_balanced),[],4)./sqrt([nnz(leftIdx),nnz(rightIdx)]),[3 2 1]);%
multchnPlot(permute(pagemtimes(W',spk_balanced_ave),[3 2 1]),0.1*(-4:4),[0],'boundwidth',stdx,...
    'plotstyle','boundplot','chnname',arrayfun(@(i) ['comp ' num2str(i)],1:18,'un',0));
sgtitle('signal projection on each component');
xlabel('t/s');ylabel('firing rate');


%reconstruct the brain map time series from the PCs
%left
leftpc = [11 12 13 14 15 16 17];
rightpc = [2 4 5 6 7];
S = pagemtimes(W',spk_balanced_ave);%pc x direction x time
recons_left = pagemtimes(W(:,leftpc),S(leftpc,:,:));%region x direction x time
recons_right = pagemtimes(W(:,rightpc),S(rightpc,:,:));
%plot the brain map of PC
figure,
generate_MIP(squeeze(recons_left(:,1,:)),'200705','','_left_recons_left');
figure,
generate_MIP(squeeze(recons_left(:,2,:)),'200705','','_left_recons_right');
figure,
generate_MIP(squeeze(recons_right(:,2,:)),'200705','','_right_recons_left');
figure,
generate_MIP(squeeze(recons_right(:,2,:)),'200705','','_right_recons_right');

[B,fitInfo] = lassoglm(spk,param_head_angle'>0,'binomial','CV',5,'NumLambda',1000);
%%
[param_head_angle_all,param_speed_all,param_head_dist_all,no_param] = samfnmultvar(@(x) x(align_with_fluo_low==1),...
    param_head_angle_all,param_speed_all,param_head_dist_all,no_param);
[param_head_angle_fluo,param_speed_fluo,param_head_dist_fluo,no_param_fluo] = samfnmultvar(@(x) x(1:5:length(x)),param_head_angle_all,param_speed_all,param_head_dist_all,no_param);
idx = removeStrangeVisuoData(no_param_fluo,param_head_dist_fluo,param_speed_fluo);
param_head_angle_visuo = param_head_angle_fluo(idx);
spikes_visuo = Spike_X_EstTrace(:,idx);
tuned2angle(param_head_angle_visuo,spikes_visuo,center);