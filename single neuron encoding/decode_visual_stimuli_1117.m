%1117 visual encoding analysis
sessionID  ='201117';fishID = '1';
load(fullfile(getpath('behavior',sessionID,fishID),'low_video_analysis_result'),'param_head_dist_all','param_head_angle_all','no_param');
load(fullfile(getpath('behavior',sessionID,fishID),'tail_swing'),'good_bout_idx');
load(fullfile(getpath('behavior',sessionID,fishID),'align_with_fluo'));
no_param = no_param(align_with_fluo_low==1);
no_param = no_param(1:5:end);
load(fullfile(getpath('neural activity',sessionID,fishID),'spike_OASIS'));
spike_raw = sMatrix_total;
delay = 3;
[numRegion,numTime] = size(spike_raw);
k = 5;

%%
hunting_bout_fluo_idx = align_behavior2neuro(good_bout_idx,align_with_fluo_low);
hunting_bout_fluo_idx_whole = arrayfun(@(i) [hunting_bout_fluo_idx(i,1)-k:hunting_bout_fluo_idx(i,2)]',1:size(hunting_bout_fluo_idx,1),'un',0);
hunting_bout_fluo_idx_whole = cat(1,hunting_bout_fluo_idx_whole{:});
%%
load(fullfile(getpath('behavior',sessionID,fishID),'tail_swing'),'bout_idx');
bout_fluo_idx = align_behavior2neuro(bout_idx,align_with_fluo_low);
%%
%spontaneous movement
mask = size(hunting_bout_fluo_idx,1);
for i=1:size(hunting_bout_fluo_idx,1)
    mask(i) = find(arrayfun(@(j) isequal(bout_fluo_idx(j,:),hunting_bout_fluo_idx(i,:)),1:size(bout_fluo_idx,1)));
end
spon_bout_fluo_idx = bout_fluo_idx;
spon_bout_fluo_idx(mask,:) = [];
spon_bout_fluo_idx_whole = arrayfun(@(i) [spon_bout_fluo_idx(i,1)-k:spon_bout_fluo_idx(i,2)]',1:size(spon_bout_fluo_idx,1),'un',0);
spon_bout_fluo_idx_whole = cat(1,spon_bout_fluo_idx_whole{:});
%%
%saccade
% load(fullfile(getpath('behavior',sessionID,fishID),'saccade'),'left_start','left_end','right_start','right_end');
% saccade = [row2col(left_start,1) row2col(left_end,1);row2col(right_start,1) row2col(right_end,1)];
% saccade_fluo_idx = align_behavior2neuro(saccade,align_with_fluo_low);
%%
mask = 1:numTime;
rest_bout_fluo_idx = setdiff(mask,[hunting_bout_fluo_idx_whole;spon_bout_fluo_idx_whole]);
%%
hunting_bout_fluo_idx_whole(no_param(hunting_bout_fluo_idx_whole)==0) = [];
spon_bout_fluo_idx_whole(no_param(spon_bout_fluo_idx_whole)==0) = [];
rest_bout_fluo_idx(no_param(rest_bout_fluo_idx)==0) = [];
%%
%the trajectory of paramecium
figure,
polarplot(deg2rad(param_head_angle_all),param_head_dist_all);
%since distance and angle is closely coupled, only one is considered here
param_head_dist_all = param_head_dist_all(align_with_fluo_low==1);
param_head_dist_fluo = param_head_dist_all(1:5:end-5);
figure,hold on;
histogram(param_head_dist_fluo(hunting_bout_fluo_idx_whole),'BinWidth',10,'Normalization','probability');hold on;
histogram(param_head_dist_fluo(spon_bout_fluo_idx_whole),'BinWidth',10,'Normalization','probability');
histogram(param_head_dist_fluo(rest_bout_fluo_idx),'BinWidth',10,'Normalization','probability');
legend({'hunting','spon move','rest'});
xlabel('prey distance');ylabel('probability');
%%
index = zeros(numTime,1);
index(hunting_bout_fluo_idx_whole) = 1;
index(spon_bout_fluo_idx_whole) = 2;
index(rest_bout_fluo_idx) = 3;
p = arrayfun(@(i) anovan(spike_raw(i,:),{index,param_head_dist_fluo},'model','interaction','continuous',2,'display','off','varnames',{'beh_states','distance'})',1:numRegion,'un',0);
p = cat(1,p{:});

[N_behstate,bin_behstate] = histcounts(p(:,1),'BinWidth',.001,'Normalization','probability');%set(gca,'XScale','log','YScale','linear');
[N_dist,bin_dist] = histcounts(p(:,2),'BinWidth',.001,'Normalization','probability');
[N_interact,bin_interact] = histcounts(p(:,3),'BinWidth',.001,'Normalization','probability');
figure,hold on;
plot(bin_behstate(2:end),N_behstate);
plot(bin_dist(2:end),N_dist);
plot(bin_interact(2:end),N_interact);
legend({'behavior states','distance','interaction'});
%%
%To-do:
%1.check the neurons sensitive to behavior states, distance and interaction
%are overlapping or not
