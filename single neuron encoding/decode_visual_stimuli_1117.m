%1117 visual encoding analysis
sessionID  ='201117';fishID = '1';
load(fullfile(getpath('behavior',sessionID,fishID),'low_video_analysis_result'),'param_head_dist_all','param_head_angle_all','no_param');
load(fullfile(getpath('behavior',sessionID,fishID),'tail_swing'),'good_bout_idx','sum_curv');
load(fullfile(getpath('behavior',sessionID,fishID),'align_with_fluo'));
load(fullfile(getpath('neural activity',sessionID,fishID),'Coherence3'),'center');

no_param = no_param(align_with_fluo_low==1);
no_param = no_param(1:5:end);
load(fullfile(getpath('neural activity',sessionID,fishID),'spike_OASIS'));
spike_raw = sMatrix_total;
discard_region = find(sum(sMatrix_total,2)==0);
spike_raw(discard_region,:) = [];
center(discard_region,:) = [];
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
% figure,
% polarplot(deg2rad(param_head_angle_all),param_head_dist_all);
%since distance and angle is closely coupled, only one is considered here
param_head_dist_all = param_head_dist_all(align_with_fluo_low==1);
param_head_dist_fluo = param_head_dist_all(1:5:end-5);
param_head_angle_fluo = param_head_angle_all(align_with_fluo_low==1);
param_head_angle_fluo = param_head_angle_fluo(1:5:end-5);
figure,hold on;
histogram(param_head_dist_fluo(hunting_bout_fluo_idx_whole),'BinWidth',10,'Normalization','probability');hold on;
histogram(param_head_dist_fluo(spon_bout_fluo_idx_whole),'BinWidth',10,'Normalization','probability');
histogram(param_head_dist_fluo(rest_bout_fluo_idx),'BinWidth',10,'Normalization','probability');
legend({'hunting','spon move','rest'});
xlabel('prey distance');ylabel('probability');
p1 = ranksum(param_head_dist_fluo(hunting_bout_fluo_idx_whole),param_head_dist_fluo(spon_bout_fluo_idx_whole))
p2 = ranksum(param_head_dist_fluo(spon_bout_fluo_idx_whole),param_head_dist_fluo(rest_bout_fluo_idx(randperm(length(rest_bout_fluo_idx),length(spon_bout_fluo_idx_whole)))))
%p2 is small because of the large sample size
%%
sum_curv = sum_curv(align_with_fluo_low==1);
sum_curv_fluo = sum_curv(1:5:end-5);
%%
index = zeros(numTime,1);
index(hunting_bout_fluo_idx_whole) = 1;
index(spon_bout_fluo_idx_whole) = 2;
index(rest_bout_fluo_idx) = 3;
p_anovan = arrayfun(@(i) anovan(spike_raw(i,:),{index,param_head_dist_fluo,sum_curv_fluo},'model','interaction','continuous',[2 3],'display','off','varnames',{'beh_states','distance','sum curv'})',1:numRegion,'un',0);
p_anovan = cat(1,p_anovan{:});

[N_behstate,bin_behstate] = histcounts(p_anovan(:,1),'BinWidth',.001/numRegion,'Normalization','probability');%set(gca,'XScale','log','YScale','linear');
[N_dist,bin_dist] = histcounts(p_anovan(:,2),'BinWidth',.001/numRegion,'Normalization','probability');
[N_interact,bin_interact] = histcounts(p_anovan(:,3),'BinWidth',.001/numRegion,'Normalization','probability');
figure,hold on;
plot(bin_behstate(2:end),N_behstate);
plot(bin_dist(2:end),N_dist);
plot(bin_interact(2:end),N_interact);
set(gca,'XScale','log');
legend({'behavior states','distance','interaction'});
xlabel('p value');ylabel('probability');
title('how many regions is affected by different variables');
%%
%generalized linear mixed effects model
pValue = zeros(numRegion,3);
for i=1:numRegion
tbl = table(spike_raw(i,:)',sum_curv_fluo,param_head_dist_fluo,param_head_angle_fluo,categorical(index),'VariableNames',{'spike','sum_curv','head_distance','head_angle','behavior_states'});
formula = 'spike ~ sum_curv + head_distance + (1|behavior_states)';
lme = fitlme(tbl,formula);
formula = 'spike ~ sum_curv + head_distance + (1|behavior_states) + (sum_curv-1|behavior_states)';
lme1 = fitlme(tbl,formula);
result = compare(lme,lme1);%,'CheckNesting',true);
pValue(i,1) = result.pValue(2);
formula = 'spike ~ sum_curv + head_distance + (1|behavior_states) + (sum_curv-1|behavior_states) + (head_distance-1|behavior_states)';
lme1 = fitlme(tbl,formula);
result = compare(lme,lme1);%,'CheckNesting',true);
pValue(i,2) = result.pValue(2);
formula = 'spike ~ sum_curv + head_distance + (1|behavior_states) + (sum_curv-1|behavior_states) + (head_distance-1|behavior_states) + (head_angle-1|behavior_states)';
lme1 = fitlme(tbl,formula);
result = compare(lme,lme1);%,'CheckNesting',true);
pValue(i,3) = result.pValue(2);
% glme = fitglme(tbl,formula,'Distribution','Gamma','link',0.5);
end
%%
%find regions for each factor
figure,hold on;
scatter3(center(:,1),center(:,2),center(:,3),10,'k');
scatter3(center(pValue(:,1)<0.001/numRegion,1),center(pValue(:,1)<0.001/numRegion,2),center(pValue(:,1)<0.001/numRegion,3),10,'r','filled');
scatter3(center(pValue(:,2)<0.001/numRegion,1),center(pValue(:,2)<0.001/numRegion,2),center(pValue(:,2)<0.001/numRegion,3),10,'b','filled');
scatter3(center(pValue(:,3)<0.001/numRegion,1),center(pValue(:,3)<0.001/numRegion,2),center(pValue(:,3)<0.001/numRegion,3),10,'g','filled');

%%
%To-do:
%1.check the neurons sensitive to behavior states, distance and interaction
%are overlapping or not
figure,hold on;
color = brewermap(size(p_anovan,2),'Set1');
[nrow,ncol] = arrange_subplots(size(p_anovan,2),[1 1]);
for i=1:size(p_anovan,2)
idx = p_anovan(:,i)<(0.001/numRegion);
subplot(nrow,ncol,i),hold on;
scatter3(center(:,1),center(:,2),center(:,3),10,'k');
scatter3(center(idx,1),center(idx,2),center(idx,3),10,color(i,:),'filled');
end
%%
figure,
imagesc(corrcoef([spike_raw(p_anovan(:,1)<(0.001/numRegion),:);spike_raw(p_anovan(:,2)<(0.001/numRegion),:);spike_raw(p_anovan(:,3)<(0.001/numRegion),:)]'))
%%
%cluster these neurons and find out the clusters with different activity
%modality
k = 3;
color = brewermap(k,'Set1');
for i=1:size(p_anovan,2)
    mask = find(p_anovan(:,i)<(0.001/numRegion));
    spike_tmp = spike_raw(mask,:)';
    corr = corrcoef(spike_tmp);%this is not good because 
    Z = linkage(squareform(abs(corr-1)),'complete');
    figure,
    subplot(1,3,1),dendrogram(Z);
    T = cluster(Z,'maxclust',k);
    subplot(1,3,2),
    imagesc(corr([find(T==1);find(T==2)],[col2row(find(T==1)) col2row(find(T==2))]));
    subplot(1,3,3),
    scatter3(center(:,1),center(:,2),center(:,3),10,'k');
    hold on;
    arrayfun(@(i) scatter3(center(mask(T==i),1),center(mask(T==i),2),center(mask(T==i),3),10,color(i,:),'filled'),1:k,'un',0);
end
%%
%check the validness of the anovan result
distance_idx = find(p_anovan(:,2)<(0.001/numRegion));
region_idx= distance_idx(randperm(length(distance_idx),1));

figure,subplot(1,2,1),hold on;
current_idx = spon_bout_fluo_idx_whole;
current_variable = param_head_angle_fluo;
scatter(current_variable(current_idx),spike_raw(region_idx,current_idx),10,'r','filled');
subplot(1,2,2),hold on;
current_variable = param_head_dist_fluo;
scatter(current_variable(current_idx),spike_raw(region_idx,current_idx),10,'r','filled');

current_idx = rest_bout_fluo_idx;
subplot(1,2,1),
current_variable = param_head_angle_fluo;
scatter(current_variable(current_idx),spike_raw(region_idx,current_idx),10,'b','filled');
subplot(1,2,2),
current_variable = param_head_dist_fluo;
scatter(current_variable(current_idx),spike_raw(region_idx,current_idx),10,'b','filled');

current_idx = hunting_bout_fluo_idx_whole;
subplot(1,2,1),
current_variable = param_head_angle_fluo;
scatter(current_variable(current_idx),spike_raw(region_idx,current_idx),10,'k','filled');
xlabel('angle');
subplot(1,2,2),
current_variable = param_head_dist_fluo;
scatter(current_variable(current_idx),spike_raw(region_idx,current_idx),10,'k','filled');
xlabel('distance');
%%
distance_idx = find(p_anovan(:,2)<(0.001/numRegion));
current_idx = spon_bout_fluo_idx_whole;
region_idx= distance_idx(randperm(length(distance_idx),1));
valid_idx= [spon_bout_fluo_idx_whole;row2col(rest_bout_fluo_idx,1);hunting_bout_fluo_idx_whole];
marker = cell(length(valid_idx),1);
marker(1:length(spon_bout_fluo_idx_whole)) = {'o'};
marker(1+length(spon_bout_fluo_idx_whole):length(spon_bout_fluo_idx_whole)+length(rest_bout_fluo_idx)) = {'s'};
marker(length(spon_bout_fluo_idx_whole)+length(rest_bout_fluo_idx)+1:end) = {'*'};

param_head_angle_tight = param_head_angle_fluo(valid_idx);
param_head_dist_tight = param_head_dist_fluo(valid_idx);
spike_raw_tight = spike_raw(:,valid_idx);
figure,colormap('jet');
for i=1:length(valid_idx)
    caxis([min(spike_raw_tight(region_idx,:)) max(spike_raw_tight(region_idx,:))]);
polarscatter(deg2rad(param_head_angle_tight(i)),param_head_dist_tight(i),10,spike_raw_tight(region_idx,i),marker{i});
hold on;
pause(0.05);
end
%%
region_idx= distance_idx(randperm(length(distance_idx),1));
figure,
subplot(1,3,1),
current_idx = spon_bout_fluo_idx_whole;
polarscatter(deg2rad(param_head_angle_fluo(current_idx)),param_head_dist_fluo(current_idx),10,spike_raw(region_idx,current_idx),'o');
caxis([min(spike_raw_tight(region_idx,:)) max(spike_raw_tight(region_idx,:))]);
subplot(1,3,2),
current_idx = rest_bout_fluo_idx;
polarscatter(deg2rad(param_head_angle_fluo(current_idx)),param_head_dist_fluo(current_idx),10,spike_raw(region_idx,current_idx),'o');
caxis([min(spike_raw_tight(region_idx,:)) max(spike_raw_tight(region_idx,:))]);
subplot(1,3,3),
current_idx = hunting_bout_fluo_idx_whole;
polarscatter(deg2rad(param_head_angle_fluo(current_idx)),param_head_dist_fluo(current_idx),10,spike_raw(region_idx,current_idx),'o');
caxis([min(spike_raw_tight(region_idx,:)) max(spike_raw_tight(region_idx,:))]);