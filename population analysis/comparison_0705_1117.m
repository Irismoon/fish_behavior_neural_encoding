%compare 0705 and 1117 eye convergence difference
%1.extract the data around eye convergence
%%
load(fullfile(getpath('neural activity','201117','1'),'spike_OASIS'));
spike_1117 = sMatrix_total;
load(fullfile(getpath('neural activity','200705','1'),'spike_OASIS'));
spike_0705 = sMatrix_total;
%%
%2.only select the segment with behavior
load(fullfile(getpath('behavior','201117','1'),'tail_swing'),'good_bout_idx');
load(fullfile(getpath('behavior','201117','1'),'align_with_fluo'));
align_good_bout_idx_to_neuro;
conv_idx = fluo_good_bout_idx(fluo_good_bout_idx(:,4)==1,:);
spike_behavior_1117 = arrayfun(@(i,j) spike_1117(:,i-5:j),conv_idx(:,1),conv_idx(:,2),'un',0);
spike_behavior_1117 = cat(2,spike_behavior_1117{:});
%%
load(fullfile(getpath('behavior','200705','1'),'conv_analysis.mat'),'conv_or_not');
load(fullfile(getpath('behavior','200705','1'),'align_with_fluo'));
good_bout_idx = reshape(find(diff(conv_or_not)),2,[])';
mask = arrayfun(@(i) all(ismember(good_bout_idx(i,:),find(align_with_fluo_low))),1:size(good_bout_idx,1));
good_bout_idx = good_bout_idx(mask,:);
good_bout_idx = reshape(arrayfun(@(i) find(i==find(align_with_fluo_low)),good_bout_idx(:)),[],2);
conv_or_not = row2col(conv_or_not(align_with_fluo_low==1),1);
idx = 1:5:length(conv_or_not);
clear conv_idx;
[~,conv_idx(:,1)] = min(abs(good_bout_idx(:,1) - col2row(idx,1)),[],2);
[~,conv_idx(:,2)] = min(abs(good_bout_idx(:,2) - col2row(idx,1)),[],2);
spike_behavior_0705 = arrayfun(@(i,j) spike_0705(:,i-5:j),conv_idx(:,1),conv_idx(:,2),'un',0);
spike_behavior_0705 = cat(2,spike_behavior_0705{:});
%%
good_region_0705 = sum(spike_behavior_0705,2)~=0;
good_region_1117 = sum(spike_behavior_1117,2)~=0;
spike_behavior_0705 = spike_behavior_0705(good_region_0705,:);
spike_behavior_1117 = spike_behavior_1117(good_region_1117,:);
%%
%2.covariance matrix
corr_0705 = corrcoef(spike_behavior_0705','Rows','all');
corr_1117 = corrcoef(spike_behavior_1117','Rows','pairwise');
%%
%3.clustering
Z_0705 = linkage((abs(corr_0705-1)),'complete');
Z_1117 = linkage((abs(corr_1117-1)),'complete');
figure,
subplot(1,2,1),
dendrogram(Z_0705);
subplot(1,2,2),
dendrogram(Z_1117);
%%
kc = 4;
label_0705 = cluster(Z_0705,'MaxClust',kc);
label_1117 = cluster(Z_1117,'MaxClust',kc); 
%%
load(fullfile(getpath('neural activity','200705','1'),'A3_judge'),'center');
color = brewermap(kc,'Set1');
figure,hold on;
center_0705 = center(good_region_0705,:);
for i=1:kc
subplot(2,kc,i)
scatter3(center_0705(label_0705==i,1),center_0705(label_0705==i,2),center_0705(label_0705==i,3),10,color(i,:),'filled');
view(0,90);
end
load(fullfile(getpath('neural activity','201117','1'),'Coherence3'),'center');
% figure,hold on;
center_1117 = center(good_region_1117,:);
for i=1:kc
subplot(2,kc,i+kc),
scatter3(center_1117(label_1117==i,1),center_1117(label_1117==i,2),center_1117(label_1117==i,3),10,color(i,:),'filled');
view(0,90);
end
%%
count_0705 = histcounts(label_0705);
count_1117 = histcounts(label_1117);
figure,
subplot(1,2,1),
pie(count_0705/sum(count_0705));
subplot(1,2,2),
pie(count_1117/sum(count_1117));
%%
[W_0705,S_0705,V_0705] = svd(spike_behavior_0705','econ');
[W_1117,S_1117,V_1117] = svd(spike_behavior_1117','econ');
S_0705 = diag(S_0705)/sum(diag(S_0705));
S_1117 = diag(S_1117)/sum(diag(S_1117));
figure,
loglog(S_0705);hold on;
loglog(S_1117);
legend({'0705','1117'});
figure,
k = 1;
subplot(1,2,1),
scatter3(center_0705(:,1),center_0705(:,2),center_0705(:,3),10,V_0705(:,k));colormap('jet');
view(0,90);
title('0705');
subplot(1,2,2),
scatter3(center_1117(:,1),center_1117(:,2),center_1117(:,3),10,V_1117(:,k));colormap('jet');
view(0,90);
title('1117');
