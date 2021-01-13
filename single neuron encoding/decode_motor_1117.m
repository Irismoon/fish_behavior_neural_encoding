sessionID = '201117';fishID = '1';
load(fullfile(getpath('behavior',sessionID,fishID),'tail_swing'),'sum_curv','good_bout_idx');
load(fullfile(getpath('behavior',sessionID,fishID),'low_video_analysis_result'),'param_head_angle_all');
align_sum_curv = sum_curv(align_with_fluo_low==1);
align_param_head_angle_all = param_head_angle_all(align_with_fluo_low==1);
idx = ismember(good_bout_idx(:,1),find(align_with_fluo_low==1));
good_bout_idx = good_bout_idx(idx,:);
load(fullfile(getpath('neural activity',sessionID,fishID),'DFoF'),'DFoF');
load(fullfile(getpath('neural activity',sessionID,fishID),'Coherence3'),'A3');
load(fullfile(getpath('neural activity','image segmentation'),'ZBB anatomy','pixel2area'));
%%
pixelIdx = tbl.(2)(ismember(tbl.(1),{'inferior olive'}));
regionIdx = arrayfun(@(i) row2col(find(A3(i,:)),1),pixelIdx,'un',0);
regionIdx = cat(1,regionIdx{:});
regionIdx = unique(regionIdx);
%%
figure,
plot(DFoF(regionIdx,:)');
yyaxis right;
plot(align_sum_curv(1:5:end));
%%
idx = find(align_with_fluo_low);
idx = idx(1:5:end);
fluo_good_bout_idx = good_bout_idx;
[~,fluo_good_bout_idx(:,1)] = min(abs(good_bout_idx(:,1) - col2row(idx,1)),[],2);
[~,fluo_good_bout_idx(:,2)] = min(abs(good_bout_idx(:,2) - col2row(idx,1)),[],2);
%%
%align to bout start
align_DFoF = arrayfun(@(i) DFoF(regionIdx,i+[-2:5])'-mean(DFoF(regionIdx,i+[-6:-3]),2)',fluo_good_bout_idx(:,1),'un',0);
align_DFoF = cat(3,align_DFoF{:});%10 x region x trial
figure,plot(mean(align_DFoF,3));
%%
m_sum_curv = arrayfun(@(x) mean(sum_curv(x+[0:3])),good_bout_idx(:,1));
m_param_angle = arrayfun(@(x) mean(param_head_angle_all(x+[-3:0])),good_bout_idx(:,1));
% varx = m_param_angle;
varx = m_sum_curv;
%%
m_DFoF = arrayfun(@(x) mean(DFoF(regionIdx,x+[-1:2])'),fluo_good_bout_idx(:,1),'un',0);
m_DFoF = cat(1,m_DFoF{:});%bout x region
figure,
for i=1:8
    subplot(2,4,i),
    scatter(varx,m_DFoF(:,i),10,'filled');
    hold on;
    xline(0,'k--');yline(0,'k--');
end
xlabel('sum curv');ylabel('neural activity');
sgtitle('inferior olive tuning to tail');