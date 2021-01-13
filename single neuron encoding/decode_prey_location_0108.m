% from the selected good bouts, decode the prey location in the visual field in a time-varying way.
sessionID = '200705';fishID = '1';
load(fullfile(getpath('behavior',sessionID,fishID),'tail_swing'),'sum_curv');
load(fullfile(getpath('neural activity',sessionID,fishID),'DFoF'),'DFoF');
load(fullfile(getpath('neural activity',sessionID,fishID),'Coherence3'),'center');
neuralSignal = DFoF;
if ismember(sessionID,{'200108','200705'})
load(fullfile(getpath('behavior',sessionID,fishID),'align_with_fluo'));
sum_curv = sum_curv(align_with_fluo_low==1);
end
sum_curv = sum_curv(1:5:end);
numRegion = size(neuralSignal,1);
%%
%simple regression on sum_curv
model = arrayfun(@(i) fitlm(neuralSignal(i,:),sum_curv),1:numRegion,'un',0);
p = cellfun(@(x) x.Coefficients.pValue(2),model);
tunedRegion = find(p<0.001/numRegion);
%%
figure,
scatter3(center(:,1),center(:,2),center(:,3),10,'k');
hold on;
scatter3(center(tunedRegion,1),center(tunedRegion,2),center(tunedRegion,3),10,'r','filled');