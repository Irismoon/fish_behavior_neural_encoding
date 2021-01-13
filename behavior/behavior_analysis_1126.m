%%
csessionID = '201117';cfishID = '1';
%%
%1.how many convergence
load(fullfile(getpath('behavior',csessionID,cfishID),'high_analysis'),'conv_or_not');
conv_idx = reshape(find(diff(conv_or_not)),2,[])'+1;
figure,
histogram(conv_idx(:,2)-conv_idx(:,1),'BinWidth',1);
title(['total convergence times: ' num2str(size(conv_idx,1))]);
%%
%2.tail direction relates to prey location
load(fullfile(getpath('behavior',csessionID,cfishID),'tail_swing'),'conv_bout_idx','angledata');
load(fullfile(getpath('behavior',csessionID,cfishID),'low_video_analysis_result'),'param_head_angle_all');
angledata = rad2deg(angledata)-90;
%%
nBout = size(conv_bout_idx,1);
boutAmp = zeros(nBout,21);
param_angle = zeros(nBout,1);
max_amp = zeros(nBout,1);
for i=1:nBout
    idx = conv_bout_idx(i,1):conv_bout_idx(i,2);
    boutAmp(i,:) = mean(angledata(idx,:));
    param_angle(i) = mean(param_head_angle_all(conv_bout_idx(i,1)+[-5:2]));
    max_amp(i) = maxabs(angledata(idx,:));
end
figure,
scatter(param_angle,max_amp);
hold on;grid on;
scatter(param_angle,boutAmp(:,end));
mdl = fitlm(param_angle,max_amp);
plot(param_angle,mdl.Fitted);
fplot(@(x) .5*x);
legend({'maximum angle','tail tip angle','fitted','should be'});
xlabel('param angle');ylabel('maximum tail angle')
%%
figure,
histogram(conv_bout_idx(:,2)-conv_bout_idx(:,1)+1);
xlabel('bout duration');ylabel('number');
%%
%how the eye angle changes during convergence
conv_idx = reshape(find(diff(conv_or_not)),2,[]);
conv_angle = arrayfun(@(i) converge_angle(conv_idx(1,i):conv_idx(2,i)),1:size(conv_idx,2),'un',0);
figure,hold on,
cellfun(@(x) plot(x),conv_angle);