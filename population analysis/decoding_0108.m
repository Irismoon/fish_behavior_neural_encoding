missing_data = mean(spk(lessIdx,:,:));%1 x region x time
missing_data = repmat(missing_data,nnz(moreIdx)-nnz(lessIdx),1,1);
spk_balanced = cat(1,spk(moreIdx,:,:),spk(lessIdx,:,:),missing_data);%
spk_balanced = permute(reshape(spk_balanced,nnz(moreIdx),2,numRegion,numT),[3,2,4,1]);%region x direction x time x trial

X = cat(3,squeeze(spk_balanced(:,1,:,:)),squeeze(spk_balanced(:,2,:,1:nnz(lessIdx))));%region x time x trial
proj_X = pagemtimes(W(:,1:3)',X);%3 x time x trial
label = [zeros(nnz(moreIdx),1);ones(nnz(lessIdx),1)];
X = reshape(proj_X,3*numT,length(label))';
% X = repmat(X,2,1);
% label = repmat(label,2,1);
numFold = 3;
c = cvpartition(label,'KFold',numFold);
accu = zeros(numFold,1);
Coeff = zeros(1+size(X,2),numFold);
for i=1:numFold
    idxTest = c.test(i);
    idxTrain = ~idxTest;
    [B,fitInfo] = lassoglm(X(idxTrain,:),label(idxTrain),'binomial','CV',5);
    idxMinDev = fitInfo.IndexMinDeviance;
    coeff = [fitInfo.Intercept(idxMinDev);B(:,idxMinDev)];
    Coeff(:,i) = coeff;
    labelhat = glmval(coeff,X(idxTest,:),'logit',fitInfo);
    accu(i) = mean((labelhat>=0.5)==label(idxTest));
end
%%
%%%%%%compare the decoding performance at different behavior states%%%%%%%
disp('running this section....');
%1.construct data
spk_hunting = Spike_X_EstTrace(:,hunting_idx_visuo);
spk_spon = Spike_X_EstTrace(:,nonhunting_idx_visuo);
param_hunting = arrayfun(@(i) mean(param_head_angle_fluo(hunting_idx_visuo(i)+[-5:0])),1:length(hunting_idx_visuo));
param_nonhunting = arrayfun(@(i) mean(param_head_angle_fluo(max([ones(1,6);nonhunting_idx_visuo(i)+[-5:0]]))),1:length(nonhunting_idx_visuo));
disp('this section done!');
%%
disp('running this section....');
%2.reduce dimension
[U,S,V] = svd(spk_hunting','econ');
S = diag(S);
include_idx_hunting = find(cumsum(S)/sum(S)>0.95,1);
proj_hunting = spk_hunting'*V(:,1:include_idx_hunting);
[Us,Ss,Vs] = svd(spk_spon','econ');
Ss = diag(Ss);
include_idx_spon = find(cumsum(Ss)/sum(Ss)>0.95,1);
proj_spon = spk_spon'*Vs(:,1:include_idx_spon);
label_hunting = discretize(param_hunting,linspace(min(param_hunting),max(param_hunting),6));
label_nonhunting = discretize(param_nonhunting,linspace(min(param_nonhunting),max(param_nonhunting),6));
disp('this section done!');
%%
%stepwise regression based on reduced dimension(1VSothers, or 1VS1VS1VS1)
%regression based on neurons with preferred angle
%other decoding methods
% dataset = [spk_hunting' row2col(label_hunting,1)];
% dataset = [spk_spon' row2col(label_nonhunting,1)];
% dataset = [proj_hunting row2col(label_hunting,1)];
% dataset = [proj_spon row2col(label_nonhunting,1)];
% dataset = [proj_hunting row2col(param_hunting,1)];
dataset = [proj_spon row2col(param_nonhunting,1)];
cvmse = crossval('mse',dataset(:,1:end-1),dataset(:,end),'Predfun',@regf,'KFold',10);
%%
disp('this section running...');
%visual encoding at single neuron level
%tuning curve for each neuron at different behavior states
%neuron x window
% m_spk_hunting = zeros(numRegion,5);
% m_spk_spon = zeros(numRegion,5);
% for i=1:5
%     m_spk_hunting(:,i) = mean(spk_hunting(:,label_hunting==i),2);
%     m_spk_spon(:,i) = mean(spk_spon(:,label_nonhunting==i));
% end
%curve fitting
% gfit = fittype('1/sqrt(2*pi)/sigma*exp(-((x-mu)/sigma)^2)+c','dependent',{'y'},'independent',{'x'},'coefficients',{'mu','sigma','c'});
% gfit = fittype('a*cos(x)+b','dependent',{'y'},'independent',{'x'},'coefficients',{'a','b'});
gfit = fittype('a*x+b','dependent',{'y'},'independent',{'x'},'coefficients',{'a','b'});
[cfit_hunting,gof_hunting] = arrayfun(@(i) fit(spk_hunting(i,:)',row2col(param_hunting,1),gfit),1:numRegion,'un',0);
[cfit_spon,gof_spon] = arrayfun(@(i) fit(spk_spon(i,:)',row2col(param_nonhunting,1),gfit),1:numRegion,'un',0);
rsq_hunting = cellfun(@(x) x.rsquare,gof_hunting);
rsq_spon = cellfun(@(x) x.rsquare,gof_spon);
figure,histogram(rsq_spon);
hold on;histogram(rsq_hunting);
legend({'spontaneous','hunting'});
xlabel('rsquare');ylabel('number of regions');
disp('this section done!');
%%
%12.15
disp('running this section......');
%search for all neurons activated differently during hunting and nonhunting
%trials
numTrial = size(spike_focus_conv_start,3);
%1.construct data, use resting state and select some as non-response trials
%should be far from hunting start or movement
rest_idx_trial = rest_idx_visuo(arrayfun(@(i) all(~ismember(rest_idx_visuo(i)+[1:duration+10],[hunting_idx_wholecourse;bout_idx_aligned_wholecourse])),1:length(rest_idx_visuo)));
rest_idx_trial = rest_idx_trial(rest_idx_trial+duration<size(Spike_X_EstTrace,2));
selc_region = zeros(numRegion,50);
for iPerm=1:50
rest_idx_trial_tmp = rest_idx_trial(randperm(length(rest_idx_trial),numTrial));
spike_rest_trial = arrayfun(@(i) Spike_X_EstTrace(:,rest_idx_trial_tmp(i)+[0:duration-1]) - Spike_X_EstTrace(:,rest_idx_trial_tmp(i)-1),1:length(rest_idx_trial_tmp),'un',0);
spike_rest_trial = cat(3,spike_rest_trial{:});
%2.only the part before convergence start/J-turn, 
%3.ttest or permutation test to find out the query neurons
[h,p] = arrayfun(@(i) ranksum(squeeze(mean(spike_rest_trial(i,:,:),2)),squeeze(mean(spike_focus_conv_start(i,1:3,:),2))),1:numRegion);
selc_region(:,iPerm) = h<0.01;
end
selc_region = all(selc_region,2);
generate_MIP(selc_region+0.5,csessionID,cfishID,'_hunting_initiation');
disp('this section done!');
%%
%1215
%to check if the habenula activity is related to eye movement
%1.extract eye movement start time
tmp = diff(right_eye_angle_fluo);
label = [false false false arrayfun(@(i) abs(tmp(i))>4*max(abs(tmp([i-3:i-1 i+1:i+3]))),4:length(tmp)-3) false false false];
% figure,
% plot(label);yyaxis right;
% plot(right_eye_angle_fluo);
% roi = [];
for regionIdx = roi(1:end)%col2row(find(selc_region),1)
conv_start = find(diff(conv_or_not_fluo)>0);
conv_end = find(diff(conv_or_not_fluo)<0);
right_eye_start = [16 72 201 263 336 479 614 671 825 839 915 1005 1112 1374 1521 1551 1729 1954 2090 2115 2157 2220];
right_eye_end = unique([97 149 234 327 395 487 702 900 947 1020 1157 1405 1541 1617 2165 2183 2209 col2row(conv_start,1)]);
left_eye_start = unique([17 72 100 201 260 328 394 479 612 671 791 839 902 1007 1111 1205 1367 1521 1544 1617 1729 1957 2088 2105 2148 2209 col2row(conv_start,1)]);
left_eye_end = [150 233 401 825 915 950 1039 1171 1212 1409 1552 1627 1773 1970 ];
%

%eye movement without convergence also induces neural activation?
left_saccade_start = left_eye_start(arrayfun(@(i) ~any(ismember(left_eye_start(i)+[-5:5],conv_start)),1:length(left_eye_start)));
left_saccade_end = left_eye_end(arrayfun(@(i) ~any(ismember(left_eye_end(i)+[-5:5],conv_start)),1:length(left_eye_end)));
right_saccade_start = right_eye_start(arrayfun(@(i) ~any(ismember(right_eye_start(i)+[-5:5],conv_start)),1:length(right_eye_start)));
right_saccade_end = right_eye_end(arrayfun(@(i) ~any(ismember(right_eye_end(i)+[-5:5],conv_start)),1:length(right_eye_end)));
%2.threshold the neural activity
spike_roi = Spike_X_EstTrace(regionIdx,:);
spike_roi_tmp = spike_roi;
% spike_roi_tmp = spike_roi(randperm(length(spike_roi)));
spike_or_not_roi = spike_roi_tmp>quantile(spike_roi,0.75);
%3.
p_right_eye_start = mean(arrayfun(@(i) any(spike_or_not_roi(right_eye_start(i)+[-3:0])),1:length(right_eye_start)))
p_right_eye_end = mean(arrayfun(@(i) any(spike_or_not_roi(right_eye_end(i)+[-3:0])),1:length(right_eye_end)))
p_left_eye_start = mean(arrayfun(@(i) any(spike_or_not_roi(left_eye_start(i)+[-3:0])),1:length(left_eye_start)))
p_left_eye_end = mean(arrayfun(@(i) any(spike_or_not_roi(left_eye_end(i)+[-3:0])),1:length(left_eye_end)))
p_conv_start = mean(arrayfun(@(i) any(spike_or_not_roi(conv_start(i)+[-3:0])),1:length(conv_start)))
p_conv_end = mean(arrayfun(@(i) any(spike_or_not_roi(conv_end(i)+[-3:0])),1:length(conv_end)))
p_left_saccade_start = mean(arrayfun(@(i) any(spike_or_not_roi(left_saccade_start(i)+[-3:0])),1:length(left_saccade_start)))
p_left_saccade_end = mean(arrayfun(@(i) any(spike_or_not_roi(left_saccade_end(i)+[-3:0])),1:length(left_saccade_end)))
p_right_saccade_start = mean(arrayfun(@(i) any(spike_or_not_roi(right_saccade_start(i)+[-3:0])),1:length(right_saccade_start)))
p_right_saccade_end = mean(arrayfun(@(i) any(spike_or_not_roi(right_saccade_end(i)+[-3:0])),1:length(right_saccade_end)))
p_left_saccade_idx = mean(arrayfun(@(i) any(spike_or_not_roi(left_saccade_idx(i)+[-3:0])),1:length(left_saccade_idx)))
p_right_saccade_idx = mean(arrayfun(@(i) any(spike_or_not_roi(right_saccade_idx(i)+[-3:0])),1:length(right_saccade_idx)))

figure,
bar([p_left_eye_start p_conv_start p_left_saccade_start]);
set(gca,'XTick',1:3,'XTickLabel',{'left increase','conv start','pure left increase'});
ylabel('spike probability');title(['region ' num2str(regionIdx) ' activation prob under eye movement states']);
legend({'original','shuffle'});
amp_left_eye_start = arrayfun(@(i) mean(spike_roi_tmp(left_eye_start(i)+[-3:0])),1:length(left_eye_start));
amp_conv_start = arrayfun(@(i) mean(spike_roi_tmp(conv_start(i)+[-3:0])),1:length(conv_start));
amp_left_saccade_start = arrayfun(@(i) mean(spike_roi_tmp(left_saccade_start(i)+[-3:0])),1:length(left_saccade_start));
figure,boxplot2({amp_left_eye_start,amp_conv_start,amp_left_saccade_start}','notch','on');
ylabel('mean firing rate');
figure,
plot(Spike_X_EstTrace(regionIdx,:));hold on;
yyaxis right;stem(conv_start,ones(nnz(conv_start),1));
% if p_conv_start>.85
%     roi = [roi regionIdx];
% end
end
%%
%how these rois activated during nonhunting state
spike_roi = Spike_X_EstTrace(roi,:);
spike_or_not_roi = spike_roi>quantile(spike_roi,0.75,2);
% spike_or_not_roi = [all(spike_or_not_roi([1],:),1);all(spike_or_not_roi([2:4],:),1);all(spike_or_not_roi([5:9],:),1);all(spike_or_not_roi(end,:),1)];
figure,
around_conv_idx_start = (focus_idx_conv_start+[-4:0])';
bar([mean(spike_or_not_roi(:,around_conv_idx_start(:)),2) mean(spike_or_not_roi(:,hunting_idx_visuo),2) mean(spike_or_not_roi(:,rest_idx_visuo),2) mean(spike_or_not_roi(:,simul_idx_visuo),2)])
%%
figure,
scatter3(center(:,1),center(:,2),center(:,3));hold on;
scatter3(center(roi(1),1),center(roi(1),2),center(roi(1),3),[],'r','filled');
scatter3(center(roi([2 3 4 end]),1),center(roi([2:4 end]),2),center(roi([2:4 end]),3),[],'g','filled');
scatter3(center(roi([5:9]),1),center(roi([5:9]),2),center(roi([5:9]),3),[],'y','filled');
arrayfun(@(i) text(center(roi(i),1),center(roi(i),2),center(roi(i),3),num2str(roi(i))),1:length(roi),'un',0);
%%
% spike_focus = arrayfun(@(i) Spike_X_EstTrace(roi,i+[-6:5]),col2row(focus_idx_conv_start,1),'un',0);
% spike_focus = arrayfun(@(i) CalTrace3_original(roi,i+[-6:5]) - mean(CalTrace3_original(roi,i+[-10:-5]),2),col2row(focus_idx_conv_start,1),'un',0);
spike_focus = arrayfun(@(i) CalTrace3_original(roi,i+[-6:5]) - mean(CalTrace3_original(roi,i+[-10:-5]),2),col2row(focus_idx_conv_start,1),'un',0);
spike_focus = cat(3,spike_focus{:});%region x time x trial
figure,hold on;color = brewermap(10,'Paired');
% arrayfun(@(i) plot(-0.6:0.1:0.5,median(spike_focus(i,:,:),3),'Color',color(i,:)),1:length(roi),'un',0);
[hl,hp] = arrayfun(@(i) boundedline(-0.6:0.1:0.5,median(spike_focus(i,:,:),3),permute(std(spike_focus(i,:,:),[],3)/sqrt(size(spike_focus,3)),[2 1 3]),...
    'cmap',color(i,:),'alpha'),1:length(roi),'un',0);
legend(cat(1,hl{:}),cellfun(@(x) num2str(x),num2cell(roi),'un',0));
%%
%granger causality
X = spike_focus;%roi x time x trial
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
%%
