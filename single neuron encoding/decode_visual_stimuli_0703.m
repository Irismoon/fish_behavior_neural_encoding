%0703 receptive field
%%
cSession = '200703';cFish = '1';
load('Coherence3','A3');
%%
load(fullfile(getpath('neural activity',cSession,cFish),'spikes.mat'), 'spikes_OASIS');
spikes_OASIS = smoothdata(spikes_OASIS,2,'movmean',5);
[nRegion,nTime] = size(spikes_OASIS);
%%
load(fullfile(getpath('behavior',cSession,cFish),'low_video_analysis_result'),'param_head_angle_all','param_head_dist_all','param_speed_all','no_param');
no_param_fluo = row2col(ismember(1:5:length(no_param),find(no_param==0)),1);
%%
%no param
param_head_angle_fluo = [piecewisefn(@mean,param_head_angle_all,'edge',1:5:length(param_head_angle_all));param_head_angle_all(end)];
param_head_dist_fluo = [piecewisefn(@mean,param_head_dist_all,'edge',1:5:length(param_head_dist_all));param_head_dist_all(end)];
param_speed_fluo = [piecewisefn(@mean,param_speed_all,'edge',1:5:length(param_speed_all));param_speed_asll(end)];
%param too far
toofar = row2col(param_head_dist_fluo>60,1);
%not moving
notmoving = row2col(param_speed_fluo<9,1);
%
idx = ~(no_param_fluo | toofar | notmoving);
[param_head_angle_visuo,param_head_dist_visuo,param_speed_visuo] = samfnmultvar(@(x) x(idx),param_head_angle_fluo,param_head_dist_fluo,param_speed_fluo);
spikes_visuo = spikes_OASIS(:,idx);
%%
%plain anova
disp('running current section...');
param_focus = param_head_angle_visuo;
delay = 1:8;
p = zeros(nRegion,3,length(delay));
rsquare = zeros(nRegion,length(delay));
rmse = zeros(length(delay),1);
for k=1:length(delay)
    for i=1:nRegion
        [p(i,:,k),tb1] = anovan(spikes_visuo(i,k:end),{param_head_angle_visuo(1:end-k+1),param_head_dist_visuo(1:end-k+1)},'continuous',[1 2],'model','interaction','display','off');
        rsquare(i,k) = (tb1{end,2}-tb1{end-1,2})/tb1{end,2};
    end
    selcRegion = find(sum(p(:,:,k)<0.001/nRegion,2)==3);
    for cv = 1:5
        X = spikes_visuo(selcRegion,k:end);y = param_head_angle_visuo(1:end-k+1);
        testidx = randperm(size(X,2),floor(size(X,2)*0.2));
        mask = true(size(X,2),1);mask(testidx) = false;
        mdl = fitlm(X(:,mask)',y(mask));
        ypred = predict(mdl,X(:,~mask)');
        rmse(k) = rmse(k)+mean((ypred - y(~mask)).^2)/var(y(~mask));
    end
end
rmse = rmse/5;
% roi = find(sum(p<0.001/nRegion,2)==3);
figure,histogram(rsquare);
disp('this section done!');
%%
%regress neural activity to a pixel vector
%create the visual space
disp('starting this section...');
[polarcoor,distgrid,anglegrid] = pixelCoor(param_head_angle_visuo,param_head_dist_visuo,param_speed_visuo);%pixel x timepoints
%regression
delay_est = 1;
D = zeros(size(polarcoor,1),length(selcRegion));
rsquare = zeros(size(polarcoor,1),1);
for ipixel = 1:size(polarcoor,1)
    mdl = fitlm(spikes_visuo(selcRegion,:)',polarcoor(ipixel,:)');
    D(ipixel,:) = mdl.Coefficients.Estimate(2:end);
    rsquare(ipixel) = mdl.Rsquared.Ordinary;
end
figure,imagesc(reshape(rsquare,length(distgrid),length(anglegrid)));
set(gca,'XTick',5:5:length(anglegrid),'XTickLabel',anglegrid(5:5:end));
set(gca,'YTick',2:2:length(distgrid),'YTickLabel',distgrid(2:2:end));
xlabel('angle');ylabel('dist');c = colorbar;title(c,'rsquare');
title('regression result: how much variance explained at each pixel');
figure,
imagesc(D);xlabel('selected region');ylabel('# pixel');colorbar;
title('weight matrix');
disp('this section done!');
%%
%find the best fitted pixel and corresponding weight vector, look if the
%region's receptive field is that pixel
[~,I] = max(rsquare);
[~,roi] = maxk(abs(D(I,:)),10);
figure,
scatter3(center(:,1),center(:,2),center(:,3),10,'filled');hold on;
scatter3(center(selcRegion(roi),1),center(selcRegion(roi),2),center(selcRegion(roi),3),30,'filled');
%%
%regress dist and angle using the selected regions
% mdl = fitlm(spikes_OASIS(roi,:)',param_head_angle_fluo);
%pure linear model is not enough to consider the delay between calcium
%signal and visual stimulus
param_focus = param_head_dist_visuo;
feature =  delaymatrix(spikes_visuo(roi,:)',3);
mdl = fitlm(feature,param_focus);
figure,scatter(param_focus,mdl.Fitted);axis equal;hold on;fplot(@(x) x);
title(['regression to param head angle: rsquare=' num2str(mdl.Rsquared.Ordinary)]);
%%
%preferred angle and dist of selected region
%split the visual space into blocks
anglelim = max(abs(param_head_angle_visuo));
anglegrid = -anglelim:5:anglelim;
anglegrid = [anglegrid anglegrid(end)+5];
distlim = [min(param_head_dist_visuo) max(param_head_dist_visuo)];
distgrid = row2col(log10(distlim(1)):0.1:log10(distlim(2)),1);
distgrid = [distgrid;distgrid(end)+0.1];
label_dist = discretize(log10(param_head_dist_visuo),distgrid);
label_angle = discretize(param_head_angle_visuo,anglegrid);
polargrid = [label_dist label_angle];
[N,distedge,angleedge] = histcounts2(log10(param_head_dist_visuo),param_head_angle_visuo,distgrid,anglegrid);
[row,col] = find(N);
spike_selected = spikes_visuo(roi,:)';%time x region
fr = arrayfun(@(i) median(spike_selected(polargrid(:,1)==row(i)&polargrid(:,2)==col(i),1)),1:length(row));
fr = fr/sum(fr);
prob_fr = nan(size(N));
for i=1:length(row)
    prob_fr(row(i),col(i)) = fr(i);
end
%%
load('tail_swing','left_tail_swing','right_tail_swing','sum_curv');
tail_swing = -1*left_tail_swing + right_tail_swing;
tail_swing_fluo = tail_swing_fluo(1:5:length(tail_swing));
sum_curv_fluo = [piecewisefn(@mean,sum_curv,'edge',1:5:length(sum_curv));sum_curv(end)];
swing_direct = swing2whichside(tail_swing,bout_idx);%-1:left,1:right,2:forward
swing_direct_fluo = swing_direct(1:5:length(swing_direct));
%%
%motor regression
delay = [1 2 3 4 5 6 7 8];
p = zeros(nRegion,3,length(delay));
rsquare = zeros(nRegion,length(delay));
for k=delay
    for i=1:nRegion
        [p(i,:,k),tb1] = anovan(spikes_OASIS(i,(k+1):end)',{swing_direct_fluo(1:end-k),sum_curv_fluo(1:end-k)},'continuous',[2],'model','interaction','display','off');
        rsquare(i,k) = 1 - tb1{end-1,2}/tb1{end,2};
    end
end
motorTunedRegion = squeeze(sum(p<0.001/nRegion/length(delay),2)==3);
nTunedRegion = sum(motorTunedRegion);%region x delay
figure,plot(delay,nTunedRegion);
xlabel('delay');ylabel('# regions tuned to motor');
figure,
scatter3(center(:,1),center(:,2),center(:,3),10,'filled');
hold on;
idx = 5;
scatter3(center(motorTunedRegion(:,idx),1),center(motorTunedRegion(:,idx),2),center(motorTunedRegion(:,idx),3),20,'filled');
%%
%fit to motor
rsquare = zeros(length(delay),1);
for i=delay(1:3)
    feature = spikes_OASIS(motorTunedRegion(:,i),:)';%time x region
    % feature = spikes_OASIS(:,:)';
    testidx = randperm(nTime-i,floor((nTime-i)/5));
    mask = true(nTime-i,1);mask(testidx)=false;
    X = feature(i+1:end,:);y = sum_curv_fluo(1:end-i);
    mdl = fitlm(X(mask,:),y(mask));
    ypred = predict(mdl,X(~mask,:));
    rsquare(i) = mean((y(~mask) - ypred).^2)/var(ypred);
end
%%
%check the receptive field of selected regions 
%only consider angle and no delay
%the population average firing rate at each angle bin 
anglelim = [min(param_head_angle_visuo) max(param_head_angle_visuo)];
anglegrid = 5*[-19:2:19];
label_angle = discretize(param_head_angle_visuo,anglegrid);
N = histcounts(param_head_angle_visuo,anglegrid);
spikes_anglebin = arrayfun(@(i) median(spikes_visuo(:,label_angle==i),2),1:length(anglegrid)-1,'un',0);
std_spikes_anglebin = arrayfun(@(i) std(spikes_visuo(:,label_angle==i),[],2),1:length(anglegrid)-1,'un',0);
del_mask = find(N<50);
label_angle(ismember(label_angle,del_mask)) = [];
anglegrid = arrayfun(@(i) mean(anglegrid(i:i+1)),1:length(anglegrid)-1);
anglegrid(del_mask) = [];
spikes_anglebin(del_mask) = [];
std_spikes_anglebin(del_mask) = [];
spikes_anglebin = cat(2,spikes_anglebin{:});
std_spikes_anglebin = cat(2,std_spikes_anglebin{:});
figure,imagesc(spikes_anglebin);
set(gca,'XTick',1:1:length(anglegrid),'XTickLabel',floor(anglegrid(1:1:length(anglegrid))));
[fr_angle,prefangle] = max(spikes_anglebin,[],2);
figure,histogram(prefangle);
Npref = accumarray(prefangle,1);
set(gca,'XTick',1:1:length(anglegrid),'XTickLabel',floor(anglegrid(1:1:length(anglegrid))));
xlabel('preferred angle');ylabel('# regions');
figure,
imagesc(corrcoef(spikes_anglebin));colorbar;
set(gca,'XTick',1:1:length(anglegrid),'XTickLabel',floor(anglegrid(1:1:length(anglegrid))));
set(gca,'YTick',1:1:length(anglegrid),'YTickLabel',floor(anglegrid(1:1:length(anglegrid))));
color = colormap(parula(length(anglegrid)));
figure,hold on;
idx = unique(prefangle);
for i=1:length(idx)
    scatter3(center(prefangle==idx(i),1),center(prefangle==idx(i),2),center(prefangle==idx(i),3),10,color(i,:),'filled');
end
cbar = colorbar;cbar.Ticks = linspace(0,1,length(idx));cbar.TickLabels = floor(anglegrid);
xlabel('x');ylabel('y');zlabel('z');