% from the selected good bouts, decode the prey location in the visual field in a time-varying way.
csessionID = '200108';cfishID = '';
numT = 11;
%%
load(fullfile(getpath('behavior',csessionID,cfishID),'tail_swing'),'good_bout_idx','angledata');
load(fullfile(getpath('behavior',csessionID,cfishID),'align_with_fluo'),'align_with_fluo_low');
load(fullfile(getpath('behavior',csessionID,cfishID),'low_video_analysis_result'),'no_param','param_head_angle_all','param_head_dist_all','param_speed_all');
load(fullfile(getpath('behavior',csessionID,cfishID),'leftright_eyes_to_param_angle_dist'),'left_eye_to_param_angle','right_eye_to_param_angle','left_eye_to_param_angle_visible', 'right_eye_to_param_angle_visible');
load(fullfile(getpath('neural activity',csessionID,cfishID),'spikes_new'),'spikes_OASIS');
load(fullfile(getpath('neural activity',csessionID,cfishID),'Coherence3_correct'),'center');
spikes_OASIS(:,1651:1779) = [];
Spike_X_EstTrace = smoothdata(spikes_OASIS,2,"gaussian",5);
%0108 good bout idx was already aligned to fluo data because its high and
%low fluo is very different from each other
angledata = angledata(align_with_fluo_low==1,:);
angledata_fluo = angledata(1:5:length(angledata),:);
[param_head_angle_all,param_head_dist_all,param_speed_all,no_param] = samfnmultvar(@(x) x(align_with_fluo_low==1),param_head_angle_all,param_head_dist_all,param_speed_all,no_param);
[left_eye_to_param_angle,right_eye_to_param_angle,left_eye_to_param_angle_visible,right_eye_to_param_angle_visible] = ...
    samfnmultvar(@(x) x(align_with_fluo_low==1),left_eye_to_param_angle,right_eye_to_param_angle,left_eye_to_param_angle_visible,right_eye_to_param_angle_visible);
[param_head_angle_fluo,param_head_dist_fluo,param_speed_fluo,no_param_fluo] = samfnmultvar(@(x) x(1:5:length(param_head_angle_all)),param_head_angle_all,param_head_dist_all,param_speed_all,no_param);
[lefteye_angle_fluo,righteye_angle_fluo,lefteye_angle_vis_fluo,righteye_angle_vis_fluo] = samfnmultvar(@(x) x(1:5:length(param_head_angle_all)),...
    left_eye_to_param_angle,right_eye_to_param_angle,left_eye_to_param_angle_visible,right_eye_to_param_angle_visible);
param_head_angle_hunting = mean(reshape(param_head_angle_all(reshape((good_bout_idx(:,1)-[1:20])',[],1)),20,[]),1);
mask = 1:5:nnz(align_with_fluo_low);
[~,hunting_idx_start] = min(abs(row2col(good_bout_idx(:,1),1) - col2row(mask)),[],2);
[~,hunting_idx_end] = min(abs(row2col(good_bout_idx(:,2),1) - col2row(mask)),[],2);
% fluo_idx = randperm(size(Spike_X_EstTrace,2)-numT,nnz(lia));
leftIdx = param_head_angle_hunting>0;
rightIdx = param_head_angle_hunting<0;
[numRegion,numWholeTime] = size(Spike_X_EstTrace);
numTrial = max(nnz(leftIdx),nnz(rightIdx));
if nnz(leftIdx)>=nnz(rightIdx)
    moreIdx = leftIdx;
    lessIdx = rightIdx;
    disp('left right');
    direction = {'left','right'};
else
    moreIdx = rightIdx;
    lessIdx = leftIdx;
    disp('right left');
    direction = {'right','left'};
end
%%
spk = zeros(length(hunting_idx_start),size(Spike_X_EstTrace,1),numT);%bouts x region x t
for i = 1:length(hunting_idx_start)
    spk(i,:,:) = Spike_X_EstTrace(:,hunting_idx_start(i)+(-floor(numT/2):floor(numT/2)));%region x 5
end
m_spk_left = squeeze(mean(spk(param_head_angle_hunting>0,:,:),1));
m_spk_right = squeeze(mean(spk(param_head_angle_hunting<0,:,:),1));
%map them on MIP
generate_MIP(m_spk_left+0.01*max(m_spk_left(:)),csessionID,cfishID,'_meanFr_left');
generate_MIP(m_spk_right+0.01*max(m_spk_right(:)),csessionID,cfishID,'_meanFr_right');
%%
missing_data = mean(spk(lessIdx,:,:));%1 x region x time
missing_data = repmat(missing_data,nnz(moreIdx)-nnz(lessIdx),1,1);
spk_balanced = cat(1,spk(moreIdx,:,:),spk(lessIdx,:,:),missing_data);%
spk_balanced = permute(reshape(spk_balanced,nnz(moreIdx),2,numRegion,numT),[3,2,4,1]);%region x direction x time x trial
% [W,V,whichMarg] = drtoolbox.techniques.dpca(mean(spk_balanced,4),5,'combinedParams',{{1},{2},{1,2}});
spk_balanced_ave = mean(spk_balanced,4);%neuron x direction x time
X = spk_balanced_ave(:,:);
% X = X - mean(X,2);
[W,explVar,V ] = svd(X,'econ');
explVar = diag(explVar)/sum(diag(explVar));
W = W(:,1:numT*2);

% drtoolbox.techniques.dpca_plot(permute(spk_balanced,[1 4 2 3]), W, W, @drtoolbox.techniques.dpca_plot_default);
figure,
hold on;plot(explVar/sum(explVar),'-*');axis tight;xlabel('comp');ylabel('explained variance');
figure,
plot(cumsum(explVar/sum(explVar)),'-*');axis tight;xlabel('comp');ylabel('cumulative explained variance');
stdx = permute(std(pagemtimes(W',spk_balanced),[],4)./sqrt([nnz(moreIdx),nnz(lessIdx)]),[3 2 1]);%
multchnPlot(permute(pagemtimes(W',spk_balanced_ave),[3 2 1]),0.1*(-floor(numT/2):floor(numT/2)),[0],'boundwidth',stdx,...
    'plotstyle','boundplot','chnname',arrayfun(@(i) ['comp ' num2str(i)],1:numT*2,'un',0));
sgtitle('signal projection on each component');
xlabel('t/s');ylabel('firing rate');
legend(direction);

generate_MIP(abs(W(:,2)),csessionID,cfishID,'_2ndPC_average_hunting');

%%
%investigate the raw trace of region #356 which contains most 2nd PC time
%series
% load('tail_swing','bout_idx');
[~,bout_idx_aligned_start] = min(abs(row2col(bout_idx(:,1),1) - mask),[],2);
[~,bout_idx_aligned_end] = min(abs(row2col(bout_idx(:,2),1) - mask),[],2);
spk_focus = Spike_X_EstTrace(356,:);
bout_idx_aligned_wholecourse = arrayfun(@(i) (bout_idx_aligned_start(i):bout_idx_aligned_end(i))',1:length(bout_idx),'un',0);
bout_idx_aligned_wholecourse = cat(1,bout_idx_aligned_wholecourse{:});
nonbout_idx = setdiff(1:numWholeTime,bout_idx_aligned_wholecourse);
figure,
group = zeros(numWholeTime,1);group(bout_idx_aligned_wholecourse) = 1;group(nonbout_idx)=2;
subplot(2,1,1)
p = ranksum(spk_focus(bout_idx_aligned_wholecourse),spk_focus(nonbout_idx));
boxplot(spk_focus,group);
sigstar([1 2],p);
set(gca,'XTickLabel',{'bout','nonbout'});
title(['Comparison of neural activity between bout and no bouts: p=' num2str(p)]);
%control
spk_focus_ctrl = spk_focus(randperm(length(spk_focus)));
%and rerun the code above
subplot(2,1,2),
p = ranksum(spk_focus_ctrl(bout_idx_aligned_wholecourse),spk_focus_ctrl(nonbout_idx));
boxplot(spk_focus_ctrl,group);
sigstar([1 2],p);
set(gca,'XTickLabel',{'bout','nonbout'});
title(['control: p=' num2str(p)]);
%hunting bouts and nonhunting bouts
hunting_idx_wholecourse = arrayfun(@(i) (hunting_idx_start(i):hunting_idx_end(i))',1:length(hunting_idx_start),'un',0);
hunting_idx_wholecourse = cat(1,hunting_idx_wholecourse{:});
p = ranksum(spk_focus(hunting_idx_wholecourse),spk_focus(setdiff(bout_idx_aligned_wholecourse,hunting_idx_wholecourse)));
figure,
subplot(2,1,1),
group = 2*ones(length(bout_idx_aligned_wholecourse),1);
group(ismember(bout_idx_aligned_wholecourse,hunting_idx_wholecourse)) = 1;
boxplot(spk_focus(bout_idx_aligned_wholecourse),group);
hold on;
sigstar([1 2],p);
set(gca,'XTickLabel',{'hunting bouts','nonhunting bouts'});
subplot(2,1,2),
group = 2*ones(length(bout_idx_aligned_wholecourse),1);
group(ismember(bout_idx_aligned_wholecourse,hunting_idx_wholecourse)) = 1;
boxplot(spk_focus(bout_idx_aligned_wholecourse),group(randperm(length(bout_idx_aligned_wholecourse))));
hold on;
idxytmp=randperm(length(bout_idx_aligned_wholecourse),length(hunting_idx_wholecourse));
p = ranksum(spk_focus(bout_idx_aligned_wholecourse(idxytmp)),spk_focus(bout_idx_aligned_wholecourse(setdiff(1:length(bout_idx_aligned_wholecourse),idxytmp))));
sigstar([1 2],p);
set(gca,'XTickLabel',{'hunting bouts','nonhunting bouts'});
title('control');
sgtitle('comparison between hunting and nonhuting bouts');
set(gca,'XTickLabel',{'hunting','nonhunting'});
%compare nonhunting bouts and nonbouts
figure,
x1 = spk_focus(setdiff(1:numWholeTime,bout_idx_aligned_wholecourse))';
x2 = spk_focus(setdiff(bout_idx_aligned_wholecourse,hunting_idx_wholecourse))';
x = cat(1,x2,x1);
group = [ones(length(x2),1);2*ones(length(x1),1)];
p = ranksum(x1,x2);
subplot(2,1,1),
boxplot(x,group);hold on;
sigstar([1 2],p);
set(gca,'XTickLabel',{'nonhunting bouts','nonbouts'});
subplot(2,1,2),
randidx = randperm(length(x));
p = ranksum(x(randidx(1:length(x1))),x(randidx((length(x1)+1):end)));
boxplot(x,group(randidx));hold on;
sigstar([1 2],p);
sgtitle('comparison between nonhunting bouts and nonbouts');
%compare left and right bouts
figure,
x1 = arrayfun(@(i) (hunting_idx_start(i):hunting_idx_end(i))',find(param_head_angle_hunting>0),'un',0);
x1 = spk_focus(cat(1,x1{:}))';%left
x2 = arrayfun(@(i) (hunting_idx_start(i):hunting_idx_end(i))',find(param_head_angle_hunting<0),'un',0);
x2 = spk_focus(cat(1,x2{:}))';%right
x = cat(1,x2,x1);
group = [ones(length(x2),1);2*ones(length(x1),1)];
p = ranksum(x1,x2);
subplot(2,1,1),
boxplot(x,group);hold on;
sigstar([1 2],p);
set(gca,'XTickLabel',{'right','left'});
subplot(2,1,2),
randidx = randperm(length(x));
p = ranksum(x(randidx(1:length(x1))),x(randidx((length(x1)+1):end)));
boxplot(x,group(randidx));hold on;
sigstar([1 2],p);
sgtitle('comparison between left and right hunting bouts');
%%
%compare left and right vision
x1 = spk_focus(hunting_idx_wholecourse(param_head_angle_fluo(hunting_idx_wholecourse)>0))';
x2 = spk_focus(hunting_idx_wholecourse(param_head_angle_fluo(hunting_idx_wholecourse)<0))';
px = ranksum(x1,x2);
simul_bout_idx_wholecourse = setdiff(bout_idx_aligned_wholecourse,hunting_idx_wholecourse);
y1 = spk_focus(simul_bout_idx_wholecourse(param_head_angle_fluo(simul_bout_idx_wholecourse)>0))';
y2 = spk_focus(simul_bout_idx_wholecourse(param_head_angle_fluo(simul_bout_idx_wholecourse)<0))';
y = cat(1,y1,y2);
py = ranksum(y1,y2);
rest_idx_wholecourse = setdiff(1:numWholeTime,bout_idx_aligned_wholecourse);
z1 = spk_focus(rest_idx_wholecourse(param_head_angle_fluo(rest_idx_wholecourse)>0))';
z2 = spk_focus(rest_idx_wholecourse(param_head_angle_fluo(rest_idx_wholecourse)<0))';
pz = ranksum(z1,z2);
data = {x1,x2;y1,y2;z1,z2};
h = boxplot2(data);
set(gca,'XTick',[1 2 3],'XTickLabel',{'hunting','simultaneous','resting'});
sigstar([h.lwhis(1,1).XData(1) h.lwhis(2,1).XData(1)],px);
sigstar([h.lwhis(1,2).XData(1) h.lwhis(2,2).XData(1)],py);
sigstar([h.lwhis(1,3).XData(1) h.lwhis(2,3).XData(1)],pz);
sigstar([h.lwhis(1,1).XData(1) h.lwhis(1,2).XData(1)],ranksum(x1,y1));
sigstar([h.lwhis(1,2).XData(1) h.lwhis(1,3).XData(1)],ranksum(z1,y1));
sigstar([h.lwhis(2,1).XData(1) h.lwhis(2,2).XData(1)],ranksum(x2,y2));
sigstar([h.lwhis(2,2).XData(1) h.lwhis(2,3).XData(1)],ranksum(z2,y2));
title('Left Right visual encoding');
%%
%if the left/right difference comes from the distance/angle
idx_no_param = find(param_head_dist_fluo>300);
hunting_idx_wholecourse(ismember(hunting_idx_wholecourse,idx_no_param)) = [];
simul_bout_idx_wholecourse(ismember(simul_bout_idx_wholecourse,idx_no_param)) = [];
rest_idx_wholecourse(ismember(rest_idx_wholecourse,idx_no_param)) = [];
dx = param_head_dist_fluo(hunting_idx_wholecourse);
dy = param_head_dist_fluo(simul_bout_idx_wholecourse);
dz = param_head_dist_fluo(rest_idx_wholecourse);
data = {dx,dy,dz};
figure,
subplot(1,2,1),
boxplot2(data');
pxy = ranksum(dx,dy);pyz = ranksum(dy,dz);
hold on;
sigstar([1 2],pxy);sigstar([2 3],pyz);
ylabel('param distance');
set(gca,'XTick',[1 2 3],'XTickLabel',{'hunitng','simultaneous','resting'});
subplot(1,2,2),
anglex = param_head_angle_fluo(hunting_idx_wholecourse);
angley = param_head_angle_fluo(simul_bout_idx_wholecourse);
anglez = param_head_angle_fluo(rest_idx_wholecourse);
data = {anglex,angley,anglez};
boxplot2(data');

%%
%behavior state x left/right
behav_var = param_head_dist_fluo;
idx_behav7direct = cell(3,2);
idx_behav7direct{1,1} = hunting_idx_wholecourse(param_head_angle_fluo(hunting_idx_wholecourse)>0);
idx_behav7direct{1,2} = hunting_idx_wholecourse(param_head_angle_fluo(hunting_idx_wholecourse)<0);
idx_behav7direct{2,1} = simul_bout_idx_wholecourse(param_head_angle_fluo(simul_bout_idx_wholecourse)>0);
idx_behav7direct{2,2} = simul_bout_idx_wholecourse(param_head_angle_fluo(simul_bout_idx_wholecourse)<0);
idx_behav7direct{3,1} = rest_idx_wholecourse(param_head_angle_fluo(rest_idx_wholecourse)>0);
idx_behav7direct{3,2} = rest_idx_wholecourse(param_head_angle_fluo(rest_idx_wholecourse)<0);
spk_behav7direct = cellfun(@(x) spk_focus(x)',idx_behav7direct,'un',0);
dist_behav7direct = cellfun(@(x) behav_var(x),idx_behav7direct,'un',0);
xlim = [min(cellfun(@min,dist_behav7direct),[],'all') max(cellfun(@max,dist_behav7direct),[],'all')];
ylim = [min(cellfun(@min,spk_behav7direct),[],'all') max(cellfun(@max,spk_behav7direct),[],'all')];
figure,
coeff = zeros(3,2,2);GOF = zeros(3,2);
Edge = cell(3,2);M_SPK = cell(3,2);Weight = cell(3,2);
for i=1:3
    for j=1:2
        subplot(3,2,(i-1)*2+j),
%         scatter(dist_behav7direct{i,j},spk_behav7direct{i,j},10,'filled');hold on;
        [Y,edges] = discretize(dist_behav7direct{i,j},0:1:max(dist_behav7direct{i,j}));
        m_spk = arrayfun(@(k) quantile(spk_behav7direct{i,j}(Y==k),0.75),1:(length(edges)-1));
%         m_spk(isnan(m_spk)) = 0;        
        [N,xedge] = histcounts(dist_behav7direct{i,j},edges);
        m_spk(N==0) = [];
        edges(N==0) = [];
        xedge(N==0) = [];
        N(N==0) = [];
        Edge{i,j} = edges(1:end-1);
        M_SPK{i,j} = m_spk;
        Weight{i,j} = N;
%         scatter(edges(1:end-1),m_spk,N);
        [f,gof] = fit(edges(1:end-1)',m_spk','poly1','weights',N/sum(N));
        coeff(i,j,:) = [f.p1 f.p2];% f.p3];
        GOF(i,j) = gof.rsquare;
%         yest = f.p1*edges(1:end-1)+f.p2;
%         hold on;plot(edges(1:end-1),yest);
%         text(edges(1),yest(1),num2str(gof.rsquare));
%         p_obs = f.p1;
%         p_shuffle = zeros(1000,1);
        %shuffle test
%         for iperm = 1:1000
%             randidx = randperm(length(m_spk));
%             f = fit(edges(1:end-1)',m_spk(randidx)','poly1','weights',N/sum(N));
%             p_shuffle(iperm) = f.p1;
%         end
%         h = ttest(p_shuffle)
%         title(['normal? :' num2str(h) 'p value:' num2str(mean(p_shuffle<p_obs),3)]);
%         [f,gof] = fit(edges(1:end-1)',m_spk','exp1','weights',N/sum(N));
%         yest =f.a*exp(f.b*edges(1:end-1));
%         hold on;plot(edges(1:end-1),yest);
%         text(edges(end-1),yest(end-1),num2str(gof.rsquare));
%         yyaxis right;
%         plot(xedge(1:end-1),N);
%         axis tight;
%         [N,xedge,yedge] = histcounts2(dist_behav7direct{i,j},spk_behav7direct{i,j},'BinWidth',[2 0.01],'Normalization','probability','XBinLimits',xlim,'YbinLimits',ylim);
%         N = N./sum(N,2);
%         imagesc(N);
%         set(gca,'XTick',1:length(yedge),'XTickLabel',yedge,'YTick',1:length(xedge),'YTickLabel',xedge);
    end
end
xlabel('distance');ylabel('probability');yyaxis left;ylabel('median firing rate');

%demix the distance and behavior state condition
figure,
color = rgb({'darkgreen','magenta','navy'});
hold on,
lateral = 2;
for i=1:3
   scatter(Edge{i,lateral},M_SPK{i,lateral},Weight{i,lateral},color(i,:),'filled');
%    h(i) = plot(Edge{i,lateral},coeff(i,lateral,1)*(Edge{i,lateral}.^2)+coeff(i,lateral,2)*Edge{i,lateral}+coeff(i,lateral,3),'Color',color(i,:));
    h(i) = plot(Edge{i,lateral},coeff(i,lateral,1)*Edge{i,lateral}+coeff(i,lateral,2),'Color',color(i,:));
    text(Edge{i,lateral}(1),M_SPK{i,lateral}(1),num2str(GOF(i,lateral)),'Color',color(i,:));
end
legend(h,{'hunting','simultaneous','resting'});
%%
 %select the timeponits of simultaneous state which has the same distance
 %with hunting, and compare their firing rates
 idxxtmp = find(anglex<0);
 dx_right = dx(idxxtmp);
dlim = [min(dx_right) quantile(dx_right,0.75)];
idxytmp = find(dy<dlim(2)&dy>dlim(1)&angley<0);
idxytmp = idxytmp(randperm(length(idxytmp),length(dx_right)));
dy_sel = dy(idxytmp);
spk_y_sel = spk_focus(simul_bout_idx_wholecourse(idxytmp));
[a,b] = min(abs(dy_sel-dx_right'),[],2);
thld = 0.2;
dy_sam = dy_sel(a<thld);
pair = [find(a<thld) b(a<thld)];
figure,
subplot(3,1,1),
boxplot([dy_sel dx_right],'notch','on');
p = ranksum(dx_right,dy_sel);
hold on;sigstar([1 2],p);
arrayfun(@(i) plot([1+randn*0.03 2+randn*0.03],[dy_sel(pair(i,1)) dx_right(pair(i,2))],'--ok','MarkerSize',3),1:size(pair,1),'un',0);
ylabel('distance');
subplot(3,1,2),
spk_x = spk_focus(hunting_idx_wholecourse);
spk_x = spk_x(idxxtmp);
boxplot([spk_y_sel spk_x],'notch','on');
p = ranksum(spk_y_sel,spk_x);
hold on;
sigstar([1 2],p);
arrayfun(@(i) plot([1+randn*0.03 2+randn*0.03],[spk_y_sel(pair(i,1)) spk_x(pair(i,2))],'--ok','MarkerSize',3),1:size(pair,1),'un',0);
set(gca,'XTickLabel',{'randomly selected simultaneous','hunting'});
ylabel('Firing Rate');
sgtitle('Compare firing rate at same distance but diff behavior');
subplot(3,1,3),
boxplot([angley(idxytmp) anglex(idxxtmp)]);
hold on;
p = ranksum(angley(idxytmp),anglex(idxxtmp));
sigstar([1 2],p);
incre = arrayfun(@(i) spk_x(pair(i,2))-spk_y_sel(pair(i,1)),1:size(pair,1));
cmap = colormap(jet(length(incre)));
[~,I] = sort(incre,'ascend');
cmap = cmap(I,:);
angley_sel = angley(idxytmp);
arrayfun(@(i) plot([1+randn*0.03 2+randn*0.03],[angley_sel(pair(i,1)) anglex(idxxtmp(pair(i,2)))],'o-','MarkerSize',3,'Color',cmap(i,:)),1:size(pair,1),'un',0);
ylabel('angle');
%%
figure,
% polarscatter(deg2rad(anglex),ones(length(anglex),1)+randn(length(anglex),1)*0.05,5,spk_focus(hunting_idx_wholecourse),'filled');colormap('jet');
% polarscatter(deg2rad(angley),ones(length(angley),1)+randn(length(angley),1)*0.05,5,spk_focus(simul_bout_idx_wholecourse),'filled');colormap('jet');
% polarscatter(deg2rad(anglez),ones(length(anglez),1)+randn(length(anglez),1)*0.05,5,spk_focus(rest_idx_wholecourse),'filled');colormap('jet');
spkx = spk_focus(hunting_idx_wholecourse);
spky = spk_focus(simul_bout_idx_wholecourse);
spkz = spk_focus(rest_idx_wholecourse);
spk3 = cat(1,spkx,spky,spkz);
angle3 = cat(1,anglex,angley,anglez);
spkx_right = spkx(anglex<0);
spky_right = spky(angley<0);
spkz_right = spkz(anglez<0);
spk3_right = cat(1,spkx_right,spky_right,spkz_right);
angle3_right = cat(1,anglex(anglex<0),angley(angley<0),anglez(anglez<0));
[N,Xedges,Yedges] = histcounts2(spk3,angle3,'BinWidth',[0.02 2]);
% imagesc(N./sum(N,1));
% imagesc(N);
spkm = arrayfun(@(i) quantile(spk3(angle3<Yedges(i+1)&angle3>Yedges(i)),.75),1:(length(Yedges)-1));
scatter(Yedges(1:end-1),spkm,10,'filled');
% xlabel('angle');ylabel('firingRate');
% set(gca,'XTick',5:5:length(Xedges),'XTickLabel',Xedges(5:5:length(Xedges)),'YTick',5:5:length(Yedges),'YTickLabel',Yedges(5:5:length(Yedges)));
%%
%comprehensively consider distance, angle at right side, and behavior
%states and see which is important for the neural activity
dx_right = dx(anglex<0);
dy_right = dy(angley<0);
dz_right = dz(anglez<0);
d3_right = cat(1,dx_right,dy_right,dz_right);
behavior_states = [repmat({'hunting'},length(dx_right),1);repmat({'simultaneous'},length(dy_right),1);repmat({'rest'},length(dz_right),1)];
[p,tb1] = anovan(spk3_right,{behavior_states,d3_right,angle3_right},'continuous',[2 3],'model','interaction','varnames',{'behavior states','distance','angle'});
%%
%consider the relative angle between prey and eye rather than that between
%the prey and head
leftanglex = left_eye_to_param_angle_fluo(hunting_idx_wholecourse);
rightanglex = right_eye_to_param_angle_fluo(hunting_idx_wholecourse);
leftangley = left_eye_to_param_angle_fluo(simul_bout_idx_wholecourse);
rightangley = right_eye_to_param_angle_fluo(simul_bout_idx_wholecourse);
leftanglez = left_eye_to_param_angle_fluo(rest_idx_wholecourse);
rightanglez = right_eye_to_param_angle_fluo(rest_idx_wholecourse);
leftanglex_right = leftanglex(anglex<0);
rightanglex_right = rightanglex(anglex<0);
leftangley_right = leftangley(angley<0);
rightangley_right = rightangley(angley<0);
leftanglez_right = leftanglez(anglez<0);
rightanglez_right = rightanglez(anglez<0);
leftangle3_right = cat(1,leftanglex_right,leftangley_right,leftanglez_right);
rightangle3_right = cat(1,rightanglex_right,rightangley_right,rightanglez_right);
[p,tb2] = anovan(spk3_right,{behavior_states,d3_right,leftangle3_right,rightangle3_right},'continuous',[2 3 4],'model','interaction','varnames',{'behavior states','distance','leftEyeAngle','rightEyeAngle'});
%%
param_speedx = param_speed_fluo(hunting_idx_wholecourse);
param_speedy = param_speed_fluo(simul_bout_idx_wholecourse);
param_speedz = param_speed_fluo(rest_idx_wholecourse);
param_speedx_right = param_speedx(anglex<0);
param_speedy_right = param_speedy(angley<0);
param_speedz_right = param_speedz(anglez<0);
param_speed3_right = cat(1,param_speedx_right,param_speedy_right,param_speedz_right);
[p,tb2] = anovan(spk3_right,{behavior_states,d3_right,leftangle3_right,rightangle3_right,param_speed3_right},...
    'continuous',[2 3 4 5],'model','interaction','varnames',...
    {'behavior states','distance','leftEyeAngle','rightEyeAngle','speed'});
%%

%%
%the movement speed and angle relative to the right eye
% param_righteye_pos = param_righteye_dist.*[cos(right_eye_to_param_angle) sin(right_eye_to_param_angle)];
param_righteye_pos = [param_pos_all(:,2) param_pos_all(:,1)];
param_righteye_vel = diff(param_righteye_pos,1,1);
param_righteye_speed = [0;sqrt(param_righteye_vel(:,1).^2 + param_righteye_vel(:,2).^2)];
param_righteye_move_angle = [0;rad2deg(atan(param_righteye_vel(:,2)./param_righteye_vel(:,1)))];
param_righteye_speed = param_righteye_speed(align_with_fluo_low==1);
param_righteye_move_angle = param_righteye_move_angle(align_with_fluo_low==1);
param_righteye_speed_fluo = param_righteye_speed(mask);
param_righteye_move_angle_fluo = param_righteye_move_angle(mask);
param_re_speedx_right = param_righteye_speed_fluo(hunting_idx_wholecourse(anglex<0));
param_re_speedy_right = param_righteye_speed_fluo(simul_bout_idx_wholecourse(angley<0));
param_re_speedz_right = param_righteye_speed_fluo(rest_idx_wholecourse(anglez<0));
param_re_speed3_right = cat(1,param_re_speedx_right,param_re_speedy_right,param_re_speedz_right);
param_re_moveAnglex_right = param_righteye_move_angle_fluo(hunting_idx_wholecourse(anglex<0));
param_re_moveAngley_right = param_righteye_move_angle_fluo(simul_bout_idx_wholecourse(angley<0));
param_re_moveAnglez_right = param_righteye_move_angle_fluo(rest_idx_wholecourse(anglez<0));
param_re_moveAngle3_right = cat(1,param_re_moveAnglex_right,param_re_moveAngley_right,param_re_moveAnglez_right);
[p,tb1] = anovan(spk3_right,{param_re_moveAngle3_right,param_re_speed3_right},'continuous',[1 2],'model','interaction','varname',{'velocity angle','param_speed'});
%%
%decode left/right,distance,righteyeangle using time points at different behavior states
% spk3 = {spkx;spky;spkz};
spk3 = {spkx_right;spky_right;spkz_right};
d3 = {dx_right;dy_right;dz_right};
righteyeAngle3 = {rightanglex_right;rightangley_right;rightanglez_right};
% angle3 = {anglex;angley;anglez};
angle3 = {anglex(anglex<0);angley(angley<0);anglez(anglez<0)};
eval_d3 = zeros(3,1);eval_righteyeAngle = zeros(3,1);eval_angle = zeros(3,1);
for i=1:3
   [f,gof] = fit(spk3{i},d3{i},'poly1');
   eval_d3(i) = gof.rsquare;
   [f,gof] = fit(spk3{i},cos(deg2rad(righteyeAngle3{i})),'poly1');
   eval_righteyeAngle(i) = gof.rsquare;
   [f,gof] = fit(spk3{i},cos(deg2rad(angle3{i})),'poly1');
   eval_angle(i) = gof.rsquare;
end
%%
% plot the evolution trajectory along time

proj_X = W(:,[2 4])'*Spike_X_EstTrace;%2 x time
basecolor = rgb({'magenta','grassgreen'});
color = basecolor(moreIdx+2*lessIdx,:);
% plot_evolution_trajectory(proj_X','marker',reshape((fluo_idx+[-5:5])',[],1),...
%     'color',color(reshape(repmat([1:length(fluo_idx)],11,1),[],1),:));
plot_evolution_trajectory(proj_X','marker',[hunting_idx_start hunting_idx_end],'color',color);
xlabel('PC2');ylabel('PC4');

%plot all bouts' projection
load('tail_swing.mat', 'bout_idx');
[~,hunting_idx_start] = min(abs(row2col(bout_idx(:,1),1) - col2row(mask)),[],2);
[~,hunting_idx_end] = min(abs(row2col(bout_idx(:,2),1) - col2row(mask)),[],2);
moreIdx = param_head_angle_fluo(hunting_idx_start)>0;
lessIdx = param_head_angle_fluo(hunting_idx_start)<0;
basecolor = rgb({'magenta','greenyellow'});
color = basecolor(moreIdx+2*lessIdx,:);
plot_evolution_trajectory(proj_X','marker',[hunting_idx_start hunting_idx_end],'color',color);

pcoi = [1 2 3 4 5];
arrayfun(@(i) generate_MIP(W(:,i),csessionID,cfishID,['_' num2str(i) 'ndPC_random']),pcoi,'un',0);

spk = permute(spk,[2 1 3]);%area x trial x time
% dpca(spk(:,param_head_angle>0,:),);%

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

% Next, we want to find the visual information encoded in neural activity.
param_head_angle_fluo = accumarray(reshape(repmat(1:length(param_head_angle_all)/5,5,1),[],1),param_head_angle_all);
leftIdx = param_head_angle_fluo>0;
rightIdx = param_head_angle_fluo<0;
if nnz(leftIdx)>=nnz(rightIdx)
    moreIdx = leftIdx;lessIdx = rightIdx;
    disp('left right');
    direction = {'left','right'};
    
else
    moreIdx = rightIdx;lessIdx = leftIdx;
    disp('right left');
    direction = {'right','left'};
end
numTrial = nnz(moreIdx);
Spike_X_EstTrace(:,moreIdx);
missingData = repmat(mean(Spike_X_EstTrace(:,lessIdx),2),1,nnz(moreIdx)-nnz(lessIdx));
X = cat(3,Spike_X_EstTrace(:,moreIdx),[Spike_X_EstTrace(:,lessIdx),missingData]);
%region x trial x 2(more&less)
X = permute(X,[1 3 2]);%region x 2 x trial
X_ave = mean(X,3);

%first use averaged data for decomposition
[W,sigma,V] = svd(X_ave,'econ');
sigma = diag(sigma);
proj_X = pagemtimes(W',X);%2 x 2 x trial
meanx = mean(proj_X,3);%2 x 2
stdx = std(proj_X,[],3)./sqrt([nnz(moreIdx) nnz(lessIdx)]);%2 x 2, PC x class
figure,
h = barwitherr(stdx,meanx);
p = arrayfun(@(i) ranksum(squeeze(proj_X(i,1,:)),squeeze(proj_X(i,2,1:nnz(lessIdx)))),1:2);
set(gca,'XTickLabel',{'PC1','PC2'});
ylabel('projected amplitude of all timepoints')
legend(direction);
figure,
bar(sigma/sum(sigma));
ylabel('explained variance');

generate_MIP(W(:,1),csessionID,cfishID,'_visual_1stPC');
generate_MIP(W(:,2),csessionID,cfishID,'_visual_2ndPC');

[W,sigma,V] = svd(X(:,:),'econ');
sigma = diag(sigma);
W = W(:,1:15);
proj_X = pagemtimes(W',X);%pc x 2 x trial
meanx = mean(proj_X,3);%pc x 2
stdx = std(proj_X,[],3)./sqrt([nnz(moreIdx) nnz(lessIdx)]);%2 x 2, PC x class
figure,
h = barwitherr(stdx,meanx);
p = arrayfun(@(i) ranksum(squeeze(proj_X(i,1,:)),squeeze(proj_X(i,2,1:nnz(lessIdx)))),1:2);
set(gca,'XTickLabel',{'PC1','PC2'});
ylabel('projected amplitude of all timepoints')
legend(direction);
figure,
plot(sigma/sum(sigma));
ylabel('explained variance');

bout_idx_aligned_wholecourse = [2 3 4 6 7 10];
arrayfun(@(i) generate_MIP(W(:,i),csessionID,cfishID,['_visual_' num2str(i) 'thPC']),bout_idx_aligned_wholecourse);

%%
%extract the pc of hunting behavior, and project neural activity during
%other behavior states and see which dimension could explains spontaneous
%activity most while not explaining hunting activity, such that it may best
%describe the subspace of hunting behavior
[U,S,V] = svd(mean(spike_focus(:,1:10),3)','econ');
spike_spon = Spike_X_EstTrace(:,setdiff(1:size(Spike_X_EstTrace,2),hunting_idx_visuo));
proj_spon = spike_spon'*V;
proj_focus = spike_focus(:,:)'*V;%time x component
figure,
plot(std(proj_spon,[],1));hold on;
plot(std(proj_focus,[],1));
plot(diag(S));axis tight;
legend({'std(proj spon)','std(proj focus)','svd sigma'});
%%
%compare the similarity between the hunting and spontaneous subspace
S = arrayfun(@(i) spike_spon(:,i+[-5:4]),randperm(size(spike_spon,2)-15,18),'un',0);%{region x time}
spike_spon = cat(3,S{:});
[Us,Ss,Vs] = svd(mean(spike_spon,3)','econ');

%%
%angle encoding
visuo_idx = find(removeStrangeVisuoData(no_param_fluo,param_head_dist_fluo,param_speed_fluo));
%sptial tuning in different state
hunting_idx_visuo = ismember(hunting_idx_wholecourse,visuo_idx);
hunting_idx_visuo = hunting_idx_wholecourse(hunting_idx_visuo);
nonhunting_idx = setdiff(1:length(visuo_idx),hunting_idx_wholecourse);
nonhunting_idx_visuo = ismember(nonhunting_idx,visuo_idx);
nonhunting_idx_visuo = nonhunting_idx(nonhunting_idx_visuo);
simul_idx_visuo = ismember(simul_bout_idx_wholecourse,visuo_idx);
simul_idx_visuo = simul_bout_idx_wholecourse(simul_idx_visuo);
rest_idx_visuo = ismember(rest_idx_wholecourse,visuo_idx);
rest_idx_visuo = rest_idx_wholecourse(rest_idx_visuo);
idx = hunting_idx_visuo;
% idx = nonhunting_idx_visuo;
% idx = simul_idx_visuo;
% idx = rest_idx_visuo;
param_head_angle_visuo = param_head_angle_fluo(idx);
spikes_visuo = Spike_X_EstTrace(:,idx);
tuned2angle(param_head_angle_visuo,spikes_visuo,center,'single');
%%
%investigate 20 degree regional activity during hunting
focus_idx = hunting_idx_visuo(degree20_visuoidx);
spike_focus = zeros(length(region_20degree),15,length(focus_idx));
for i=1:length(focus_idx)
    spike_focus(:,:,i) = Spike_X_EstTrace(region_20degree,focus_idx(i)+[-5:9]);
end
idx = 4;
figure('Position',[2044 169 1278 420]),ax(1) = subplot(1,3,1,polaraxes),ax(2) = subplot(1,3,2),ax(3) = subplot(1,3,3),
cmap = [min(spike_focus(:,:,idx),[],'all') max(spike_focus(:,:,idx),[],'all')];
for i=1:15
%     ax = plot_brain_contour(center);
%     set(gca,ax(1));
    sgtitle(num2str(i)),
    subplot(1,3,1),cla,
    polarplot(deg2rad(param_head_angle_fluo(focus_idx(idx)+[-5:5])),ones(11,1),'-k*');
    polarplot(angledata_fluo(focus_idx(idx)+i-6,:)+pi/2,linspace(0,1,21));
    hold on;
    polarscatter(deg2rad(param_head_angle_fluo(focus_idx(idx)+i-6)),1,20,'r','o','filled');
    polarscatter([deg2rad(lefteye_angle_fluo(focus_idx(idx)+i-6)) deg2rad(righteye_angle_fluo(focus_idx(idx)+i-6))],[.5 .5],20,rgb({'green','blue'}),'o','filled');
    polarscatter([deg2rad(left_eye_angle_fluo(focus_idx(idx)+i-6)) deg2rad(right_eye_angle_fluo(focus_idx(idx)+i-6))],[.25 .25],20,rgb({'green','blue'}),'o','filled');
    title(num2str(conv_or_not_fluo(focus_idx(idx)+i-6)));
    
    subplot(1,3,2),
    scatter(center(region_20degree,1),center(region_20degree,2),20,spike_focus(:,i,idx),'filled');
    caxis(cmap);grid on;

    subplot(1,3,3),
    scatter(center(region_20degree,2),center(region_20degree,3),20,spike_focus(:,i,idx),'filled');
    caxis(cmap);grid on;
%     multiview;
    pause(.5);
end
%%
%investigate 20 degree regional activity during hunting
disp('running this section....');
duration = 15;
focus_idx = hunting_idx_visuo(60);
% focus_idx = simul_idx_visuo(267);
% focus_idx = rest_idx_visuo(140)
spike_focus = zeros(length(region_20degree),duration,length(focus_idx));
for i=1:length(focus_idx)
    spike_focus(:,:,i) = Spike_X_EstTrace(region_20degree,focus_idx(i)+[-5:duration-6]);
end
idx = 1;
writer = VideoWriter(fullfile(getpath('result'),csessionID,['functional_mapping_frame_' num2str(focus_idx(idx))]));
open(writer);
fig = figure('Position',[2044 169 606 673]);
cmap = [min(spike_focus(:,:,idx),[],'all') max(spike_focus(:,:,idx),[],'all')];
for i=1:duration
%     ax = plot_brain_contour(center);
%     set(gca,ax(1));
    clf;
    sgtitle(['frame: ' num2str(i)]),
    ax2 = axes('Position',[[0.2000 0.3000 0.6000 0.5000]]);
    scatter(center(region_20degree,1),center(region_20degree,2),20,spike_focus(:,i,idx),'filled');
    caxis(cmap);grid on;axis tight;
    ax1 = polaraxes('Position',[0.1 0.1 0.8 0.7],'Color','none');
    %prey relative to head
    polarplot(pi/2-deg2rad(param_head_angle_fluo(focus_idx(idx)+[-5:5])),4*ones(11,1),'-k*');
    hold on;
    polarscatter(pi/2-deg2rad(param_head_angle_fluo(focus_idx(idx)+i-6)),4,20,'r','o','filled');
    %tail
    polarplot(pi+angledata_fluo(focus_idx(idx)+i-6,:),1+linspace(0,1,21),'LineWidth',2);
    %two eyes
    polarscatter(pi/2-[deg2rad(lefteye_angle_fluo(focus_idx(idx)+i-6)) deg2rad(righteye_angle_fluo(focus_idx(idx)+i-6))],[3 3],20,rgb({'green','blue'}),'o','filled');
    polarscatter(pi/2-[deg2rad(left_eye_angle_fluo(focus_idx(idx)+i-6)) deg2rad(right_eye_angle_fluo(focus_idx(idx)+i-6))],[2.5 2.5],20,rgb({'green','blue'}),'o','filled');
    title(['conv or not:' num2str(conv_or_not_fluo(focus_idx(idx)+i-6))]);
    set(gca,'Color','none');
%     subplot(1,3,3),
%     scatter(center(region_20degree,2),center(region_20degree,3),20,spike_focus(:,i,idx),'filled');
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
focus_idx = hunting_idx_visuo(ismember(hunting_idx_visuo,hunting_idx_start));
left_param_mask = arrayfun(@(i) mean(param_head_angle_fluo(focus_idx(i)+[-5:1])),1:length(focus_idx))>0;
focus_idx = focus_idx(~left_param_mask);
spike_focus = zeros(length(region_20degree),duration,length(focus_idx));
for i=1:length(focus_idx)
    spike_focus(:,:,i) = Spike_X_EstTrace(region_20degree,focus_idx(i)+[-5:duration-6]);
end
spike_focus = Spike_X_EstTrace(region_20degree,750:784);
spike_focus = mean(spike_focus,3);
idx = 1;
writer = VideoWriter(fullfile(getpath('result'),csessionID,['functional_mapping_frame_' num2str(focus_idx(idx))]));
open(writer);
fig = figure('Position',[2044 169 606 673]);
cmap = [min(spike_focus(:,:,idx),[],'all') max(spike_focus(:,:,idx),[],'all')];
for i=1:duration
%     ax = plot_brain_contour(center);
%     set(gca,ax(1));
    clf;
    sgtitle(['frame: ' num2str(i)]),
    ax2 = axes('Position',[[0.2000 0.3000 0.6000 0.5000]]);
    scatter(center(region_20degree,1),center(region_20degree,2),20,spike_focus(:,i,idx),'filled');
    caxis(cmap);grid on;axis tight;
    ax1 = polaraxes('Position',[0.1 0.1 0.8 0.7],'Color','none');
    %prey relative to head
    polarplot(pi/2-deg2rad(param_head_angle_fluo(focus_idx(idx)+[-5:5])),4*ones(11,1),'-k*');
    hold on;
    polarscatter(pi/2-deg2rad(param_head_angle_fluo(focus_idx(idx)+i-6)),4,20,'r','o','filled');
    %tail
    polarplot(pi+angledata_fluo(focus_idx(idx)+i-6,:),1+linspace(0,1,21),'LineWidth',2);
    %two eyes
    polarscatter(pi/2-[deg2rad(lefteye_angle_fluo(focus_idx(idx)+i-6)) deg2rad(righteye_angle_fluo(focus_idx(idx)+i-6))],[3 3],20,rgb({'green','blue'}),'o','filled');
    polarscatter(pi/2-[deg2rad(left_eye_angle_fluo(focus_idx(idx)+i-6)) deg2rad(right_eye_angle_fluo(focus_idx(idx)+i-6))],[2.5 2.5],20,rgb({'green','blue'}),'o','filled');
    title(['conv or not:' num2str(conv_or_not_fluo(focus_idx(idx)+i-6))]);
    set(gca,'Color','none');
%     subplot(1,3,3),
%     scatter(center(region_20degree,2),center(region_20degree,3),20,spike_focus(:,i,idx),'filled');
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
focus_idx = hunting_idx_visuo(ismember(hunting_idx_visuo,hunting_idx_start));
left_param_mask = arrayfun(@(i) mean(param_head_angle_fluo(focus_idx(i)+[-5:1])),1:length(focus_idx))>0;
focus_idx = focus_idx(~left_param_mask);
spike_focus = zeros(length(region_20degree),duration,length(focus_idx));
for i=1:length(focus_idx)
    spike_focus(:,:,i) = Spike_X_EstTrace(region_20degree,focus_idx(i)+[-5:duration-6]);
end%ignore the actual difference of visual stimulus, like different angle, distance, direction, speed, etc.
%clustering the areas with very similar response 
%cross-correlation between regions to find out the info flow route
del_region = sum(all(spike_focus==0,2),3)>0.1*length(focus_idx);
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
focus_idx = randperm(size(Spike_X_EstTrace,2)-duration,length(focus_idx));
spike_focus = zeros(length(region_20degree),duration,length(focus_idx));
for i=1:length(focus_idx)
    spike_focus(:,:,i) = Spike_X_EstTrace(region_20degree,focus_idx(i)+[-5:duration-6]);
end%ignore the actual difference of visual stimulus, like different angle, distance, direction, speed, etc.
%%
spike_focus = mean(spike_focus,3);
thld = quantile(spike_focus(:),.6);
detect_region = median(spike_focus(:,4:13),2)>thld;
[activation,act_time] = max(spike_focus(detect_region,4:13),[],2);
plot_brain_contour(center,'3d');
scatter3(center(keep_region(detect_region),1),center(keep_region(detect_region),2),center(keep_region(detect_region),3),10,act_time,'filled');
colormap('jet');