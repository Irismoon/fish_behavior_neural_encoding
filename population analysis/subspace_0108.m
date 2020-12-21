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
% [W,V,whichMarg] =   .techniques.dpca(mean(spk_balanced,4),5,'combinedParams',{{1},{2},{1,2}});
spk_balanced_ave = mean(spk_balanced,4);%neuron x direction x time
X = spk_balanced_ave(:,:);%region x (direction x time)
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
% plot the evolution trajectory along time

proj_X = W(:,[2 4])'*Spike_X_EstTrace;%2 x time
basecolor = rgb({'magenta','greenyellow'});
color = basecolor(moreIdx+2*lessIdx,:);
% plot_evolution_trajectory(proj_X','marker',reshape((fluo_idx+[-5:5])',[],1),...
%     'color',color(reshape(repmat([1:length(fluo_idx)],11,1),[],1),:));
plot_evolution_trajectory(proj_X','marker',[hunting_idx_start hunting_idx_end],'color',color);
xlabel('PC2');ylabel('PC4');

%plot all bouts' projection
moreIdx = param_head_angle_fluo(bout_idx_aligned_start)>0;
lessIdx = param_head_angle_fluo(bout_idx_aligned_start)<0;
basecolor = rgb({'magenta','greenyellow'});
color = basecolor(moreIdx+2*lessIdx,:);
plot_evolution_trajectory(proj_X','marker',[bout_idx_aligned_start bout_idx_aligned_end],'color',color);
%%
%extract the pc of hunting behavior, and project neural activity during
%other behavior states and see which dimension could explains spontaneous
%activity most while not explaining hunting activity, such that it may best
%describe the subspace of hunting behavior
focus_idx = arrayfun(@(i) hunting_idx_start(i)+[-5:4]',1:length(hunting_idx_start),'un',0);
focus_idx = cat(1,focus_idx{:});
% focus_idx = focus_idx(ismember(focus_idx,visuo_idx));
spike_focus = reshape(Spike_X_EstTrace(:,focus_idx),[],10,length(hunting_idx_start));
param_focus = param_head_angle_fluo(focus_idx);
spon_idx = setdiff((1:size(Spike_X_EstTrace,2))',focus_idx);
spon_idx = spon_idx(ismember(spon_idx,visuo_idx));
spike_spon = Spike_X_EstTrace(:,spon_idx);
param_spon = param_head_angle_fluo(spon_idx);
[U,S,V] = svd(mean(spike_focus,3)','econ');
S = diag(S.^2)/sum(diag(S.^2));
proj_spon = spike_spon'*V;
proj_focus = spike_focus(:,:)'*V;%time x component
figure,
plot(sum(proj_spon.^2)/sum((spike_spon-mean(spike_spon,1)).^2,'all'));hold on;
plot(sum(proj_focus.^2)/sum((spike_focus-mean(spike_focus,1)).^2,'all'));
plot(diag(S.^2)/sum(diag(S.^2)));axis tight;%
legend({'proj spon','proj focus','svd sigma^2'});
ylabel('variance(projected response)/variance(original response)');
xlabel('# projected dimension');
%%
%compare the similarity between the hunting and spontaneous subspace
spike_spon_trial = arrayfun(@(i) spike_spon(:,i+[-5:4]),randperm(size(spike_spon,2)-15,size(spike_focus,3)),'un',0);%{region x time}
spike_spon_trial = cat(3,spike_spon_trial{:});
[Us,Ss,Vs] = svd(mean(spike_spon_trial,3)','econ');
% [Us,Ss,Vs] = svd(spike_spon_all','econ');
Ss = diag(Ss.^2)/sum(diag(Ss.^2));
princIdxs = find(cumsum(Ss)>0.9,1);
princIdx = find(cumsum(S)>0.9,1);
figure,
coef = V'*Vs;
imagesc(coef);colorbar;caxis([-1,1]);
xlabel('hunting');ylabel('spon');
title('correlation of hunting subspace and spontaneous subspace');
%%
%use the common mode to regress visual, to see if the common basis
%represents the visual stimuli in both cases
%the 1st and 2nd PC
figure,
subplot(2,1,1),
proj_focus = pagemtimes(V',spike_focus);%dim x time x trial
proj_focus = proj_focus(:,:);%dim x sample
mdl_focus = fitlm(proj_focus',param_focus);
stem(log10(mdl_focus.Coefficients.pValue(2:end)));axis tight;hold on;
proj_focus = pagemtimes(Vs',spike_focus);%dim x time x trial
proj_focus = proj_focus(:,:);%dim x sample
mdl_focus = fitlm(proj_focus',param_focus);
stem(log10(mdl_focus.Coefficients.pValue(2:end)));axis tight;
ylabel('log10(pvalue)');xlabel('dimension');
legend({'hunting pc','spon pc'});
title('hunting neural activity');
 
subplot(2,1,2),
proj_spon = pagemtimes(V',spike_spon);%dim x sample
mdl_spon = fitlm(proj_spon',param_spon);
stem(log10(mdl_spon.Coefficients.pValue(2:end)));axis tight;
hold on;
proj_spon = pagemtimes(Vs',spike_spon);%dim x sample
mdl_spon = fitlm(proj_spon',param_spon);
stem(log10(mdl_spon.Coefficients.pValue(2:end)));axis tight;
title('spontaneous neural activity');
%%
%calculate the alignment index
alignmentIdx = diag(V'*mean(spike_focus,3)*mean(spike_focus,3)'*V)/sum((S.^2))
%small alignment index indicates high orthogonal between hunting and spon
%%
%measure if the hunting activity has large variance in the found space
%while spontanous activity has small variance
DataStruct(1).A = Spike_X_EstTrace(:,hunting_idx_wholecourse)';
DataStruct(2).A = Spike_X_EstTrace(:,nonhunting_idx)';
% DataStruct(1).A = mean(spike_focus,3)';
% DataStruct(2).A = mean(spike_spon_trial,3)';
k = 10;
DataStruct(1).dim = k;
DataStruct(2).dim = k;
[QSubspaces] = maxvar_subspaces(DataStruct);
figure,
C = QSubspaces(1).Q'*QSubspaces(2).Q;
imagesc(C);caxis([-1,1]);
title(['costFn: 1:' num2str(QSubspaces(1).costFn) ' 2:' num2str(QSubspaces(2).costFn)]);
xlabel('spontaneous');ylabel('hunting');
%%
figure,
multchnPlot(reshape(proj_spon,2242,1,10),[1 2242],conv_hunting_start,'plotstyle','plot');
%%
%project original signal into the orthogonal subspaces such that projected
%signal is orthogonal and repr,esents the part uniquely belonging to each
%state
proj_focus = DataStruct(1).A*QSubspaces(1).Q;
proj_spon = DataStruct(2).A*QSubspaces(2).Q;%sample x 10
generate_MIP(abs(QSubspaces(1).Q(:,1)),csessionID,cfishID,'_1stOrtho_hunting');
generate_MIP(abs(QSubspaces(1).Q(:,2)),csessionID,cfishID,'_2ndOrtho_hunting');
%%
align_idx = zeros(k,2,2);%dim x Q(hunting/spon) x data(hunting/spon)
figure,
state = {'hunting','spon'};
for i=1:2
    align_idx(:,1,i) = align_index(DataStruct(i).A,QSubspaces(1).Q);
    align_idx(:,2,i) = align_index(DataStruct(i).A,QSubspaces(2).Q);
    subplot(2,1,i),
    imagesc(align_idx(:,:,i)');colormap;cbar = colorbar;%caxis([-1,1]);
    xlabel('dim');ax = gca;ylabel('Q');
    title(cbar,'explained var');
    ax.YTick = [1 2];ax.YTickLabel = state;
    title(state{i});
end
%%
%hunting-specific signal
moreIdx = param_head_angle_fluo(hunting_idx_start)>0;
lessIdx = param_head_angle_fluo(hunting_idx_start)<0;
basecolor = rgb({'magenta','greenyellow'});
color = basecolor(moreIdx+2*lessIdx,:);
plot_evolution_trajectory(Spike_X_EstTrace'*QSubspaces(1).Q(:,1:2),'marker',[hunting_idx_start hunting_idx_end],'color',color,'plottype','plot');
%%
%split the data into three parts: spontaneous, hunting, transition between
%them
conv_dif = diff(conv_or_not_fluo);
idx = reshape(find(conv_dif),2,[])';
conv_hunting_start = idx(:,1);conv_hunting_end = idx(:,2);
conv_hunting_wholecourse = [];
for i=1:length(conv_hunting_start)
    conv_hunting_wholecourse = [conv_hunting_wholecourse;[conv_hunting_start(i):conv_hunting_end(i)]'];
end
transition_idx = (conv_hunting_start+[-10:-1])';%time x trial
transition_idx = transition_idx(:);
spon_idx = setdiff(1:numWholeTime,[transition_idx;conv_hunting_wholecourse]);
DataStruct(1).A = Spike_X_EstTrace(:,conv_idx)';
DataStruct(2).A = Spike_X_EstTrace(:,transition_idx)';
k = 10;
DataStruct(1).dim = k;
DataStruct(2).dim = k;
[QSubspaces] = maxvar_subspaces(DataStruct);
figure,
C = QSubspaces(1).Q'*QSubspaces(2).Q;
imagesc(C);caxis([-1,1]);
title(['costFn: 1:' num2str(QSubspaces(1).costFn) ' 2:' num2str(QSubspaces(2).costFn)]);
xlabel('spontaneous');ylabel('hunting');
Q_conv = QSubspaces(1).Q;
Q_trans = QSubspaces(2).Q;
%
DataStruct(1).A = Spike_X_EstTrace(:,spon_idx)';
DataStruct(2).A = Spike_X_EstTrace(:,transition_idx)';
k = 10;
DataStruct(1).dim = k;
DataStruct(2).dim = k;
[QSubspaces] = maxvar_subspaces(DataStruct);
Q_spon = QSubspaces(1).Q;
Q_trans2 = QSubspaces(2).Q;

figure,
C = QSubspaces(1).Q'*QSubspaces(2).Q;
imagesc(C);caxis([-1,1]);
title(['costFn: 1:' num2str(QSubspaces(1).costFn) ' 2:' num2str(QSubspaces(2).costFn)]);
xlabel('spontaneous');ylabel('hunting');