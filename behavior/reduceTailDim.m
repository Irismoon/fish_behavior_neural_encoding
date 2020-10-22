function reduceTailDim(postfix)
%1.extract tail swing segment
[sessionID,fishID] = getfish();
numFish = length(sessionID);
bout_all_curv = cell(numFish,1);
bout_all_xy = bout_all_curv;
single_bout_len = cell(numFish,1);
for iFish=1:numFish
    load(fullfile(getpath('behavior',sessionID{iFish},fishID{iFish}),'tail_swing'),'curvdata','bout_idx','new_centerline_20');
    single_bout_len{iFish} = bout_idx(:,2) - bout_idx(:,1)+1;
    all_bout_len = sum(single_bout_len{iFish});
    single_fish_bout_curv = zeros(all_bout_len,20);
    single_fish_bout_xy = zeros(all_bout_len,20,2);
    for iBout = 1:length(bout_idx)
        curv = curvdata(bout_idx(iBout,1):bout_idx(iBout,2),:);
        single_fish_bout_curv(sum(single_bout_len{iFish}(1:iBout-1))+1:sum(single_bout_len{iFish}(1:iBout)),:) = curv;
        single_fish_bout_xy(sum(single_bout_len{iFish}(1:iBout-1))+1:sum(single_bout_len{iFish}(1:iBout)),:,:) = new_centerline_20(bout_idx(iBout,1):bout_idx(iBout,2),:,:);
    end
    bout_all_curv{iFish} = single_fish_bout_curv;
    bout_all_xy{iFish} = single_fish_bout_xy;
end
bout_all_curv = cat(1,bout_all_curv{:});%
bout_all_xy = cat(1,bout_all_xy{:});
single_bout_len = cat(1,single_bout_len{:});
bout_single_curv = arrayfun(@(i) bout_all_curv(sum(single_bout_len(1:i-1))+(1:single_bout_len(i)),:),1:length(single_bout_len),'un',0);
bout_single_xy = arrayfun(@(i) bout_all_xy(sum(single_bout_len(1:i-1))+(1:single_bout_len(i)),:,:),1:length(single_bout_len),'un',0);
assert(length(bout_all_curv)==sum(single_bout_len),'the length must be the same!');
%2.pca
[coeff,score,latent,~,explained] = pca(bout_all_curv);
score_single = arrayfun(@(i) score(sum(single_bout_len(1:i-1))+(1:single_bout_len(i)),:),1:length(single_bout_len),'un',0);
S = matdoc('comment','use eight fish');
save(fullfile(getpath('result'),['tail_pca' postfix]),'coeff','score','score_single','latent','explained','bout_all_curv','bout_single_curv',...
    'single_bout_len','bout_single_xy');

figure,
plot(cumsum(explained));
xlabel('# PC');ylabel('explained variance');
title('reduce 20 dimensions of tail angle');
axis tight;
figure,
plot(coeff(:,1:4));
legend({'PC 1','PC 2','PC 3','PC 4'});
axis tight;
ylabel('angle');
%%
%how does one bout project on the PC space
figure,
for i=1:4
    subplot(2,2,i)
    idx = randperm(length(single_bout_len),1);
    single_bout = score(sum(single_bout_len(1:idx-1))+1:sum(single_bout_len(1:idx)),1:3);
    single_bout(end,2:3) = nan;
    single_bout(end+1:end+2,:) = nan;
    patch(single_bout(:,1),single_bout(:,2),single_bout(:,3),1:single_bout_len(idx)+2,'EdgeColor','interp');
    view(45,45);colorbar;
end
end

