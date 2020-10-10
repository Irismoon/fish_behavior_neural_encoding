function bout_idx = detectBout(sessionID,fishID)
load(fullfile(getpath('behavior',sessionID,fishID),'tail_swing'),'left_tail_swing','right_tail_swing','sum_curv');
tail_swing = left_tail_swing+abs(right_tail_swing);
%detect continuous tail swing start and end
diff_swing = diff(tail_swing);
startIdx = find(diff_swing==1);
endIdx = find(diff_swing==-1);
%if the previous end is very close to the next start, combine them into one
bout_idx = [startIdx endIdx];
bout_idx_toy = reshape(bout_idx',[],1);
idxtmp = (diff(bout_idx_toy));
mask = find(idxtmp(2:2:end)<4);
while ~isempty(mask)
    bout_idx_new = [];
    for i=1:length(bout_idx)
        if ismember(i,mask)
            bout_idx_new = [bout_idx_new;bout_idx(i,1) bout_idx(i+1,2)];
        else
            bout_idx_new = [bout_idx_new;bout_idx(i,:)];
        end
    end
    delidx = find(diff(bout_idx_new(:,2))==0)+1;
    bout_idx_new(delidx,:) = [];
    bout_idx_toy = reshape(bout_idx_new',[],1);
    bout_idx = bout_idx_new;
    idxtmp = diff(bout_idx_toy);
    mask = find(idxtmp(2:2:end)<4);
end
save(fullfile(getpath('behavior',sessionID,fishID),'tail_swing'),'bout_idx','-append');
figure,
plot(sum_curv(bout_idx(2,1):bout_idx(4,2)));
hold on;
ax = gca;
mark_line(ax,[bout_idx(2,2) bout_idx(3,1) bout_idx(3,2) bout_idx(4,1)]-bout_idx(2,1)+1);
title('three example bouts');
end

