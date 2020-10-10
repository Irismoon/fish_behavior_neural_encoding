function boutdist = dtw_alignBout()
%1.align bouts by dtw technique
%2.get the distance from dtw as similarity metric between bouts
load(fullfile(getpath('result'),'tail_pca'),'coeff','score','latent','bout_all','single_bout_len');
ncomp = 3;
numbout = length(single_bout_len);
single_bout = arrayfun(@(iBout) bout_all(...
    sum(single_bout_len(1:iBout-1))+1:sum(single_bout_len(1:iBout)),1:ncomp)',...
    (1:numbout)','un',0);%comp x time
pair = nchoosek(1:numbout,2);
pair = [pair(:,2) pair(:,1)];
dist1 = arrayfun(@(i) dtw(single_bout{pair(i,1)},single_bout{pair(i,2)}),(1:size(pair,1))');
dist2 = arrayfun(@(i) dtw(single_bout{pair(i,1)},-single_bout{pair(i,2)}),(1:size(pair,1))');
dist = min([dist1 dist2],[],2);
mask = tril(true(numbout,numbout),-1);
boutdist = zeros(numbout,numbout);
boutdist(mask) = dist;
boutdist = boutdist+boutdist';
S = matdoc();
save(fullfile(getpath('result'),'tail_boutdist'),'boutdist','S');
end