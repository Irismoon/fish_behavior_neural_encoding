function boutdist = dtw_alignBout(postfix)
%1.align bouts by dtw technique
%2.get the distance from dtw as similarity metric between bouts
load(fullfile(getpath('result'),['tail_pca' postfix]),'coeff','score_single','latent','explained','bout_single_curv','single_bout_len');
ncomp = find(cumsum(explained)>0.98,1,'first');
numbout = length(single_bout_len);
pair = nchoosek(1:numbout,2);
pair = [pair(:,2) pair(:,1)];
single_bout = score_single;
dist1 = arrayfun(@(i) partial_dtw(single_bout{pair(i,1)}(:,1:ncomp)',single_bout{pair(i,2)}(:,1:ncomp)'),(1:size(pair,1))');
dist2 = arrayfun(@(i) partial_dtw(single_bout{pair(i,1)}(:,1:ncomp)',-single_bout{pair(i,2)}(:,1:ncomp)'),(1:size(pair,1))');
dist = min([dist1 dist2],[],2);
mask = tril(true(numbout,numbout),-1);
boutdist = zeros(numbout,numbout);
boutdist(mask) = dist;
boutdist = boutdist+boutdist';
S = matdoc();
save(fullfile(getpath('result'),['tail_boutdist' postfix]),'boutdist','S');
end
function dist = partial_dtw(x,y)
len1 = size(x,2);
len2 = size(y,2);
if len1>len2 long = x;short=y; len_long = len1; len_short = len2; else long=y;short=x; len_long = len2; len_short = len1; end
if abs(len1-len2)<=5
    dist = dtw(long,short);
else
    dist = dtw(long(:,1:len_short+5),short);%only use part of the long bout
end
end