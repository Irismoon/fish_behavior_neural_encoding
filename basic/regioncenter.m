function regioncenter(csessionID,cfishID)
%function regioncenter(csessionID,cfishID)
%find the center pixel location of each region
folderinfo = dir(fullfile(getpath('neural activity',csessionID,cfishID),'Coherence3*'));
datenum = arrayfun(@(i) folderinfo(i).datenum,1:length(folderinfo));
[~,I] = max(datenum);
filepath = fullfile(folderinfo(I).folder,folderinfo(I).name);
load(filepath,'A3');
[nPixel,nRegion] = size(A3);
I = arrayfun(@(i) find(A3(:,i)),1:nRegion,'un',0);
center = zeros(nRegion,3);
if ismember(csessionID,{'201117'})
    sz = [376 308 210];
else 
    sz = [600 600 280];
end
for i=1:nRegion
    [row,col,page] = ind2sub([sz(1) sz(2) sz(3)],I{i});
    center(i,:) = [mean(row) mean(col) mean(page)];
end
%calculate the distance between center
regiondist = squareform(pdist(center));
save(filepath,'regiondist','center','-append');
end

