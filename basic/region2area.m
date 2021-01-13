function region2area(sessionID)
load(fullfile(getpath('neural activity',sessionID),'Coherence3'),'A3');
numRegion = size(A3,2);
[areaName,perc] = arrayfun(@(i) lookuparea(i,sessionID),(1:numRegion)','un',0);
regionNo = (1:numRegion)';
tbl = table(regionNo,areaName,perc);
S = matdoc('comment','if there is more than one area included in one region, its percent of each area is indicated in the perc column');
save(fullfile(getpath('neural activity',sessionID),'region2standardArea'),'tbl','S');
end

