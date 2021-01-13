function [areaName,perc] = lookuparea(queryRegion,sessionID)
load(fullfile(getpath('neural activity','image segmentation'),'ZBB anatomy','pixel2area'));
load(fullfile(getpath('neural activity',sessionID),'Coherence3'),'A3');
pixelIdx = find(A3(:,queryRegion));
areaName = arrayfun(@(i) tbl.(1)((tbl.(2)==i)),pixelIdx,'un',0);
areaName = cat(1,areaName{:});
[areaName,ia,ic] = unique(areaName);
A = accumarray(ic,1);
perc = A/sum(A);
% [idx,I] = max(A);
% area = areaName(I);
% perc = idx/sum(A);
end

