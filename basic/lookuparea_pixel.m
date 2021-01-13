function areaName = lookuparea_pixel(querycoor,sessionID)
load(fullfile(getpath('neural activity','image segmentation'),'ZBB anatomy','pixel2area'));
pixelIdx = querycoor;
areaName = tbl.(1)((tbl.(2)==pixelIdx));
end