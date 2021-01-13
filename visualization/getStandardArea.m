function getStandardArea()
folderpath = fullfile(getpath('neural activity','image segmentation'),'ZBB anatomy','*.tif');
file = dir(folderpath);
pixelofarea = cell(length(file),1);
for idx = 1:length(file)
    imagpath = fullfile(file(idx).folder,file(idx).name);
    imag = arrayfun(@(iframe) imread(imagpath,iframe),1:length(imfinfo(imagpath)),'un',0);
    imag = cellfun(@(x) x(:,1:376)',imag,'un',0);
    imag = cat(3,imag{:});
    region = find(imag);
    pixelofarea{idx} = region;
end
pixelNo = cat(1,pixelofarea{:});
areaName = row2col({file(:).name},1);
areaName = cellfun(@(x) strsplit(x,'.'),areaName,'un',0);
areaName = cellfun(@(x) x{1},areaName,'un',0);
label = arrayfun(@(x) x*ones(length(pixelofarea{x}),1),1:length(file),'un',0);
label = cat(1,label{:});
areaName = areaName(label);
tbl = table(areaName,pixelNo);
save(fullfile(getpath('neural activity','image segmentation'),'ZBB anatomy','pixel2area'),'tbl');
end

