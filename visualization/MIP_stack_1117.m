%generate 1117 MIP stack
prepath = fullfile(getpath('imaging','201117','1'),'*MIP*');
fileinfo = dir(prepath);
t = arrayfun(@(i) regexp(fileinfo(i).name,'\d{2}','match'),1:length(fileinfo),'un',0);
t = cellfun(@(x) str2double(x{1}),t);
[t,I] = sort(t,'ascend');
fileinfo = fileinfo(I);

k=0;
stackpath = fullfile(getpath('neural activity',sessionID,fishID),'MIP_affine_stack.tif');
for ifile=1:length(fileinfo)
    filepath = fullfile(fileinfo(ifile).folder,fileinfo(ifile).name,'MIP_affineX*.tif');
    tifinfo = dir(filepath);
    t = arrayfun(@(i) regexp(tifinfo(i).name,'\d*','match'),1:length(tifinfo),'un',0);
    t = cellfun(@(x) str2double(x{1}),t);
    [t,I] = sort(t,'ascend');
    tifinfo = tifinfo(I);
    for itif = 1:length(tifinfo)
        k=k+1;
        frame_tmp = imread(fullfile(tifinfo(itif).folder,tifinfo(itif).name));
        if k==1
        imwrite(frame_tmp,stackpath);
        else
            imwrite(frame_tmp,stackpath,'WriteMode','append');
        end
    end
end