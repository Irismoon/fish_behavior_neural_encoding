function combine_image_single2stack(filepath,sessionID,fishID)
%%
%get all the subfolder
subfolder = dir(filepath);
mask_folder = arrayfun(@(i) contains(subfolder(i).name,'MIP'),1:length(subfolder));
subfolder = subfolder(mask_folder);
name = {subfolder(:).name};
idx = cellfun(@(x) str2num(regexp(x,'\d*','match','once')),name);
[~,I] = sort(idx,'ascend');
subfolder = subfolder(I);
file_out = fullfile(getpath('imaging',sessionID,fishID),'MIP_affine_stack.tif');
%%
%for each subfolder, get all files
%detect stack, if yes, replace it by interpolated frames
%write it into the single stack file
numImage = 0;
for iSubf = 1:length(subfolder)
    file = dir(fullfile(subfolder(iSubf).folder,subfolder(iSubf).name));
    mask_file = arrayfun(@(i) isfile(fullfile(file(i).folder,file(i).name)),1:length(file));
    file = file(mask_file);
    name = {file(:).name};
    mask_interp = cellfun(@(x) contains(x,'interpolation'),name);
    mask_stack = cellfun(@(x) contains(x,'stack'),name);
    whichfileIsInterp = find(mask_interp);
    name_interp = name(whichfileIsInterp);
    whichframeIsInterp = cellfun(@(x) str2num(regexp(x,'\d*','match','once')),name_interp);
    if any(mask_stack)
        idx_stack = find(mask_stack);
        stack_filepath = fullfile(file(idx_stack).folder,file(idx_stack).name);
        imginfo = imfinfo(stack_filepath);
        for iImg = 1:length(imginfo)
            single_stack_image = imread(stack_filepath,iImg);
            idxtmp = find(iImg==whichframeIsInterp);
            if ~isempty(idxtmp)
                interpfile = file(whichfileIsInterp(idxtmp));
                single_stack_image = imread(fullfile(interpfile.folder,interpfile.name));
            end
            numImage = numImage + 1;
            if numImage==1
                imwrite(single_stack_image,file_out);
            elseif numImage>1
                imwrite(single_stack_image,file_out,'WriteMode','append');
            end
        end
    else
    %single images and interp 
       whichfileIsSingle = setdiff(1:length(file),whichfileIsInterp);
       [~,I] = sort(cellfun(@(x) str2num(regexp(x,'\d*','match','once')),name(whichfileIsSingle)),'ascend');
       whichfileIsSingle = whichfileIsSingle(I);
       name_single = name(whichfileIsSingle);
       for iImg = 1:length(whichfileIsSingle)
           single_stack_image = imread(fullfile(file(1).folder,name_single{iImg}));
           idxtmp = find(iImg==whichframeIsInterp);
           if ~isempty(idxtmp)
               interpfile = fullfile(file(1).folder,file(whichfileIsInterp(idxtmp)).name);
               single_stack_image = imread(interpfile);
           end
           numImage = numImage + 1;
           if numImage==1
               imwrite(single_stack_image,file_out);
           elseif numImage>1
               imwrite(single_stack_image,file_out,'WriteMode','append');
           end
       end
    end
end

