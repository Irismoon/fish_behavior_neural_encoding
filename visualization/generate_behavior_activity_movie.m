%% Put the figure of activity or some behavior parameters such as eye angle
%% and the videos (high and low magnification videos and MIP of fluorescence images) together
%% to form a video to monitor the relation between activity an behavior.

%% input: sequences, length_sequences, bad_sequences, FiringRate, A3,
%% path_high, path_low, path_MIP_ori, path_MIP_aligned, path_out,
%% output: a video saved at path_out
function generate_behavior_activity_movie(sessionID,fishID,behavior_var)
% function generate_behavior_activity_movie(sessionID,fishID,behavior_var)
%behavior_var could be a char vector or a cell array containing char arrays
%it could be any behavior variable you are interested in 
if ~iscell(behavior_var)
    behavior_var = {behavior_var};
end
%%
parentpath = getpath('behavior',sessionID,fishID);
load(fullfile(parentpath,'align_with_fluo'));
%%
fileinfo = dir(fullfile(parentpath,'*high_align_with_fluo.avi'));
fileName_high = fullfile(fileinfo.folder,fileinfo.name);
obj_high = VideoReader(fileName_high);

fileinfo = dir(fullfile(parentpath,'*low_align_with_fluo.avi'));
fileName_low = fullfile(fileinfo.folder,fileinfo.name);
obj_low = VideoReader(fileName_low);

nframe = 0;
vidFrame_high = [];
while hasFrame(obj_high)
    nframe = nframe+1;
    vidtmp = readFrame(obj_high);
    if mod(nframe,5)==1
        vidFrame_high = cat(4,vidFrame_high,vidtmp);
    end
end
nframe = 0;
vidFrame_low = [];
while hasFrame(obj_low)
    nframe = nframe+1;
    vidtmp = readFrame(obj_low);
    if mod(nframe,5)==1
        vidFrame_low = cat(4,vidFrame_low,vidtmp);
    end
end

fileName_MIP_DFoF = fullfile(getpath('imaging',sessionID,fishID),'MIP_DFoF_Stack_contrastAdjusted.tif');
fileName_MIP_affine = fullfile(getpath('imaging',sessionID,fishID),'MIP_affine_stack.tif');
%%
idx_keep_frame_high = find(align_with_fluo_high);
idx_keep_frame_high = idx_keep_frame_high(1:5:length(idx_keep_frame_high));
filename = behvar2filename(behavior_var);
fig = figure('Visible','off');
hold on;
for i = 1:length(behavior_var)
    load(fullfile(getpath('behavior',sessionID,fishID),filename{i}),behavior_var{i});
    eval([behavior_var{i} '=' behavior_var{i} '(idx_keep_frame_high);']);
    eval([behavior_var{i} ' = zscore(' behavior_var{i} ');']);
    tagname = behavior_var{i};
    eval(['plot(' behavior_var{i} ',"tag",tagname);']);
end
xlabel('frame(0.1s)');
legends = behavior_var;
set(gcf,'position',[0,0,2440,600]);axis tight;
%%
path_out = getpath('behavior',sessionID,fishID);
fileName_out = fullfile(path_out,['behavior_activity_movie_' strjoin(behavior_var,'&') '.avi']);
obj_out = VideoWriter(fileName_out);
open(obj_out);
%% prepare the input and output files
%% for each frame, read the pictures, plot, and put them together, write into a vedio
imagInfo = imfinfo(fullfile(getpath('imaging',sessionID,fishID),'MIP_DFoF_Stack_contrastAdjusted.tif'));
numt = length(imagInfo);
fprintf('starting writing into movie...');
for it = 1:numt
    if mod(it,500)==0 disp(num2str(it)); end
    frame_high = vidFrame_high(:,:,:,it);%video is five times faster than imaging
    frame_low = vidFrame_low(:,:,:,it);
    frame_high = flip(frame_high,1);
    frame_low = flip(frame_low,1);
    frame_MIP_DFoF = imread(fileName_MIP_DFoF,it);
    frame_MIP_DFoF = repmat(uint8(double(frame_MIP_DFoF)/65535*256),1,1,3);%contrast adjusted
%     frame_MIP_DFoF = uint8(repmat(double(frame_MIP_DFoF),1,1,3)*900/65535); % 3 channels, uint16 to uint8(all data above 255 will saturate, this will increase the brightness and decrease the contrast)
    frame_MIP_affine = imread(fileName_MIP_affine,it);
    frame_MIP_affine = repmat(frame_MIP_affine,1,1,3);
    
    for iBeh = 1:length(behavior_var)
        h = findobj(fig,'type','line','tag',behavior_var{iBeh}); % request the color of the relevant line
        eval(['plot(it,' behavior_var{iBeh} '(it),"o","Color",h.Color,"tag","mark_' behavior_var{iBeh} '");']); % draw a mark for the present frame for each line
        eval(['text(it,' behavior_var{iBeh} '(it),sprintf("(%.3f)",' behavior_var{iBeh} '(it)),"VerticalAlignment","bottom","Color",h.Color);']);
    end
    legend(legends,'Interpreter','none'); % This should be done for every single frame to avoid adding legends for the marks.
    title(['frame ',int2str(it)]);
    frame_plot = getframe(gcf);
    frame_plot = frame_plot.cdata;
    
    %% put them together to form a picture
    width_out = size(frame_low,2) + size(frame_MIP_DFoF,2) + size(frame_MIP_affine,2); % the positions of these panels are a little complex
    if size(frame_plot,2)>width_out
        width_out = size(frame_plot,2);
    end
    height_out = size(frame_MIP_DFoF,1) + size(frame_plot,1);
    frame_out = uint8(zeros(height_out,width_out,3));
    %----------------------------------------%
    %|                 (plot)              |
    %-----------------------------------------
    %|      (low)     |          |          |
    %|----------------| MIP_DFoF |MIP_affine|
    %|     (high)     |          |          |
    [hplot,wplot,~] = size(frame_plot);
    [hlow,wlow,~] = size(frame_low);
    [hhigh,whigh,~] = size(frame_high);
    [hMIP_d,wMIP_d,~] = size(frame_MIP_DFoF);
    [hMIP_a,wMIP_a,~] = size(frame_MIP_affine);
    frame_out(1:hplot,1:wplot,:) = frame_plot; % top: frame_plot
    frame_out(hplot+1:hplot+hlow,1:wlow,:) = frame_low; % bottom left first row: frame_low
    frame_out(hplot+hlow+1:hplot+hlow+hhigh,1:whigh,:) = frame_high; % bottom left second row: frame_high
    frame_out(hplot+1:hplot+hMIP_d,max(wlow,whigh)+1:max(wlow,whigh)+wMIP_d,:) = frame_MIP_DFoF; % bottom right: frame_MIP_DFoF
    frame_out(hplot+1:hplot+hMIP_a,max(wlow,whigh)+wMIP_d+1:max(wlow,whigh)+wMIP_d+wMIP_a,:) = frame_MIP_affine; % bottom middle: frame_MIP_aligned
    writeVideo(obj_out,frame_out);
    
    %delete the marker at each frame such that it won't remain in next
    for iBeh = 1:length(behavior_var)
        h = findobj('type','line','tag',['mark_' behavior_var{iBeh}]);
        delete(h);
    end
    h = findobj('type','text'); % request the texts and delete them
    delete(h);
end
frame_out = uint8(zeros(height_out,width_out,3)); % add a blank frame to the end of each segment to indicate that thers may be a break
writeVideo(obj_out,frame_out);
close(obj_out);
close(fig);
disp('done!');
end