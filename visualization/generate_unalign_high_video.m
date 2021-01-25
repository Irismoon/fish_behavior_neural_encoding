function [outputArg1,outputArg2] = generate_unalign_high_video(sessionID,fishID)
prepath = fullfile(getpath('server_orig',sessionID,fishID),'*.avi');
fileinfo = dir(prepath);
t = arrayfun(@(i) regexp(fileinfo(i).name,'\d{2}_processed','match'),1:length(fileinfo),'un',0);
idx = cellfun(@(x) ~isempty(x),t);
t = t(idx);
fileinfo = fileinfo(idx);
t = cellfun(@(x) str2double(x{1}(1:2)),t);
[t,I] = sort(t,'ascend');
fileinfo = fileinfo(I);

writer = VideoWriter(fullfile(getpath('behavior',sessionID,fishID),'high_unalign.avi'));
open(writer);
for ifile = 1:length(fileinfo)
    disp(['processing file ' num2str(ifile)]);   
    filepath = fullfile(fileinfo(ifile).folder,fileinfo(ifile).name);
    V = VideoReader(filepath);
    while hasFrame(V)
        writeVideo(writer,readFrame(V));
    end
end
close(writer);
end