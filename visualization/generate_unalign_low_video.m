function [outputArg1,outputArg2] = generate_unalign_low_video(sessionID,fishID,flag)
prepath = fullfile(getpath('server_behavior',sessionID,fishID),'*.avi');
fileinfo = dir(prepath);
t = arrayfun(@(i) regexp(fileinfo(i).name,'dpf_\d{2}_low','match'),1:length(fileinfo),'un',0);
t = cellfun(@(x) str2double(x{1}(5:6)),t);
[t,I] = sort(t,'ascend');
fileinfo = fileinfo(I);

writer = VideoWriter(fullfile(getpath('behavior',sessionID,fishID),'low_unalign.avi'));
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