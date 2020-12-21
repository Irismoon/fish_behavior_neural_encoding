function [outputArg1,outputArg2] = generate_original_low_video(sessionID,fishID)
prepath = fullfile(getpath('server_orig_neural',sessionID,fishID),'fish_track_*');
fileinfo = dir(prepath);
t = arrayfun(@(i) regexp(fileinfo(i).name,'\s\d{2}\s\d{2}\s\d{2}','match'),1:length(fileinfo),'un',0);
t = cellfun(@(x) str2double(strrep(x,' ','')),t);
[t,I] = sort(t,'ascend');
fileinfo = fileinfo(I);

writer = VideoWriter(fullfile(getpath('behavior',sessionID,fishID),'low_video.avi'));
open(writer);
for ifile = 1:length(fileinfo)
    disp(['processing file ' num2str(ifile)]);   
    filepath = fullfile(fileinfo(ifile).folder,fileinfo(ifile).name);
    info = imfinfo(filepath);
    imag = arrayfun(@(i) imread(filepath,'Index',i,'Info',info),1:length(info),'un',0);
    imag = uint8(cat(4,imag{:}));
    writeVideo(writer,imag);
end
close(writer);
end

