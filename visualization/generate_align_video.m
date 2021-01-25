function generate_low_align_video(sessionID,fishID)
vidLow = VideoReader(fullfile(getpath('behavior',sessionID,fishID),'low_all.avi'));
vidHigh = VideoReader(fullfile(getpath('behavior',sessionID,fishID),'high_all.avi'));

load(fullfile(getpath('behavior',sessionID,fishID),'align_with_fluo'));

vidLow_writer = VideoWriter(fullfile(getpath('behavior',sessionID,fishID),'low_align_with_fluo.avi'));
vidHigh_writer = VideoWriter(fullfile(getpath('behavior',sessionID,fishID),'high_align_with_fluo.avi'));

kLow = 0;
kHigh = 0;
open(vidLow_writer);open(vidHigh_writer);
 while hasFrame(vidLow)
     frame_tmp = readFrame(vidLow);
     kLow = kLow+1;
     if align_with_fluo_low(kLow)==1
        writeVideo(vidLow_writer,frame_tmp);
    end
end

 while hasFrame(vidHigh)
     frame_tmp = readFrame(vidHigh);
     kHigh = kHigh+1;
     if align_with_fluo_high(kHigh)==1
        writeVideo(vidHigh_writer,frame_tmp);
     end
end
close(vidLow_writer);close(vidHigh_writer);