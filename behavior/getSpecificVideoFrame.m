function frame_out = getSpecificVideoFrame(filename,numFrame)
vidObj = VideoReader(filename);
k=0;
while hasFrame(vidObj)
    k = k+1;
    tmp = readFrame(vidObj);
    if k==numFrame
        frame_out = tmp;
        return;
    end 
end
end

