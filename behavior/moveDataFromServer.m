function moveDataFromServer(str)
[sessionID,fishID] = getfish;
numfish = length(sessionID);
prepath = 'Z:\\Fish-Brain-Behavior-Analysis\analysis_result';
for i=1:numfish
    folderpath = fullfile(prepath,[sessionID{i} '_fish' num2str(fishID{i}) '*']);
    folderinfo = dir(folderpath);
    returnpath = fullfile(folderinfo.folder,folderinfo.name);
    load(fullfile(returnpath,'low_video_analysis_result.mat'),'centerline_all');
    load(fullfile(returnpath,'tail_swing.mat'),'curvdata');
    centerline = centerline_all;
    save(fullfile(getpath('behavior',sessionID{i},fishID{i}),'tail_swing'),'centerline','curvdata','-append');
end
end

