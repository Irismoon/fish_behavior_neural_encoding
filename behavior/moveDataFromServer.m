function moveDataFromServer(sessionID,fishID)
prepath = 'Z:\\Fish-Brain-Behavior-Analysis\analysis_result';
folderpath = fullfile(prepath,[sessionID '_fish' num2str(fishID) '*']);
folderinfo = dir(folderpath);
returnpath = fullfile(folderinfo.folder,folderinfo.name);
load(fullfile(returnpath,'low_video_analysis_result.mat'),'centerline_all');
load(fullfile(returnpath,'tail_swing.mat'),'curvdata');
centerline = centerline_all;
save(fullfile(getpath('behavior',sessionID,fishID),'tail_swing'),'centerline','curvdata','-append');
end

