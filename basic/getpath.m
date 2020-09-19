function returnpath = getpath(datatype,sessionID,fishID)
%datatype could be {neural activity,activity}. or {'behavior'}
%sessionID could be '200705'
%fishID could be '1'/'2'/''
if ismember(datatype,{'neural activity','activity'})
    returnpath = fullfile('F:\Fish-Brain-Behavior-Analysis\results\',sessionID);
elseif ismember(datatype,{'behavior'})
    folderpath = fullfile('F:\Fish-Brain-Behavior-Analysis\results\behaviour',...
        [sessionID '_fish' num2str(fishID) '*']);
    folderinfo = dir(folderpath);
    returnpath = fullfile(folderinfo.folder,folderinfo.name);
elseif ismember(datatype,{'imaging'})
    returnpath = fullfile('F:\Fish-Brain-Behavior-Analysis\results\imaging',sessionID,fishID);
else
    error('wrong query, datatype could be {neural activity,activity}. or {behavior} sessionID could be "200705" and fishID could be "1"/"2"/""');
end
end

