function returnpath = getpath(datatype,sessionID,fishID)
%datatype could be {neural activity,activity}. or {'behavior'}
%sessionID could be '200705'
%fishID could be '1'/'2'/''
[~,computer_name] = system('hostname');
if ispc
    prepath = 'D:';
elseif isunix
    prepath = '/home/hdd1';
end
if ismember(datatype,{'neural activity','activity'})
    returnpath = fullfile(prepath,'Fish-Brain-Behavior-Analysis','results',sessionID);
elseif ismember(datatype,{'behavior'})
    folderpath = fullfile(prepath,'Fish-Brain-Behavior-Analysis','results','behaviour',...
        [sessionID '_fish' num2str(fishID) '*']);
    folderinfo = dir(folderpath);
    returnpath = fullfile(folderinfo.folder,folderinfo.name);
elseif ismember(datatype,{'imaging'})
    returnpath = fullfile(prepath,'Fish-Brain-Behavior-Analysis','results',sessionID,fishID,'1st section');
elseif ismember(datatype,{'data'})
    returnpath = fullfile(prepath,'Mango_Wang','Data');
elseif ismember(datatype,{'code','Code'})
    returnpath = fullfile(prepath,'Mango_Wang','Code');
elseif ismember(datatype,{'result','results'})
    returnpath = fullfile(prepath,'Mango_Wang','Result');
else
    error('wrong query, datatype could be {neural activity,activity}. or {behavior} sessionID could be "200705" and fishID could be "1"/"2"/""');
end
end

