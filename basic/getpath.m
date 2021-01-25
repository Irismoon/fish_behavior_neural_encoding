function returnpath = getpath(datatype,sessionID,fishID)
%function returnpath = getpath(datatype,sessionID,fishID)
%datatype could be {neural activity,activity}. or {'behavior'}
%sessionID could be '200705'
%fishID could be '1'/'2'/''
[~,computer_name] = system('hostname');
if ispc
    prepath = 'D:';
    serverpath = 'Z:';
elseif isunix
    prepath = '/home/hdd1';
end
if ismember(datatype,{'neural activity','activity'})
    returnpath = fullfile(prepath,'Fish-Brain-Behavior-Analysis','results',sessionID);
elseif ismember(datatype,{'behavior'})
    folderpath = fullfile(prepath,'Fish-Brain-Behavior-Analysis','results','behaviour',...
        [sessionID '_fish' num2str(fishID) '*']);
    folderinfo = dir(folderpath);
    if length(folderinfo)>1
        idx = input(['There is ' length(folderinfo) ' fishes, please select which fish to use']);
    else
        idx = 1;
    end
    returnpath = fullfile(folderinfo(idx).folder,folderinfo(idx).name);
elseif ismember(datatype,{'imaging'})
    returnpath = fullfile(prepath,'Fish-Brain-Behavior-Analysis','results',sessionID,fishID,'1st section');
elseif ismember(datatype,{'data'})
    returnpath = fullfile(prepath,'Mango_Wang','Data');
elseif ismember(datatype,{'code','Code'})
    returnpath = fullfile(prepath,'Mango_Wang','Code');
elseif ismember(datatype,{'result','results'})
    returnpath = fullfile(prepath,'Mango_Wang','Result');
elseif ismember(datatype,{'server_original_neural','server_orig_neural'})
    folderpath = fullfile(serverpath,'Fish-Brain-Behavior-Analysis',[sessionID '_fish' num2str(fishID) '*']);
    folderinfo = dir(folderpath);
    assert(length(folderinfo)==1,'please specify the fish number');
    returnpath = fullfile(folderinfo.folder,folderinfo.name);
elseif ismember(datatype,{'server_analysis_result','server_analysis','server_behavior'})
    folderpath = fullfile(serverpath,'Fish-Brain-Behavior-Analysis','analysis_result',[sessionID '_fish' num2str(fishID) '*']);
    folderinfo = dir(folderpath);
    returnpath = fullfile(folderinfo.folder,folderinfo.name);
elseif ismember(datatype,{'server_orig'})
    folderpath = fullfile(serverpath,'Fish-Brain-Behavior-Analysis',[sessionID '_fish' num2str(fishID) '*']);
    folderinfo = dir(folderpath);
    returnpath = fullfile(folderinfo.folder,folderinfo.name);
else
    error('wrong query, datatype could be {neural activity,activity}. or {behavior} sessionID could be "200705" and fishID could be "1"/"2"/""');
end
end

