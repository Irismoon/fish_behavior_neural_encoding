function align_all_with_param(sessionID,fishID)
load(fullfile(getpath('behavior',sessionID,fishID),'aligned_behavior'));
C = whos();
file_out = fullfile(getpath('behavior',sessionID,fishID),'with_param_aligned_behavior.mat');
numPoint = length(no_param);
for iVar = 1:length(C)
    if eval(['abs(length(' C(iVar).name ')-numPoint)<0.001'])
        eval([C(iVar).name '=' C(iVar).name '(no_param==1);']);
        if ~(exist(file_out,'file')==2)
            save(file_out,C(iVar).name);
        else
            save(file_out,C(iVar).name,'-append');
        end
    end
end
end

