function align_behavior_to_neural(sessionID,fishID)
%function align_behavior_to_neural(sessionID,fishID)
T = getBehaviorVariableName();
numVar = height(T);
load(fullfile(getpath('behavior',sessionID,fishID),'align_with_fluo'));
file_out = fullfile(getpath('behavior',sessionID,fishID),'aligned_behavior.mat');
for iVar = 1:numVar
    load(fullfile(getpath('behavior',sessionID,fishID),T.filename{iVar}),T.behavior_variables{iVar});
    if eval(['length(' T.behavior_variables{iVar} ')==length(align_with_fluo_high)'])
        eval([T.behavior_variables{iVar} '=' T.behavior_variables{iVar} '(align_with_fluo_high==1);']);
        eval([T.behavior_variables{iVar} '=' T.behavior_variables{iVar} '(1:5:length(' T.behavior_variables{iVar} '));']);
        if ~(exist(file_out,'file')==2)
            save(file_out,T.behavior_variables{iVar});
        else
            save(file_out,T.behavior_variables{iVar},'-append');
        end
    end
end
end

