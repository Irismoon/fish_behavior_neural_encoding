function filename = behvar2filename(behavior_var)
%function behvar2filename(behavior_var)
%to find the file name of a behavior variable, like the filename of
%left_eye_angle is high_analysis
%e.g., for behavior_var as
%{'left_eye_angle','converge_angle','left_tail_swing','left_start','conv_or_not'},
%filename is {'high_analysis','tail_swing','saccade'}, index is [1 1 2 3 1]
if ~iscell(behavior_var)
    behavior_var = {behavior_var};
end
T = getBehaviorVariableName();
[lia,locb] = cellfun(@(x) ismember(x,T.behavior_variables),behavior_var);
assert(all(lia),'some variables not found');
filename = T.filename(locb);
end

