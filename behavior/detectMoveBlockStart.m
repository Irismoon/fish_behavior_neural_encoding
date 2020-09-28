function startFrame = detectMoveBlockStart(sessionID,fishID)
load(fullfile(getpath('behavior',sessionID,fishID),'tail_swing'),'sum_curv','left_tail_swing','right_tail_swing');
%one movement block is defined as a series of continuous bouts.
swing_or_not = abs(left_tail_swing) + abs(right_tail_swing);
swing_idx = find(swing_or_not);
swing_sequence = diff(swing_idx);
% [y,edge] = histcounts(swing_sequence,1:20:max(swing_sequence));
startFrame = [swing_idx(1);swing_idx(find(swing_sequence>50)+1)];
contin_or_not = arrayfun(@(i) mean(swing_or_not(startFrame(i)+(1:5)))>.5,1:length(startFrame));
%this block at least continues more than 5 frames in next 10 frames
startFrame = startFrame(contin_or_not);
end

