function [outputArg1,outputArg2] = select_bouts(inputArg1,inputArg2)
[sessionID,fishID] = getfish;
numfish = length(sessionID);
for iFish=1:numfish
    load(fullfile(getpath('behavior',sessionID{iFish},fishID{iFish}),'high_analysis'),'conv_or_not');
    load(fullfile(getpath('behavior',sessionID{iFish},fishID{iFish}),'tail_swing'),'sum_curv','curvdata','bout_idx','new_centerline_20');
    load(fullfile(getpath('behavior',sessionID{iFish},fishID{iFish}),'low_video_analysis_result'),'param_head_angle_all','param_head_dist_all','no_param');
    if (length(conv_or_not)-length(curvdata))>-5 && (length(conv_or_not)-length(curvdata))<=0
        conv_or_not = [conv_or_not;zeros(length(curvdata)-length(conv_or_not),1)];
    elseif (length(conv_or_not)-length(curvdata))<5 && (length(conv_or_not)-length(curvdata))>0
        conv_or_not(end:end-(length(conv_or_not) - length(curvdata))+1) = [];
    else
        load(fullfile(getpath('behavior',sessionID{iFish},fishID{iFish}),'align_with_fluo'));
        conv_or_not = conv_or_not(align_with_fluo_high==1);
        [curvdata,new_centerline_20,param_head_angle_all,param_head_dist_all,no_param] = samfnmultvar(@(x) x(align_with_fluo_low==1,:,:),curvdata,new_centerline_20,...
            param_head_angle_all,param_head_dist_all,no_param);
    end
    %only select bouts with eye converge
    conv_mask = arrayfun(@(i) mean(conv_or_not(bout_idx(i,1):bout_idx(i,2))),1:size(bout_idx,1))>0;%
    bout_idx = bout_idx(conv_mask,:);
    %only select bouts with prey in the visual field
    prey_mask = arrayfun(@(i) mean(no_param(bout_idx(i,1):bout_idx(i,2))),1:size(bout_idx,1))>0;%
    bout_idx = bout_idx(prey_mask,:);
    %only select bouts responding correctly to the prey when prey is within
    %the visual field
    startFrame = bout_idx(:,1);
    param_head_angle_move = zeros(length(startFrame),1);param_head_dist_move = param_head_angle_move;
    for iframe=1:length(startFrame)
        if startFrame(iframe)>20
            [param_head_angle_move(iframe),param_head_dist_move(iframe)] = samfnmultvar(@(x) mean(x(startFrame(iframe)-20:startFrame(iframe)-1)),...
                param_head_angle_all,param_head_dist_all);
        else
            disp('first frame detected!');
            [param_head_angle_move(iframe),param_head_dist_move(iframe)] = samfnmultvar(@(x) x(startFrame(iframe)),...
                param_head_angle_all,param_head_dist_all);
        end
    end
    [sign_curv,sum_curv_move] = arrayfun(@(i) maxabsk(sum_curv(startFrame(i):startFrame(i)+5),2),1:length(startFrame));
    %we keep forward and corect turns and delete those obviously wrong
    %turns, 
%     [bout_idx,sum_curv_move,param_head_angle_move] = samfnmultvar(@(x) x(sign_curv>0,:),bout_idx,row2col(sum_curv_move,1),row2col(param_head_angle_move));
    %angle<0 when at right, angle>0 when at left, curv<0 towards right
    %while curv>0 towards left
    idx = (row2col(sum_curv_move,1) .* row2col(param_head_angle_move,1))>0;
    good_bout_idx = bout_idx(idx | (row2col(sign_curv,1)<0),:);%same direction between prey and tail swing, or forward
    save(fullfile(getpath('behavior',sessionID{iFish},fishID{iFish}),'tail_swing'),'good_bout_idx','-append');
end
end
function [y,yy] = maxabsk(x,k)
[~,I] = maxk(abs(x),k);
y = prod(x(I));
yy = x(I(1));
end
