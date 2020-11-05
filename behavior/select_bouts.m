function [outputArg1,outputArg2] = select_bouts(inputArg1,inputArg2)
[sessionID,fishID] = getfish;
numfish = length(sessionID);
for iFish=1:numfish
    load(fullfile(getpath('behavior',sessionID{iFish},fishID{iFish}),'high_analysis'),'conv_or_not');
    load(fullfile(getpath('behavior',sessionID{iFish},fishID{iFish}),'tail_swing'),'angledata','sum_curv','curvdata','bout_idx','new_centerline_20');
    load(fullfile(getpath('behavior',sessionID{iFish},fishID{iFish}),'low_video_analysis_result'),'param_head_angle_all','param_head_dist_all','no_param');
    if (length(conv_or_not)-length(curvdata))>-5 && (length(conv_or_not)-length(curvdata))<=0
        conv_or_not = [conv_or_not;zeros(length(curvdata)-length(conv_or_not),1)];
    elseif (length(conv_or_not)-length(curvdata))<5 && (length(conv_or_not)-length(curvdata))>0
        conv_or_not(end:end-(length(conv_or_not) - length(curvdata))+1) = [];
    else
        load(fullfile(getpath('behavior',sessionID{iFish},fishID{iFish}),'align_with_fluo'));
        conv_or_not = conv_or_not(align_with_fluo_high==1);
        [angledata,sum_curv,curvdata,new_centerline_20,param_head_angle_all,param_head_dist_all,no_param] = samfnmultvar(@(x) x(align_with_fluo_low==1,:,:),...
            angledata,sum_curv,curvdata,new_centerline_20,...
            param_head_angle_all,param_head_dist_all,no_param);
    end
    %only select bouts with eye converge
    conv_mask = arrayfun(@(i) mean(conv_or_not(bout_idx(i,1):bout_idx(i,2))),1:size(bout_idx,1))>0.5;%
    bout_idx = bout_idx(conv_mask,:);
    %only select bouts with prey in the visual field
    prey_mask = arrayfun(@(i) mean(no_param(bout_idx(i,1):bout_idx(i,2))),1:size(bout_idx,1))>0;%
    bout_idx = bout_idx(prey_mask,:);
    %only select bouts responding correctly to the prey when prey is within
    %the visual field
    startFrame = bout_idx(:,1);endFrame = bout_idx(:,2);
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
    %first feature of J-turn
    %consistent swing direction and prey location
    [sign_curv,absmax_sum_curv_move] = arrayfun(@(i) maxabsk(sum_curv(startFrame(i):endFrame(i)),3),1:length(startFrame));
    idx_J = (row2col(absmax_sum_curv_move,1) .* row2col(param_head_angle_move,1))>0 & row2col(sign_curv,1);%same direction with prey and continuous same direction swing
    %2nd feature of J-turn
    %no very large opposite directional swing 
    min_sum_curv_move = arrayfun(@(i) min(sum_curv(startFrame(i):endFrame(i))),1:length(startFrame));
    max_sum_curv_move = arrayfun(@(i) max(sum_curv(startFrame(i):endFrame(i))),1:length(startFrame));
    for iFrame=1:length(startFrame)
    if sign(absmax_sum_curv_move(iFrame))>0&&sign(min_sum_curv_move(iFrame))<0
        idx_J(iFrame) = idx_J(iFrame)&(abs(min_sum_curv_move(iFrame))<0.8);
    elseif sign(absmax_sum_curv_move(iFrame))<0&&sign(max_sum_curv_move(iFrame))>0
        idx_J(iFrame) = idx_J(iFrame)&(max_sum_curv_move(iFrame)<0.8);
    end
    end
    %3rd feature of J-turn
    %no very large left->right/right->left beat between continuous frame
    for iFrame=1:length(startFrame)
        tmp = sum_curv(startFrame(iFrame):endFrame(iFrame));
        tmp = [tmp(1:end-1) tmp(2:end)];
        idx_J(iFrame) = idx_J(iFrame)&(max(abs(tmp(:,2) - tmp(:,1)))<1.5);
    end
    %4th feature of J-turn %ToDo
    %left and right turn angle not too large
%     idx2 = swing_angle<45;
    %1st feature of forward.
    %left/right small-angle beat
    swing_angle = rad2deg(arrayfun(@(i) max(angledata(startFrame(i):endFrame(i),end)-angledata(startFrame(i):endFrame(i),1)),1:length(startFrame)));
    
    idx_forward = (row2col(sign_curv,1)==0&row2col(swing_angle,1)<20);%
    good_bout_idx = bout_idx(idx_J | idx_forward,:);%same direction between prey and tail swing, or forward
    save(fullfile(getpath('behavior',sessionID{iFish},fishID{iFish}),'tail_swing'),'good_bout_idx','-append');
end
end

