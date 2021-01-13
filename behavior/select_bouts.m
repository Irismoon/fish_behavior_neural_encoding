function [outputArg1,outputArg2] = select_bouts(sessionID,fishID)

load(fullfile(getpath('behavior',sessionID,fishID),'high_analysis'),'conv_or_not');
load(fullfile(getpath('behavior',sessionID,fishID),'tail_swing'),'angledata','sum_curv','curvdata','bout_idx','new_centerline_20','angle_data');
load(fullfile(getpath('behavior',sessionID,fishID),'low_video_analysis_result'),'param_head_angle_all','param_head_dist_all','no_param');
if exist('angle_data','var')
    angledata = angle_data;
end
if (length(conv_or_not)-length(curvdata))>-5 && (length(conv_or_not)-length(curvdata))<=0
    conv_or_not = [conv_or_not;zeros(length(curvdata)-length(conv_or_not),1)];
elseif (length(conv_or_not)-length(curvdata))<5 && (length(conv_or_not)-length(curvdata))>0
    conv_or_not(end:end-(length(conv_or_not) - length(curvdata))+1) = [];
else
    load(fullfile(getpath('behavior',sessionID,fishID),'align_with_fluo'));
    conv_or_not = conv_or_not(align_with_fluo_high==1);
    [angledata,sum_curv,curvdata,new_centerline_20,param_head_angle_all,param_head_dist_all,no_param] = samfnmultvar(@(x) x(align_with_fluo_low==1,:,:),...
        angledata,sum_curv,curvdata,new_centerline_20,...
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
startFrame = bout_idx(:,1);endFrame = bout_idx(:,2);
param_head_angle_move = zeros(length(startFrame),1);param_head_dist_move = param_head_angle_move;
param_move_direction = param_head_angle_move;
for iframe=1:length(startFrame)
    if startFrame(iframe)>10
        [param_head_angle_move(iframe),param_head_dist_move(iframe)] = samfnmultvar(@(x) mean(x(startFrame(iframe)-10:startFrame(iframe)-1)),...
            param_head_angle_all,param_head_dist_all);
         tmp = sign(diff(param_head_angle_all(startFrame(iframe)+[-5:0])));
         if all(tmp==1)
            param_move_direction(iframe) = 1;
         elseif all(tmp==-1)
             param_move_direction(iframe) = -1;
         else
             param_move_direction(iframe) = 0;
         end
    else
        disp('first frame detected!');
        [param_head_angle_move(iframe),param_head_dist_move(iframe)] = samfnmultvar(@(x) x(startFrame(iframe)),...
            param_head_angle_all,param_head_dist_all);
        tmp = sign(diff(param_head_angle_all(1:startFrame(iframe))));
        if all(tmp==1)
            param_move_direction(iframe) = 1;
         elseif all(tmp==-1)
             param_move_direction(iframe) = -1;
         else
             param_move_direction(iframe) = 0;
         end
    end
end
%first feature of J-turn
%consistent swing direction and prey location
[sign_curv,absmax_sum_curv_move] = arrayfun(@(i) sign_maxabsk(sum_curv(startFrame(i):endFrame(i)),2),1:length(startFrame));
idx_J = (row2col(absmax_sum_curv_move,1) .* row2col(param_head_angle_move,1))>0 & row2col(sign_curv,1);%same direction with prey and continuous same direction swing
idx_J = idx_J | ((row2col(absmax_sum_curv_move,1) .* row2col(param_move_direction,1))>0 & row2col(sign_curv,1));%or same direction as the moving direction of preyy
%2nd feature of J-turn
%no very large opposite directional swing
min_sum_curv_move = arrayfun(@(i) min(sum_curv(startFrame(i):endFrame(i))),1:length(startFrame));
max_sum_curv_move = arrayfun(@(i) max(sum_curv(startFrame(i):endFrame(i))),1:length(startFrame));
for iFrame=1:length(startFrame)
    if sign(absmax_sum_curv_move(iFrame))>0&&sign(min_sum_curv_move(iFrame))<0%if it swing to left maximally
        idx_J(iFrame) = idx_J(iFrame)&(abs(min_sum_curv_move(iFrame))<0.8);%its right swing should be small
    elseif sign(absmax_sum_curv_move(iFrame))<0&&sign(max_sum_curv_move(iFrame))>0%if it swing to right maximally
        idx_J(iFrame) = idx_J(iFrame)&(max_sum_curv_move(iFrame)<0.8);%its left swing should be small
    end
end
%3rd feature of J-turn
%no very large left->right/right->left beat between continuous frame
for iFrame=1:length(startFrame)
    tmp = sum_curv(startFrame(iFrame):endFrame(iFrame));
    tmp = [row2col(tmp(1:end-1),1) row2col(tmp(2:end),1)];
    tmp = ~(any( sign(tmp(:,2)).*sign(tmp(:,1))<0   &   abs(tmp(:,2))>0.8   &   abs(tmp(:,1))>0.8 ));
    idx_J(iFrame) = idx_J(iFrame)&tmp;
end
%4th feature of J-turn %ToDo
%left and right turn angle not too large
%     idx2 = swing_angle<45;
%1st feature of forward.
%left/right small-angle beat
swing_angle = rad2deg(arrayfun(@(i) max(angledata(startFrame(i):endFrame(i),end)-angledata(startFrame(i):endFrame(i),1)),1:length(startFrame)));

idx_forward = (row2col(sign_curv,1)==0&row2col(swing_angle,1)<20);%
good_bout_idx = bout_idx(idx_J | idx_forward,:);%same direction between prey and tail swing, or forward
idx = idx_J + 2*idx_forward;
good_bout_idx(:,3) = idx(idx~=0);

m_conv_mask = arrayfun(@(i) mean(conv_or_not(good_bout_idx(i,1)+[-2:0])),1:size(good_bout_idx,1));%
good_bout_idx(:,4) = m_conv_mask<1;%<1 is simultaneous, ==1 is already converged
S = matdoc('comment','good_bout_idx: 3rd column: 1 is turn, 2 is forward; 4th column; true is convergence start, false is already converged');
save(fullfile(getpath('behavior',sessionID,fishID),'tail_swing'),'good_bout_idx','S','-append');
end
function [sign_x,x] = sign_maxabsk(X,k)
 [tmp,I] = maxk(abs(X),k);
 tmp = X(I);
 sign_x = sign(tmp);
 x = tmp(1);
 sign_x = prod(sign_x)>0;
end