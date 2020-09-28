%%
load(fullfile(getpath('behavior',sessionID,fishID),'high_analysis'),'converge_angle','conv_or_not','focusAngle');
load(fullfile(getpath('behavior',sessionID,fishID),'tail_swing'),'sum_curv','left_tail_swing','right_tail_swing');
load(fullfile(getpath('behavior',sessionID,fishID),'leftright_eyes_to_param_angle_dist'),'left_eye_to_param_angle',...
    'param_lefteye_dist','right_eye_to_param_angle','param_righteye_dist','left_eye_to_param_angle_visible',...
    'right_eye_to_param_angle_visible');
load(fullfile(getpath('behavior',sessionID,fishID),'align_with_fluo'),'align_with_fluo_high','align_with_fluo_low');
[converge_angle,conv_or_not,focusAngle,sum_curv,left_tail_swing,...
    right_tail_swing,left_eye_to_param_angle,param_lefteye_dist,...
    right_eye_to_param_angle,param_righteye_dist] = samfnmultvar(@(x) x(align_with_fluo_high==1),converge_angle,conv_or_not,focusAngle,sum_curv,left_tail_swing,...
    right_tail_swing,left_eye_to_param_angle,param_lefteye_dist,...
    right_eye_to_param_angle,param_righteye_dist);
%%
figure,
color = rgb({'salmon','steelblue'});
subplot(1,3,1)
hold on;
plot(converge_angle,'Color',color(1,:));
ylim = [min(converge_angle) max(converge_angle)];
% scatter(find(conv_or_not),converge_angle(conv_or_not==1),20,...
%     'filled','MarkerEdgeColor','k','MarkerFaceColor',color(1,:),'LineWidth',2);
%%
yyaxis right;
plot(sum_curv,'Color',color(2,:));
swing_or_not = left_tail_swing+right_tail_swing;
% scatter(find(swing_or_not),sum_curv(swing_or_not==1),10,...
%     'filled','MarkerEdgeColor','k','MarkerFaceColor',color(2,:),'LineWidth',1);
title('behavior variables');axis tight;
legend({'converge_angle','sum_curv'},'Interpreter','none');
%%
subplot(1,3,2)
hold on;
plot(left_eye_to_param_angle,'Color',color(1,:));
plot(right_eye_to_param_angle,'-','Color',color(2,:));
title('visual variables--angle');axis tight;
legend({'left eye to param','right eye to param'});
%%
subplot(1,3,3),hold on;
plot(param_lefteye_dist,'-','Color',color(1,:));
plot(param_righteye_dist,'-','Color',color(2,:));
title('visual variables--dist');axis tight;
axis tight;
