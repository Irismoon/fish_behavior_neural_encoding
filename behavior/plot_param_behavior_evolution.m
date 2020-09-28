function plot_param_behavior_evolution(aligned_param_head_dist,aligned_param_head_angle,move_dir)
%%
if max(aligned_param_head_angle)>7
    aligned_param_head_angle = deg2rad(aligned_param_head_angle);
end
figure,hold on;
ax=gca;
ax.XLim = [-83,83];ax.YLim = [-83,83];
for iPoint = 1:length(move_dir)
    scatter(aligned_param_head_dist(iPoint)*cos(aligned_param_head_angle(iPoint)),aligned_param_head_dist(iPoint)*sin(aligned_param_head_angle(iPoint)),'r*');
    scatter(aligned_param_head_dist(iPoint+1)*cos(aligned_param_head_angle(iPoint+1)),aligned_param_head_dist(iPoint+1)*sin(aligned_param_head_angle(iPoint+1)),'bs');
    quiver(aligned_param_head_dist(iPoint)*cos(aligned_param_head_angle(iPoint)),aligned_param_head_dist(iPoint)*sin(aligned_param_head_angle(iPoint)),...
        real(move_dir(iPoint)),imag(move_dir(iPoint)));
    title(['frame ' num2str(iPoint)]);
    pause;
    cla;
end
end

