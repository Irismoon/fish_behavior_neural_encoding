function [pixelgrid,distgrid,anglegrid] = pixelCoor(param_head_angle_visuo,param_head_dist_visuo,param_speed_visuo)
%function [pixelgrid,distgrid,anglegrid] = pixelCoor(param_head_angle_visuo,param_head_dist_visuo,param_speed_visuo)
anglelim = max(abs(param_head_angle_visuo));
anglegrid = -anglelim:5:anglelim;
anglegrid = [anglegrid anglegrid(end)+5];
distlim = [min(param_head_dist_visuo) max(param_head_dist_visuo)];
distgrid = row2col(log10(distlim(1)):0.1:log10(distlim(2)),1);
distgrid = [distgrid;distgrid(end)+0.1];
label_dist = discretize(log10(param_head_dist_visuo),distgrid);
label_angle = discretize(param_head_angle_visuo,anglegrid);
coor = [label_angle label_dist];
pixelgrid = zeros(length(distgrid),length(anglegrid),length(param_speed_visuo));
for i=1:size(coor,1)
    pixelgrid(coor(i,2),coor(i,1),i) = param_speed_visuo(i);
end
pixelgrid = reshape(pixelgrid,[],length(param_speed_visuo));