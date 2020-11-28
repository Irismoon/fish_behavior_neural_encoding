function idx = removeStrangeVisuoData(no_param_fluo,param_head_dist_fluo,param_speed_fluo)
no_param = row2col(no_param_fluo==0,1);
%param too far
toofar = row2col(param_head_dist_fluo>100,1);
%not moving
notmoving = row2col(param_speed_fluo<9,1);
%
idx = ~(no_param | toofar | notmoving);
end

