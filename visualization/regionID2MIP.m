function mip_map = regionID2MIP(value,sessionID,fishID)
%mip_map is 880 x 880 matrix
load(fullfile(getpath('neural activity',sessionID,fishID),'Coherence3'),'A3');
value_max = max(value);
value = row2col(value,1);
value = value/value_max*65535;
map_value = reshape(A3*value,[600 600 280]);
mip_map = [squeeze(max(map_value,[],3)) squeeze(max(map_value,[],2));squeeze(max(map_value,[],1))' zeros(280,280)];
end

