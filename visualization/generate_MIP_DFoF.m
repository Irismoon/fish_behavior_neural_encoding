function generate_MIP_DFoF(sessionID,fishID)
imagefile_path = getpath('imaging',sessionID,fishID);
%%
%load DFoF and A3
load(fullfile(getpath('neural activity',sessionID,fishID),'DFoF'));
load(fullfile(getpath('neural activity',sessionID,fishID),'Coherence3'),'A3');
DFoF_keep = DFoF;
DFoF_keep_max = max(DFoF_keep(:));
DFoF_keep = DFoF_keep/DFoF_keep_max*65535;
for it=1:size(DFoF_keep,2)
	temp = A3 * DFoF_keep(:,it);
	temp = reshape(temp,600,600,280);
	MIPs=[max(temp,[],3) squeeze(max(temp,[],2));squeeze(max(temp,[],1))' zeros(size(temp,3),size(temp,3))];
	MIP=uint16(MIPs);
	if it==1
		imwrite(MIP,fullfile(imagefile_path,'MIP_DFoF_Stack.tif'));
	elseif it>1
		imwrite(MIP,fullfile(imagefile_path,'MIP_DFoF_Stack.tif'),'WriteMode','append');
	end
end