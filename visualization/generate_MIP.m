function generate_MIP(X,sessionID,fishID,postfix)
imagefile_path = getpath('result');
if ~exist('postfix','var')
    postfix = '';
end
X = row2col(X,min(size(X)));
%%
%load DFoF and A3
if strcmp(sessionID,'200705')
    load(fullfile(getpath('neural activity',sessionID,fishID),'static_seg'),'A');
elseif strcmp(sessionID,'200108')
    load(fullfile(getpath('neural activity',sessionID,fishID),'Coherence3_correct'),'A3');
    A  = A3;
end
X_keep = X;
if min(X_keep)<0
    warning('The negative number will be treat as zero!');
end
X_keep_max = max(X_keep(:));
X_keep = X_keep/X_keep_max*65535;
for it=1:size(X_keep,2)
	temp = A * X_keep(:,it);
	temp = reshape(temp,600,600,280);
	MIPs=[max(temp,[],3) squeeze(max(temp,[],2));squeeze(max(temp,[],1))' zeros(size(temp,3),size(temp,3))];
	MIP=uint16(MIPs);
	if it==1
		imwrite(MIP,fullfile(imagefile_path,[sessionID ' MIP' postfix '.tif']));
	elseif it>1
		imwrite(MIP,fullfile(imagefile_path,[sessionID ' MIP' postfix '.tif']),'WriteMode','append');
	end
end