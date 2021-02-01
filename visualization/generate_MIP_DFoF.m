function generate_MIP_DFoF(sessionID,fishID)
%function generate_MIP_DFoF(sessionID,fishID)
%generate DFoF MIP throughout the recording
imagefile_path = getpath('imaging',sessionID,fishID);
%%
%load DFoF and A3
load(fullfile(getpath('neural activity',sessionID,fishID),'DFoF'));
load(fullfile(getpath('neural activity',sessionID,fishID),'Coherence3'),'A3');
DFoF_keep = DFoF;
if min(DFoF_keep(:))<0
    DFoF_keep = DFoF_keep + abs(min(DFoF_keep(:)));
end
DFoF_keep(isnan(DFoF_keep)) = 0;
DFoF_keep(isinf(DFoF_keep)) = 0;
% if isinf(DFoF_keep_max)
%     DFoF_keep_tmp = DFoF_keep(:);
%     [DFoF_keep_tmp_descend,idx] = sort(DFoF_keep_tmp,'descend');
%     DFoF_keep(isinf(DFoF_keep)) = DFoF_keep_tmp_descend(find(~isinf(DFoF_keep_tmp_descend),1));
% end
% DFoF_keep = DFoF_keep.^0.5;
[N,bin] = histcounts(DFoF_keep(:),'BinWidth',1);
thld = bin(find((cumsum(N)/sum(N))>0.99,1)+1);
DFoF_keep(DFoF_keep>thld) = thld;
DFoF_keep_max = max(DFoF_keep(:));
DFoF_keep = DFoF_keep/DFoF_keep_max*65535;
while nnz(DFoF_keep>50000)<100   
    DFoF_keep(DFoF_keep>50000) = 0;
    DFoF_keep_max = max(DFoF_keep(:));
    DFoF_keep = DFoF_keep/DFoF_keep_max*65535;
end
if ismember(sessionID,{'201117'})
    sz = [376 308 210];
else
    sz = [600 600 280];
end

contour = reshape(full(sum(A3,2)~=0),sz);
contour = [max(contour,[],3) squeeze(max(contour,[],2));squeeze(max(contour,[],1))' zeros(size(contour,3),size(contour,3))];
for it=1:size(DFoF_keep,2)
    if mod(it,500)==0 disp('500 frames done!'); end
    temp = A3 * DFoF_keep(:,it);
    temp = reshape(temp,sz(1),sz(2),sz(3));
    MIPs=[max(temp,[],3) squeeze(max(temp,[],2));squeeze(max(temp,[],1))' zeros(size(temp,3),size(temp,3))];
    MIPs = MIPs + max(MIPs(:))*0.001*contour;
    MIP=uint16(MIPs);
    if it==1
        % 		imwrite(MIP,fullfile(imagefile_path,'MIP_DFoF_Stack.tif'));
        imwrite(MIP,'D:\MIP_DFoF_Stack.tif');
    elseif it>1
        % 		imwrite(MIP,fullfile(imagefile_path,'MIP_DFoF_Stack.tif'),'WriteMode','append');
        imwrite(MIP,'D:\MIP_DFoF_Stack.tif','WriteMode','append');
    end
end