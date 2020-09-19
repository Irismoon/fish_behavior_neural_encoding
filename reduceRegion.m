function fr_dimReduced = reduceRegion(fr)
%function reduceRegion(fr)
%fr should be time x region shape
%correlation coefficient, and reduce the number of regions
%%
[numt,numRegion] = size(fr);
%%
[rho,p] = corrcoef(fr);
rho = triu(rho,-1);
p = triu(p,-1);
[row,col] = find(abs(rho)>0.5 & p<0.05);
%%
%delete the highly correlated regions
del_region = unique(col);
remain_region = setdiff(1:numRegion,del_region);
%%
disp('run this section');
fr_small = fr(:,remain_region);
[coeff,score,latent,tsquared,explained,mu] = pca(fr_small);
%%
figure,plot(cumsum(explained));
title('acumulated explained variance');
ylabel('percent');
xlabel('index of component');
%%
accum_var = cumsum(explained);
keep_comp_idx = find(accum_var>90);
fr_dimReduced = score(:,1:keep_comp_idx);
%%
figure,
histogram(explained(1:keep_comp_idx),300);
title('explained variance');
%%
figure,
imagesc(coeff(:,1:300)');colorbar;
title('coeff matrix of pca');
xlabel('region');ylabel('pc');
%%
%this is about ppca
[coeff_p,score_p,pcvariance,mu_p,v] = ppca(fr_small,30);
subspace(coeff(:,1:30),coeff_p)
figure,
plot(cumsum(pcvariance)/sum(pcvariance)*100);
%%
%plot the evolution of pc activity with time 
filename = fullfile('F://Mango_Wang/Result','neural_PC evolution.gif');
fig=figure;hold on;
ax = gca;
ax.XLim = [min(score(:,1)),max(score(:,1))];
ax.YLim = [min(score(:,2)),max(score(:,2))];
ax.ZLim = [min(score(:,3)),max(score(:,3))];
% view(0,0);%pc 1 and 3
% view(0,90);%pc 1 and 2
view(90,0);%pc 2 and 3
xlabel('PC 1');ylabel('PC 2');zlabel('PC 3');
for it = 1:numt-1
    plot3([score(it,1) score(it+1,1)],[score(it,2) score(it+1,2)],[score(it,3),score(it+1,3)],'Color','k');
    drawnow;
    frame = getframe(fig);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if it==1
        imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',0);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0);
    end
end
end

