function ax = plot_brain_contour(center,flag)
%function ax = plot_brain_contour(center,flag)
%flag could be '3d' (one scatter3 plot) or '2d' (two scatter subplot)
x = center(:,1);y = center(:,2);
xlim = [min(x) max(x)];
ylim = [min(y) max(y)];
xgrid = xlim(1):15:xlim(2)+5;ygrid = ylim(1):15:ylim(2)+5;
xygrid = cell(length(xgrid)-1,length(ygrid)-1);
labelx = discretize(x,xgrid);labely = discretize(y,ygrid);
for ix = 1:length(xgrid)-1
    for iy = 1:length(ygrid)-1
        idx = find(labelx==ix & labely==iy);
        [ylow,regionlow] = min(center(idx,3));
        [yhigh,regionhigh] = max(center(idx,3));
        if isempty(ylow) ylow = nan; end
        if isempty(yhigh) yhigh = nan; end
        xygrid{ix,iy} = cat(2,repmat([mean(xgrid(ix:ix+1)) mean(ygrid(iy:iy+1))],2,1),[ylow;yhigh]);
    end
end
xygrid =xygrid(:);
xygrid = reshape(permute(cat(3,xygrid{:}),[2 1 3]),3,[])';%(grid*2) x 3
idx = isnan(xygrid(:,3));
xygrid(idx,:) = [];
figure,
if strcmp(flag,'3d')
    scatter3(xygrid(:,1),xygrid(:,2),xygrid(:,3),10,'k');hold on;
elseif strcmp(flag,'2d')
    ax(1) = subplot(1,2,1),
    scatter(xygrid(:,1),xygrid(:,2),10,'k');hold on;
    ax(2) = subplot(1,2,2),
    scatter(xygrid(:,2),xygrid(:,3),10,'k');hold on;
end
end
