function [outputArg1,outputArg2] = tuned2angle(param_head_angle_visuo,spikes_visuo,center,flag)
if length(param_head_angle_visuo)~=size(spikes_visuo,2)
    warning('wrong sample number');
end
assert(size(spikes_visuo,1)==size(center,1),'region number not match!');
anglelim = [min(param_head_angle_visuo) max(param_head_angle_visuo)];
anglegrid = 5*[-25:2:25];
label_angle = discretize(param_head_angle_visuo,anglegrid);
N = histcounts(param_head_angle_visuo,anglegrid);
spikes_anglebin = arrayfun(@(i) nanmedian(spikes_visuo(:,label_angle==i),2),1:length(anglegrid)-1,'un',0);
std_spikes_anglebin = arrayfun(@(i) std(spikes_visuo(:,label_angle==i),[],2),1:length(anglegrid)-1,'un',0);
del_mask = find(N<0.05*length(param_head_angle_visuo));
label_angle(ismember(label_angle,del_mask)) = [];
label_flag = unique(label_angle);
for i=1:length(label_flag)
    label_angle(label_angle==label_flag(i)) = i;
end
anglegrid = arrayfun(@(i) mean(anglegrid(i:i+1)),1:length(anglegrid)-1);
anglegrid(del_mask) = [];
spikes_anglebin(del_mask) = [];
std_spikes_anglebin(del_mask) = [];
N(del_mask) = [];
disp(['Number of samples in each angle bin:' num2str(N)]);
spikes_anglebin = cat(2,spikes_anglebin{:});
std_spikes_anglebin = cat(2,std_spikes_anglebin{:});
figure,imagesc(spikes_anglebin);
set(gca,'XTick',1:1:length(anglegrid),'XTickLabel',floor(anglegrid(1:1:length(anglegrid))));
[fr_angle,prefangle] = max(spikes_anglebin,[],2);
figure,histogram(prefangle);
Npref = accumarray(prefangle,1);
set(gca,'XTick',1:1:length(anglegrid),'XTickLabel',floor(anglegrid(1:1:length(anglegrid))));
xlabel('preferred angle');ylabel('# regions');
figure,
imagesc(corrcoef(spikes_anglebin));colorbar;
set(gca,'XTick',1:1:length(anglegrid),'XTickLabel',floor(anglegrid(1:1:length(anglegrid))));
set(gca,'YTick',1:1:length(anglegrid),'YTickLabel',floor(anglegrid(1:1:length(anglegrid))));
color = colormap(jet(length(anglegrid)));
figure,hold on;
if strcmp(flag,'all')
    idx = unique(prefangle);
    for i=1:length(idx)
        scatter3(center(prefangle==idx(i),1),center(prefangle==idx(i),2),center(prefangle==idx(i),3),10,color(i,:),'filled');
    end
    colormap('jet');cbar = colorbar;cbar.Ticks = linspace(0,1,length(idx));cbar.TickLabels = floor(anglegrid);
    xlabel('x');ylabel('y');zlabel('z');
elseif strcmp(flag,'single')
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
    %     surface(xgrid,ygrid,xygrid(:,3));
    idx = unique(prefangle);
    [row,col] = arrange_subplots(length(idx));
    for i=1:length(idx)
        subplot_tight(row,col,i),
        scatter3(xygrid(:,1),xygrid(:,2),xygrid(:,3),10,'k');hold on;
        scatter3(center(prefangle==idx(i),1),center(prefangle==idx(i),2),center(prefangle==idx(i),3),10,color(i,:),'filled');
        xlabel('x');ylabel('y'),zlabel('z');
        title(num2str(anglegrid(idx(i))));
        view(0,90);
    end
    colormap('jet');cbar = colorbar;cbar.Ticks = linspace(0,1,length(idx));cbar.TickLabels = floor(anglegrid);
    xlabel('x');ylabel('y');zlabel('z');
end
end