function plot_continuous_bouts(ax,bouts)
%function plot_continuous_bouts(ax,bouts)
%bouts: t x tail x 2 segment
bouts = permute(bouts,[2 1 3]);%tail x t x 2
bouts(:,:,1) = bouts(:,:,1) - mean2(bouts(:,:,1));
[numTail,numTime,numSide] = size(bouts);
color = parula(numTime);
hold on;
if numSide==1
    arrayfun(@(i) plot(ax,bouts(:,i),numTail:-1:1,'Color',color(i,:),'LineWidth',5),1:numTime,'un',0);
    maxx = max(abs(bouts),[],'all');
elseif numSide==2
    arrayfun(@(i) plot(ax,bouts(:,i,1),bouts(:,i,2),'Color',color(i,:),'LineWidth',5),1:numTime,'un',0);
    maxx = max(abs(bouts(:,:,1)),[],'all');
else
    warning('page are too much！');
end
axis tight;
set(gca,'XLim',[-maxx maxx]);
end

