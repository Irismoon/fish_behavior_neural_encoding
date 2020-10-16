function plot_continuous_bouts(ax,bouts)
%function plot_continuous_bouts(ax,bouts)
%bouts: t x tail segment
bouts = bouts';%tail x t
[numTail,numTime] = size(bouts);
color = parula(numTime);
hold on;
arrayfun(@(i) plot(ax,bouts(:,i),numTail:-1:1,'Color',color(i,:),'LineWidth',5),1:numTime,'un',0);
axis tight;
maxx = max(abs(bouts),[],'all');
set(gca,'XLim',[-maxx maxx]);
end

