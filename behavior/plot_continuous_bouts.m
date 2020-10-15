function plot_continuous_bouts(ax,bouts)
%function plot_continuous_bouts(ax,bouts)
%bouts: t x tail segment
bouts = bouts';%tail x t
[~,numT] = size(bouts);
plot(ax,bouts+(1:numT),'Color','k');
axis tight;
end

