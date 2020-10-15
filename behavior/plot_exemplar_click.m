function plot_exemplar_click(X,sz,bouts)
figure,hold on;
% scatter(X(:,1),X(:,2),sz,1:length(sz),'filled','UserData',1:length(sz),'ButtonDownFcn',{@plot_exemplar,bouts});
arrayfun(@(i) scatter(X(i,1),X(i,2),sz(i),i,'filled','UserData',i,'ButtonDownFcn',{@plot_exemplar,bouts}),1:length(sz),'un',0);
colormap('jet');colorbar;

end
function plot_exemplar(src,env,bouts)
    idx = src.UserData;
    figure,
    plot_continuous_bouts(gca,bouts{idx});
    title(idx);
end

