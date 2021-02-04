function [bin_center,spike_max_bin] = map_3dTO2d(center,spike_focus)
[N_x,bin_x] = discretize(center(:,1),[min(center(:,1)):5:(max(center(:,1))+5)]);
[N_y,bin_y] = discretize(center(:,2),[min(center(:,2)):5:(max(center(:,2))+5)]);
tmp = arrayfun(@(i) arrayfun(@(j) N_x==i&N_y==j,unique(N_y),'un',0),unique(N_x),'un',0);
tmp = cellfun(@(x) cat(2,x{:}),tmp,'un',0);
tmp = permute(cat(3,tmp{:}),[3 2 1]);%region x y x x
[row,col] = find(any(tmp,3));
bin_center = zeros(size(row,1),2);
bin_region = cell(size(bin_center,1),1);
for i=1:size(bin_center,1)
    bin_center(i,:) = [bin_x(row(i)) bin_y(col(i))];
    bin_region{i} = find(squeeze(tmp(row(i),col(i),:)));
end
[spike_max_bin,I] = arrayfun(@(iT) arrayfun(@(i) max(spike_focus(bin_region{i},iT)),(1:size(bin_region,1))'),1:size(spike_focus,2),'un',0);
spike_max_bin = spike_max_bin{1};%bin x time
tmp = arrayfun(@(i) bin_region{i}(I{1}(i)),1:length(bin_region));
spike_max_bin(:,2) = spike_focus(tmp,2);
end

