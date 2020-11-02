%
n = 100;
fish_pos = randn(n,2);
fish_direction = rand(n,1)*2*pi;
prey_pos = arrayfun(@(i) generate_prey_location(fish_pos(i,:)',fish_direction(i)),1:n,'un',0);
prey_pos = cat(2,prey_pos{:})';
%%
figure,
for i=1:n
    clf;
    quiver(fish_pos(i,1)-cos(fish_direction(i)),fish_pos(i,2)-sin(fish_direction(i)),cos(fish_direction(i)),sin(fish_direction(i)));
    hold on;
    scatter(prey_pos(i,1),prey_pos(i,2),20,'filled');
    set(gca,'XLim',[min(fish_pos(:,1)) max(fish_pos(:,1))]);
    set(gca,'YLim',[min(fish_pos(:,2)) max(fish_pos(:,2))]);
    pause; 
end