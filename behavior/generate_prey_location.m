function prey_pos_abs = generate_prey_location(fish_head_pos,fish_direction)
%the ogrigin of the coordianate is fish's head, and its direction is its
%angle (from tail to head as the vector) relative to the xaxis of the
%image(from left to right as a vector)
persistent k
if isempty(k)
    k = 0;
end

delta = .1;%control the velocity of the prey
%generate a prey at the right of fish and moving downwards
prey_pos_fish = [1;1-k*delta];
if prey_pos_fish(2)<-3
    k = 0;
end
%generate the absolute position of prey
prey_pos_abs = [sin(fish_direction) cos(fish_direction);-cos(fish_direction) sin(fish_direction)]*prey_pos_fish + fish_head_pos;
k = k+1;
end

