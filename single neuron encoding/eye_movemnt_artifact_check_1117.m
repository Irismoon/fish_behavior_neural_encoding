%eye movement artifact check
%1.find all frames with single eye movement
%2.average across those frames
%3.compare to the param pos
sessionID = '201117';
fishID = '1';
load(fullfile(getpath('behavior',sessionID,fishID),'high_analysis'),'left_eye_angle','right_eye_angle');
figure,plot(left_eye_angle);hold on;plot(right_eye_angle);
left_eye_angle = left_eye_angle(align_with_fluo_low==1);
left_eye_angle_fluo = left_eye_angle(1:5:end-5);
right_eye_angle = right_eye_angle(align_with_fluo_low==1);
right_eye_angle_fluo = right_eye_angle(1:5:end-5);
left_eye_speed = [0;diff(left_eye_angle_fluo)];
right_eye_speed = [0;diff(right_eye_angle_fluo)];
figure,hold on;
histogram(left_eye_speed);
histogram(right_eye_speed);
%%
%compare the left and right eye movement speed when convergence starts
eye_speed = [arrayfun(@(i) mean(left_eye_speed(conv_idx(i,1)+[-1:1])),1:size(conv_idx,1));arrayfun(@(i) mean(right_eye_speed(conv_idx(i,1)+[-1:1])),1:size(conv_idx,1))]';
figure,
subplot(2,1,1),
histogram(eye_speed(any(bin(:,1:2),2),1));
hold on;
histogram(eye_speed(any(bin(:,1:2),2),2));
title('bin 1-2');
legend({'left eye','right eye'});
subplot(2,1,2),
histogram(eye_speed(any(bin(:,3:5),2),1));
hold on;
histogram(eye_speed(any(bin(:,3:5),2),2));
title('bin 3-5');
xlabel('movement speed');
%%
eye_speed = [left_eye_speed right_eye_speed];
eye_move = abs(eye_speed)>8;
left_eye_move = [false(5,1);eye_move(6:end-5,1)& ~any([eye_move(1:end-10,2) eye_move(2:end-9,2) eye_move(3:end-8,2) eye_move(4:end-7,2) eye_move(5:end-6,2) eye_move(6:end-5,2) eye_move(7:end-4,2) eye_move(8:end-3,2) eye_move(9:end-2,2) eye_move(10:end-1,2)],2);false(5,1)];
right_eye_move = [false(5,1);eye_move(6:end-5,2)& ~any([eye_move(1:end-10,1) eye_move(2:end-9,1) eye_move(3:end-8,1) eye_move(4:end-7,1) eye_move(5:end-6,1) eye_move(6:end-5,1) eye_move(7:end-4,1) eye_move(8:end-3,1) eye_move(9:end-2,1) eye_move(10:end-1,1)],2);false(5,1)];
figure,plot(left_eye_angle_fluo(1:1000));hold on;plot(right_eye_angle_fluo(1:1000));
yyaxis right;
stem(left_eye_move(1:1000));
title('only left eye');
figure,plot(left_eye_angle_fluo(1:1000));hold on;plot(right_eye_angle_fluo(1:1000));
yyaxis right;
stem(right_eye_move(1:1000));
title('only right eye');
[mean(abs(left_eye_speed(left_eye_move))) mean(abs(right_eye_speed(left_eye_move))) mean(abs(left_eye_speed(right_eye_move))) mean(abs(right_eye_speed(right_eye_move)))]
%%
%param position lateralization
left_eye_prey_pos = [arrayfun(@(i) median(param_head_angle_fluo(i+[-5:0])),find(left_eye_move)) arrayfun(@(i) mean(param_head_dist_fluo(i+[-5:0])),find(left_eye_move))];
right_eye_prey_pos = [arrayfun(@(i) median(param_head_angle_fluo(i+[-5:0])),find(right_eye_move)) arrayfun(@(i) mean(param_head_dist_fluo(i+[-5:0])),find(right_eye_move))];
figure,
polarscatter(deg2rad(left_eye_prey_pos(:,1)), left_eye_prey_pos(:,2), 10,'r');
hold on;
polarscatter(deg2rad(right_eye_prey_pos(:,1)), right_eye_prey_pos(:,2), 10,'b');
legend({'only left eye','only right eye'});
disp('comparing angle:');
ranksum(left_eye_prey_pos(:,1),right_eye_prey_pos(:,1))
disp('comparing distance:');
ranksum(left_eye_prey_pos(:,2),right_eye_prey_pos(:,2))
%%
figure,
histogram(left_eye_prey_pos(:,1),'BinWidth',10);hold on;
histogram(right_eye_prey_pos(:,2),'BinWidth',10);
%select only 10-60 degree
ranksum(left_eye_prey_pos(left_eye_prey_pos(:,1)<60&left_eye_prey_pos(:,1)>10,1),right_eye_prey_pos(right_eye_prey_pos(:,1)<60&right_eye_prey_pos(:,1)>10,1))
ranksum(left_eye_prey_pos(left_eye_prey_pos(:,1)<60&left_eye_prey_pos(:,1)>10,2),right_eye_prey_pos(right_eye_prey_pos(:,1)<60&right_eye_prey_pos(:,1)>10,2))
left_eye_move = find(left_eye_move);
right_eye_move = find(right_eye_move);
mask_left = left_eye_prey_pos(:,1)<60&left_eye_prey_pos(:,1)>10;
mask_right = right_eye_prey_pos(:,1)<60&right_eye_prey_pos(:,1)>10;
left_eye_prey_pos = left_eye_prey_pos(mask_left,:);
right_eye_prey_pos = right_eye_prey_pos(mask_right,:);
left_eye_move = left_eye_move(mask_left);
right_eye_move = right_eye_move(mask_right);
left_eye_move = left_eye_move(randperm(length(left_eye_move),length(right_eye_move)));
%%
left_eye_spike = arrayfun(@(i) mean(spike_raw(:,i+[-5:5]) - mean(spike_raw(:,i+[-8:-5]),2),2),left_eye_move,'un',0);
left_eye_spike = mean(cat(2,left_eye_spike{:}),2);
right_eye_spike = arrayfun(@(i) mean(spike_raw(:,i+[-5:5]) - mean(spike_raw(:,i+[-8:-5]),2),2),right_eye_move,'un',0);
right_eye_spike = mean(cat(2,right_eye_spike{:}),2);
tmp = zeros(size(sMatrix_total,1),1);
tmp(remain_region) = abs(left_eye_spike);
generate_MIP(tmp,'201117','1','_only_left');
tmp(remain_region) = abs(right_eye_spike);
generate_MIP(tmp,'201117','1','_only_right');