function converted_angle = angle_midline2normal(orig_angle,formflag)
%convert the angle in the fish midline as x axis coordiante into that in
%the normal polar coordinate
disp('default is degree rather than radian!')
if nargin<2
    formflag = 'deg';
end
if strcmp(formflag,'deg')
    orig_angle = deg2rad(orig_angle);
end
orig_angle = -orig_angle;%to correct the original data recording from left is negative angle while right is positive angle to opposite
converted_angle = angle(sin(orig_angle)+1i*cos(orig_angle));
end

