function bout_label = swing2whichside(tail_swing,bout_idx)
%in tail_swing, left is -1, right is 1, no move is 0
%if alternating left and right, forward
%if mostly left or right, turn to one side
[nBout,~] = size(bout_idx);
bout_label = zeros(size(tail_swing));
for i=1:nBout
    swingtmp = tail_swing(bout_idx(i,1):bout_idx(i,2));
    label = (nnz(swingtmp<0)-nnz(swingtmp>0));
    if abs(label)/length(swingtmp)<0.1
        bout_label(bout_idx(i,1):bout_idx(i,2)) = 2;
    else
        bout_label(bout_idx(i,1):bout_idx(i,2)) = -sign(label);
    end
end

