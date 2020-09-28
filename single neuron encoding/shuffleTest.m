function [p,h_normal] = shuffleTest(cond1,cond2)
%
num1 = length(cond1);
num2 = length(cond2);
cond1 = row2col(cond1,1);cond2 = row2col(cond2,1);
comb = [cond1;cond2];
numPerm = 2000;
obs_diff = mean(cond1) - mean(cond2);
rand_diff = zeros(numPerm,1);
for iPerm = 1:numPerm
    I = randperm(num1+num2);
    newComb = comb(I);
    newcond1 = newComb(1:num1);
    newcond2 = newComb(num1+1:end);
    rand_diff(iPerm) = mean(newcond1) - mean(newcond2);
end
h_normal = kstest(rand_diff/std(rand_diff));
if obs_diff>0
    p = mean(rand_diff>obs_diff);
elseif obs_diff<0
    p = mean(rand_diff<obs_diff);
end
end

