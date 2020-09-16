function reduceRegion(fr)
%function reduceRegion(fr)
%fr should be time x region shape
%correlation coefficient, and reduce the number of regions
%%
[rho,p] = corrcoef(fr);
rho = triu(rho,-1);
p = triu(p,-1);
[row,col] = find(abs(rho)>0.5 & p<0.05);
%%
%delete the highly correlated regions
del_region = unique(col);
end

