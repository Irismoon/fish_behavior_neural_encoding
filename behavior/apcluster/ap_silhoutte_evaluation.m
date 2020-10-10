function [NC,Sil,Silmin,NCfix] = ap_silhoutte_evaluation(M,labels,NC,NCfix,cut)
% function [NC,Sil,Silmin,NCfix] = ap_silhoutte_evaluation(M,labels,NC,NCfix,cut)
%input:
%M must be a distance matrix
%output:
%Sil, mean silhouette across all samples for each clsuter group
%Silmin, different from Sil in terms that it's the minimum mean silhouette
%across clusters in this cluster group
%e.g., for labels=[1 1 1;2 1 2;1 2 2;2 2 3;1 3 4;2 3 4], Sil is dealing all
%samples for each column, whereas
%Silmin(1)=min(mean(silhouette(label(:,1)==1)),mean(label(:,1)==2)),Silmin(2)=min(mean(label(:,2)==1),mean(label(:,2)==2),mean(label(:,2)==3)),
%if one cluster group, i.e., labels(:,i), has one cluster holding less than
%cut samples, this cluster group is deleted since it divide samples into
%too small clusters
mask = tril(true(size(M)),-1);
M = M(mask)';
dim = length(NC);
Sil = zeros(dim,1);
Sildelete = Sil;
Silmin = Sil;
parfor i = 1:dim
    Y = labels(:,i);
    Smax = silhouette([], Y, M);
    dn = isfinite(Smax);
    Sil(i) = mean(Smax(dn));
    [C, Y, dmax] = ind2cluster(Y);
    dmax = min(dmax);
    Sildelete(i) = dmax < cut;
    Q =zeros(length(C),1);
    for j = 1:length(C)
        R = C{j};
        R = Smax(R);
        dn = isfinite(R);
        Q(j) = mean(R(dn));
    end
    Silmin(i) = min(Q);
end

Sildelete = find(Sildelete);
if length(Sildelete) < dim
    Sil(Sildelete) = [];
    Silmin(Sildelete) = [];
    NC(Sildelete) = [];
    NCfix(Sildelete) = [];
end
end

