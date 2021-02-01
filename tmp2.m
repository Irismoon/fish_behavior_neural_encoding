n = 1000;
k = 200;
A1 = ones(k,k);
A2 = zeros(k,n-k);
A3 = zeros(n-k,k);
A4 = ones(n-k,n-k);

A = [A1 A2;A3 A4];
[eigenVector,eigenValue] = eig(A);
%%
eigenValue = diag(eigenValue);
figure,
plot(eigenValue);