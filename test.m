% test for getting an idea of the lambda_min of CC^T, where C is the
% concatenation matrix resulting from the newton sketch method

clear all

n = 4;
s = 2;
r = nchoosek(n,s);
idx = nchoosek(1:n,s);
C  = zeros(n,n);
for i=1:r
    D{i} = zeros(s,n);
    ind = idx(i,:);
    D{i}(1,ind(1)) = 1;
    D{i}(2,ind(2)) = 1;
    C = C + D{i}'*D{i};
end

C