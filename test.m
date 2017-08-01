% test for getting an idea of the lambda_min of CC^T, where C is the
% concatenation matrix resulting from the newton sketch method

clear all

n = 6;
s = 3;
r = nchoosek(n,s);
idx = nchoosek(1:n,s);
C  = zeros(n,n); % represents C^T * C;
for i=1:r
    D{i} = zeros(s,n); %  construction of I_C_i
    ind = idx(i,:);
    for k = 1:s
    D{i}(k,ind(k)) = 1;
    end
  %  C = C + D{i}'*D{i};
end
A = zeros(n,n);
for i = 1:n
    for j=1:n
        A(i,j)=j+(i-1)*n;
    end
end
A = A + eye(n); % making A invertible
B = zeros(s*r,s*r);
for i=1:r
    for j=1:r
        
        B(1+(i-1)*s:s*i,1+(j-1)*s:s*j) = D{i}*A*D{j}';
    end
end
        
        