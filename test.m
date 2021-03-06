% test for getting an idea of the lambda_min of CC^T, where C is the
% concatenation matrix resulting from the count-min sketch method

clear all

n = 5;
s = 2;
%r = nchoosek(n,s);

r = s^n;

%idx = nchoosek(1:n,s);
idx = permn(1:s,n);

C  = zeros(n,n); % represents C^T * C;
id = eye(s);
for i=1:r
    D{i} = zeros(s,n); %  construction of I_C_i
    ind = idx(i,:);
    for j = 1:n
    D{i}(:,j) = id(:,ind(j));
    end
    C = C + D{i}'*D{i};
end
C = C/(s^(n-1))
% A = zeros(n,n);
% for i = 1:n
%     for j=1:n
%         A(i,j)=j+(i-1)*n;
%     end
% end
% A = A + eye(n); % making A invertible
% B = zeros(s*r,s*r);
% for i=1:r
%     for j=1:r
%         
%         B(1+(i-1)*s:s*i,1+(j-1)*s:s*j) = D{i}*A*D{j}';
%     end
% end
        
        