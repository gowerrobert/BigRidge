% test for getting an idea of the lambda_min of CC^T in the shuffle case, where C is the
% concatenation matrix resulting from the shuffle sketch method

clear all

n = 8;
A = rand(n,n); 
%A = [[1;2;3] [4;5;6] [7;8;9]];
D = zeros(n,n);
perm = perms(1:n);
for i = 1:factorial(n)
   
    P{i} = zeros(n,n);
    vect = perm(i,:);
    for l = 1:n
    P{i}(l,vect(l)) = 1;
    end
    D = D + P{i}'*A*P{i};
end
E = (sum(A(:))-trace(A))*ones(n,n);
E = E + ((n-1)*trace(A) - (sum(A(:))-trace(A)) )*eye(n);
E = factorial(n-2)*E;
D
E