function Countstream(A,s, ll)
    
m,n = size(A);
sgn = reshape(sample(1:2,prob.n) * 2 - 3, 1,prob.n); # one half are +1 and the rest are −1


x = broadcast(*, A, sgn); # flip the signs of each column w.p. 50%


C = zeros(m, s); # initialize C
for j=1:n
    
    C[:, ll[j]] = C[:, ll[j]] +  x[:, j];
    
end
return C;
end