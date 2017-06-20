function Countstream(A, s)
m,n = size(A);
sgn = reshape(sample(1:2,n) * 2 - 3, 1,n); # one half are +1 and the rest are −1
#A = broadcast(*, A, sgn); # flip the signs of each column w.p. 50%
ll = sample(1:s, n); # sample n items from [s] with replacement
C = zeros(m, s); # initialize C 
for j=1:n
    C[:, ll[j]] = C[:, ll[j]] + A[:, j];
end
    return C;
end