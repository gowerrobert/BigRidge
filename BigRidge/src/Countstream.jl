function Countstream(A,b,)
    
m,n = size(A);


#x = broadcast(*, A, sgn); # flip the signs of each column w.p. 50%


C = zeros(m, s); # initialize C
for j=1:n
    
    C[:, ll[j]] = C[:, ll[j]] +  A[:, j];
    
end
return C;
end