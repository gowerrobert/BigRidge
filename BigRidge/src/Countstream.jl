function Countstream(A,s, ll)
m,n = size(A);

    
C = zeros(m, s); # initialize C 
for j=1:n
#   sgn = 2*sample(1:2,1)-3;

        C[:, ll[j]] = C[:, ll[j]] +  A[:, j];

        # C[:, ll[j]] = C[:, ll[j]] + sgn.*A[:, j];
end
    return C;
end