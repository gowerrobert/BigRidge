function hada(A, idx,sgn)

# Sx for x a matrix and S our sketch matrix


n  = size(A,1)
p  = size(A,2)
s = length(idx);

    # This allocates a matrix of size A at every iteration!
x = broadcast(*, A, sgn); # flip the signs of each column w.p. 50%


if x  == 0
    return zeros(s,p);
end

if floor(log2(n)) == log2(n)
    
    m = Int(floor(n/2));
    
    idx_new = idx - m*(idx .> m);
    
    if n > 3
        
        x1 = x[1:m,:];
        x2 = x[m+1:end,:];
        
        
        vect1 = hada(x1,idx_new,sgn[1:m]);
        vect2 = hada(x2,idx_new,sgn[1:m]);
        
        return vect1+vect2 -2*vect2.*repmat(idx.>m,1,p);
        
    elseif n == 3
        
        x1 = x[1:2,:];
        x3 = reshape(x[3,:],1,p);
        vect1 = hada(x1,idx_new,sgn[1:2]);
        vect2 = repmat(x3,s,1);
        
        return vect1+vect2 -2*vect2.*repmat((idx .== 3),1,p);
        
        
    else
        x1 = reshape(x[1,:],1,p);
        x2 = reshape(x[2,:],1,p);
        return repmat(x1+x2,s,1) - 2*repmat(x2,s,1).*repmat(idx-1,1,p);
        
    end
else
    integ = Int(ceil(log2(n)));
    N = 2^integ;
    X = zeros(N,p);
    X[1:n,:] = A;
        sgn1 = zeros(N,1);
        sgn1[1:n] = sgn;
    return hada(X,idx,sgn1);
end


end