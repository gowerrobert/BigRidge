function hada(x, idx)

# Sx for x a matrix and S our sketch matrix


n  = size(x,1)
p  = size(x,2)
#s = length(idx);


if x  == 0
    return zeros(s,p);
end

if floor(log2(n)) == log2(n)
    
    m = Int(floor(n/2));
    
   # idx_new = idx - m*(idx .> m);
        
    
    if n > 1
        
        x1 = x[1:m,:];
        x2 = x[m+1:end,:];
       # idx1 = idx .* (idx .<= m);
       # idx2 = (idx-m) .* (idx .> m);

        
      #  vect1 = hada(x1,idx_new);
      #  vect2 = hada(x2,idx_new);
            
        
       # return vect1+vect2 -2*vect2.*repmat(idx.>m,1,p);
        
            return hada(x1+x2,idx .* (idx .<= m)) + hada(x1-x2,(idx-m) .* (idx .> m));
        
        
    else
    #    x1 = reshape(x[1,:],1,p);
     #   x2 = reshape(x[2,:],1,p);
      #  return repmat(x1+x2,s,1) - 2*repmat(x2,s,1).*repmat(idx-1,1,p);
       
            return x.*idx;
    end
else
    integ = Int(ceil(log2(n)));
    N = 2^integ;
    X = zeros(N,p);
    X[1:n,:] = x;
        
    return hada(X,idx);
end


end