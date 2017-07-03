function hadamardsVect(x, idx)

# SHDx for x a vectot and S our sketch matrix
    
    n = length(x);
    
   # n = Int(2.^(1+floor(log2(n1))));

   # A = zeros(n,p);
  #  A[1:n1,:] = AA;

m = Int(floor(n/2));


if m > 1
    
        x1 = x[1:m];
    x2 = x[m+1:end];
    
    
    idx_new = idx.*(idx.<=m)+(idx-m).*(idx.>m);
  #  
    vect1 = hadamardsVect(x1,idx_new);
    vect2 = hadamardsVect(x2,idx_new);
    
    return vect1+vect2 -2*vect2.*(idx.>m);
    
else
      #  x1 = reshape(A[1,:],1,p);
     #   x2 = reshape(A[2,:],1,p);
  #  return repmat(A1+A2,s,1) - 2*repmat(A2,s,1).*repmat(idx-1,1,p);
 return x[1] + x[2] - 2*x[2].*(idx-1);
        
end

end

function hadamards(A,idx)
    
    if ndims(A) == 1
        return hadamardsVect(A,idx);
    else
    
    n,p = size(A);
    s = length(idx);
    tab = zeros(s,p);

for j =1:p
        x = A[:,j];
        
        tab[:,j] = hadamardsVect(x,idx);
    end
    return tab;
    end
end

