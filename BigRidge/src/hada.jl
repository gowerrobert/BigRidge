function hada(x, idx)

# SHDx for x a matrix and S our sketch matrix

n,p = size(x);
s = length(idx);

m = Int(floor(n/2));

   # idx_new = idx.*(idx.<=m)+(idx-m).*(idx.>m);

    idx_new = idx - m*(idx .> m); 

    
if n > 3
    
    x1 = x[1:m,:];
    x2 = x[m+1:end,:];
    
    
    
    vect1 = hada(x1,idx_new);
    vect2 = hada(x2,idx_new);
    
    return vect1+vect2 -2*vect2.*repmat(idx.>m,1,p);
#     return vect1+vect2 -2*vect2.*broadcast(*,ones(1,p),idx.>m);


        
    elseif n == 3
        
        x1 = x[1:2,:];
        x3 = reshape(x[3,:],1,p);
        vect1 = hada(x1,idx_new);
        vect2 = repmat(x3,s,1);
        
        return vect1+vect2 -2*vect2.*repmat((idx .== 3),1,p);
#         return vect1+vect2 -2*vect2.*broadcast(*,ones(1,p),idx .== 3) ;

       
    else 
        x1 = reshape(x[1,:],1,p);
       x2 = reshape(x[2,:],1,p);
      return repmat(x1+x2,s,1) - 2*repmat(x2,s,1).*repmat(idx-1,1,p);
#       return repmat(x1+x2,s,1) - 2*repmat(x2,s,1).*broadcast(*,ones(1,p),idx-1);


        
end

end