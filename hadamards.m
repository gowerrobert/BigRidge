function [result] = hadamards(x,idx)

% SHDx for x a vectot and S our sketch matrix

n = length(x);

m = floor(n/2);



if m > 1
    
    x1 = x(1:m);
x2 = x(m+1:end);

%idx_new = idx;
%idx_new(idx > m) = idx(idx > m)-m;
idx_new = idx.*(idx<=m)+(idx-m).*(idx>m);

    
       vect1 = hadamards(x1,idx_new);
     vect2 = hadamards(x2,idx_new);
    
      % result = (vect1+vect2).*(idx<=m)+(vect1-vect2).*(idx>m) ;
          result = vect1+vect2 -2*vect2.*(idx>m) ;

else
    result = x(1) + x(2) - 2*x(2)*(idx-1) ;
end

