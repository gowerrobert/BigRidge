function result = hadamards(x,idx)

% SHDx for x a vectot and S our sketch matrix


n = length(x);
m = floor(n/2);


if m > 1
    
idx1 = idx( idx <= m );
idx2 = idx( idx > m) - m;
lx = length(idx1);
    x1 = x(1:m);
    x2 = x(m+1:end);
    vect1 = hadamards(x1,[idx1;idx2]);
    vect2 = hadamards(x2,[idx1;idx2]);
   
    
    result = [vect1(1:lx) + vect2(1:lx); vect1(lx+1:end) - vect2(lx+1:end)];
    
else
    vect = [x(1)+x(2);x(1)-x(2)] ;
    result = vect(idx);
end
