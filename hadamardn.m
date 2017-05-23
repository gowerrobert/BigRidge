function result = hadamardn(x)

% HDx

n = length(x);
m = floor(n/2);

if m > 1
    
    x1 = x(1:m);
x2 = x(m+1:end);



    
     vect1 = hadamardn(x1);
     vect2 = hadamardn(x2);
    
    % O(nlog(n)) complexity
       result = [vect1 + vect2; vect1 - vect2];
    
else
    result = [x(1)+x(2);x(1)-x(2)] ;
end

