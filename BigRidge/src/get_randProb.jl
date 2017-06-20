function get_randProb(n::Int64)
    # defining a random prob based on gaussian
    prob = Prob([],[],[],0,0,string("rand-",n));
    prob.n=  n;
    X = floor(3*rand(n,n));
    prob.b = X'*rand(n);
    lambda = maximum(sum(X.^2,1))/(4.0*n);
    prob.A = X'*X + lambda*eye(n)
    prob.xsol = prob.A\prob.b;
    prob.lambda = lambda;  
    return prob;

end


        
        