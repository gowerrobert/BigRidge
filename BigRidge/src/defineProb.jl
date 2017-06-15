function defineProb(dataset,method_name)
    
    # defining a prob based on the name of the chosen method and the dataset

X,y  = loadDataset(dataset);

prob = Prob([],[],[],0,0,dataset);

    n = minimum(size(X)); #the number of data points is not the minimum dimension of X.
p = maximum(size(X));
lambda = 1/p;

if method_name == "Hadamard"
    
    N = Int(2.^(1+floor(log2(n)))); # make the dimension of our problem a power of 2
    mat = eye(N)*lambda;
    mat[1:n,1:n] = X*X'+eye(n)*lambda; #What about X? This doesn't change the dimensions of X.
    prob.A = mat;
    
    vect = zeros(N);
    vect[1:n] = X*y;
    prob.b = vect;
    
    prob.xsol = prob.A\prob.b;
    prob.lambda = lambda;
    
    prob.n = N;
    
    
    
    
else
    
    sX = size(X);
    prob.n=  sX[1];
    lambda = 1/sX[2] ;
    prob.A = X*X'+eye(prob.n)*lambda;
    
    prob.b = X*y;
    prob.xsol = prob.A\prob.b;
    prob.lambda = lambda;
    
end

return prob;

end


        
        