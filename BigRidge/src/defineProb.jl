function defineProb(dataset)
    
    # defining a prob based on the name of the chosen method and the dataset

X,y  = loadDataset(dataset);

prob = Prob([],[],[],0,0,dataset);

n,p = size(X); 
lambda = 1/p;

# if method_name == "Hadamard"
#     
#     N = Int(2.^(1+floor(log2(n)))); # make the dimension of our problem a power of 2
#     mat = eye(N)*lambda;
#     mat[1:n,1:n] = X*X'+eye(n)*lambda; #What about X? This doesn't change the dimensions of X.
#     prob.A = mat;
#     
#     vect = zeros(N);
#     vect[1:n] = X*y;
#     prob.b = vect;
#     
#     prob.xsol = prob.A\prob.b;
#     prob.lambda = lambda;
#     
#     prob.n = N;
#     
#     
#     
#     
# else
    
    prob.n=  n;
    lambda = maximum(sum(X.^2,1))/(4.0*p);# This is the dimensionally homogeanous choice for lambda  #1/p ;
    prob.A = X*X'+eye(n)*lambda;
    
    prob.b = X*y;
    prob.xsol = prob.A\prob.b;
    prob.lambda = lambda;
    

return prob;

end


        
        