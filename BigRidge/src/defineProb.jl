function defineProb(dataset)
    # defining a prob based on the dataset
    X,y  = loadDataset(dataset);
    prob = Prob([],[],[],0,0,dataset);
    n,p = size(X); 
    # centering and scaling the data
    Xmean = mean(X,1);
    Xstd = std(X,1);  
    ind = (0.==Xstd);Xstd[ind] =1.0;  #replace 0 in std by 1 incase there is a constant feature
    X= (X.-Xmean)./Xstd; # Centering and scaling
    
    prob.n=  n;
    lambda = maximum(sum(X.^2,1))/(4.0*p);# This is the dimensionally homogeanous choice for lambda  #1/p ;
    prob.A = X*X'+eye(n)*lambda;
    prob.b = X*y;

    try 
       prob.xsol = load("../data/$(dataset)-xsol.jld", "xsol"
    catch loaderror
       println(loaderror)
       prob.xsol = prob.A\prob.b;
       save("../data/$(dataset)-xsol.jld", "xsol", prob.xsol)
    end
    prob.lambda = lambda;  
    return prob;

end


        
        