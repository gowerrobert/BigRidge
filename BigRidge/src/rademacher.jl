function boot_rademacher(prob::Prob,options::MyOptions)
    flopsperiter = (options.sketchsize)^3; # Re-think this
    name = "rademacher";
    stepmethod = step_rademacher
    method = Method(flopsperiter,name,step_rademacher,boot_rademacher)
    return method;
end

#function step_countmin(as::Array{Float64,2},bs::Array{Float64,2},sas::Array{Float64,2},S::Array{Float64,2},prob::Prob, x::Array{Float64}, options::MyOptions )

function step_rademacher(prob::Prob, x::Array{Float64}, options::MyOptions )

    s = options.sketchsize;
    rho = convert(Int64,floor(n/s)); #hard coded density of rows

    SA = zeros(s,prob.n); 
   # ind = sample(1:prob.n,prob.n,replace=false);
   # ind = reshape(ind,s,convert(Int64,n/s))
    ind =  Array{Int64}(s,rho);#zeros(s,rho); # matrix of indices
    for i =1:s
         ind[i,:] = sample(1:prob.n,rho,replace=false);
    end

    SA = zeros(s,prob.n);
    Sb = zeros(s);
    for i =1:s
          SA[i,:] = sum(prob.A[ind[i,:],:],1);
          Sb[i] = sum(prob.b[ind[i,:]])
    end
    SAS = zeros(s,s); 
    for i =1:s
         SAS[i,:] = sum(SA[:,ind[i,:]],2);
    end
    vect = SA*x-Sb; # 
 #   try
    y = SAS\vect;   # solving (S^TAS) y = (S^TAx-S^Tb)  
  #  catch somethgsdssd

 #   end
    for i =1:s #adding on S^T y
       x[ind[i,:]] = x[ind[i,:]]-y[i];
    end


end
        
        
        
        
        
        
        
