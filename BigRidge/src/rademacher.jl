function boot_rademacher(prob::Prob,options::MyOptions)
    s  =options.sketchsize;                   
    flopsperiter = s^3+  #  solving  SAS\(SA*x-Sb)  
    (s+1)*((prob.n+1+s)*convert(Int64,ceil(log(floor(prob.n/s)))))+ #  calculating sketch SA, Sb and SAS
    prob.n*(s+1); #computing SA*x-Sb 
    name = string("rademacher-",s);
    if(mod(prob.n,s)!=0) 
        s= s+1;
    end
    SA = zeros(s,prob.n); # saving space for sketched matrix
    stepmethod = step_rademacher;
    method = Method(flopsperiter,name,step_rademacher,boot_rademacher,SA)
    return method;
end

#function step_countmin(as::Array{Float64,2},bs::Array{Float64,2},sas::Array{Float64,2},S::Array{Float64,2},prob::Prob, x::Array{Float64}, options::MyOptions )

function step_rademacher(prob::Prob, x::Array{Float64}, options::MyOptions, method::Method )

    s = options.sketchsize;
    rho = convert(Int64,floor(prob.n/s)); #hard coded density of rows
    ind = sample(1:prob.n,prob.n,replace=false); # shuffle all indices
    modns = mod(prob.n,s);
    divi = prob.n-modns;
    sold =s;
    if(modns!=0) 
        s= s+1;
    end
    indM = Array{Int64}(s,rho);# build a matrix indices with approx equally distributed number of indices over s rows
    indM[1:sold,:] =reshape(ind[1:divi], sold,rho);
    if(modns!=0)
        indM[s,1:modns] = ind[divi+1:end];
        indM[s,modns+1:end] =sample(1:prob.n, rho-modns,replace=false);
    end
    
#   ind =  Array{Int64}(s,rho);#zeros(s,rho); # matrix of indices
#   for i =1:s
#        ind[i,:] = sample(1:prob.n,rho,replace=false);
#   end
  # SA = zeros(s,prob.n);  # this is called method.DATA now
    Sb = zeros(s);
    sigs = sample(1:2,prob.n,replace=true).*2.-3;
    for i =1:s
           method.DATA[i,:] = sum(sigs[indM[i,:]].*prob.A[indM[i,:],:],1);
           Sb[i] = sum(sigs[indM[i,:]].*prob.b[indM[i,:]])
    end
    SAS = zeros(s,s); 
     for i =1:s
          SAS[i,:] = sum(sigs[indM[i,:]]'.*method.DATA[:,indM[i,:]],2);
     end
     y = SAS\(method.DATA*x-Sb);   # solving (S^TAS) y = (S^TAx-S^Tb)  
     for i =1:s #adding on S^T y
        x[indM[i,:]] = x[indM[i,:]]-sigs[indM[i,:]].*y[i];
     end
end
        
        
        
        
        
        
        
