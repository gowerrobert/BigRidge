function boot_rademacher(prob::Prob,options::MyOptions)
    s  =options.sketchsize;                   
    name = string("rademacher-",s,"-",Int(options.AUX[2]));
    if(Int(options.AUX[2])*options.sketchsize > prob.n) # check density size
          println("change rho density to n/s: ", Int(floor(prob.n/options.sketchsize)) )
          options.AUX[2] = Int(floor(prob.n/options.sketchsize));
    end
    rho = Int(options.AUX[2]);    # density
    if(Int(options.AUX[1])==1) # with sign flipping
        stepmethod =  step_rademacher_sign_flip;
        sigs = sample(1:2,rho*s,replace=true).*2.-3;
        name = string(name,"-sgn");
    else
        stepmethod = step_rademacher;
        sigs =[];
    end  
    SA = zeros(s,prob.n); # saving space for sketched matrix SA
    SAS = zeros(s,s);
    Sb = zeros(s);
    ind = zeros(rho*s);
    #ind = sample(1:prob.n,rho*s,replace=false); 
    flopsperiter = s^3+  #  solving  SAS\(SA*x-Sb)  
    (s+1)*(prob.n+s)*Int(ceil(log(rho)))+ #  calculating sketch SA, Sb and SAS 
    prob.n*(s+1); #computing SA*x-Sb 
    method = SketchMethod(flopsperiter,name,stepmethod,boot_rademacher,SA,SAS,Sb,ind,sigs)
    return method;
end

    
function step_rademacher_sign_flip(prob::Prob, x::Array{Float64}, options::MyOptions, method::SketchMethod )
    rho = Int(options.AUX[2]); # the density per row
    s = options.sketchsize;
    sample!(1:prob.n, method.ind ; replace=false)  # get a sample of rows 
    method.ind =reshape(method.ind, s,rho);
    sample!(1:2,method.sigs; replace=true);
    method.sigs[:] =method.sigs.*2.-3;
    method.sigs = reshape(method.sigs, s,rho);
    for i =1:s
           method.SA[i,:] = sum(method.sigs[i,:].*prob.A[method.ind[i,:],:],1);
           method.Sb[i] = sum(method.sigs[i,:].*prob.b[method.ind[i,:]])
    end
    for i =1:s
          method.SAS[i,:] = sum(method.sigs[i,:]'.*method.SA[:,method.ind[i,:]],2);
    end
    method.Sb[:] =method.SA*x-method.Sb;
    method.Sb[:] = method.SAS\method.Sb;   # solving (S^TAS) y = (S^TAx-S^Tb)  
    for i =1:s #adding on S^T y
       x[method.ind[i,:]] -= method.sigs[i,:].*method.Sb[i];
    end
end
 
function step_rademacher(prob::Prob, x::Array{Float64}, options::MyOptions, method::SketchMethod )
# With no sign swapping
    rho = Int(options.AUX[2]); # the density per row
    s = options.sketchsize;
    sample!(1:prob.n, method.ind ; replace=false)  # get a sample of rows 
    method.ind =reshape(method.ind, s,rho);
    for i =1:s
           method.SA[i,:] = sum(prob.A[method.ind[i,:],:],1);
           method.Sb[i] = sum(prob.b[method.ind[i,:]])
    end
     for i =1:s
          method.SAS[i,:] = sum(method.SA[:,method.ind[i,:]],2);
     end
     method.Sb[:] = method.SAS\(method.SA*x-method.Sb);   # solving (S^TAS) y = (S^TAx-S^Tb)  
     for i =1:s #adding on S^T y
        x[method.ind[i,:]] -= method.Sb[i];
        #x[method.ind[i,:]] = x[method.ind[i,:]]-method.Sb[i];
     end
end
        

function step_rademacher_full_indexed(prob::Prob, x::Array{Float64}, options::MyOptions, method::SketchMethod )
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
        
        
        
