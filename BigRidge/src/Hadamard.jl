include("hada.jl");
function boot_Hadamard(prob::Prob,options::MyOptions)

    name = "Hadamard";
    stepmethod = step_Hadamard
    #              cost of linear solve     #cost of calculating SA                                 #cost of adding x and subtracting Sb.
    flopsperiter = (options.sketchsize)^3 + prob.n*options.sketchsize*convert(Int64,ceil(log(options.sketchsize))) +2*prob.n ; 
    println("flopsperiter: ", flopsperiter)

    method = Method(flopsperiter,name,step_Hadamard,boot_Hadamard)
    return method;
end


  function step_Hadamard(prob::Prob, x::Array{Float64}, options::MyOptions )

sgn = reshape(sample(1:2,prob.n) * 2 - 3, 1,prob.n);
    prob.A = broadcast(*, prob.A, sgn);
    
    s = options.sketchsize;

    
 idx = sample(1:prob.n,s,replace = false);
    prob.b = reshape(prob.b,prob.n,1);
        
    sa = hada(prob.A,idx); # mat[:,1:end-1];  # S * A
    sb = hada(prob.b,idx) # mat[:,end];   # S * b
    sas = hada(sa',idx);  # SAS^T    
    S = hada(eye(prob.n),idx);  # NEED a more efficient function for calculating S
  
    vect = sa*x-sb; # 
    y = sas\vect;   # solving (S^TAS) y = (S^TAx-S^Tb)    
   
     STy = S'*y;     # calculating STy 
    x[:] = x[:] -STy; 
    
    
    
# y_n = zeros(prob.n,1);
 #   y_n[idx] = y;
  #STy =  hada(y_n,1:prob.n);
   #     x[:] = x[:] -STy; 
    
    
end
        
        
        
        
        
        
        

