include("hada.jl")
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
    idx = sample(1:prob.n,options.sketchsize,replace=false);
        
    #M = [prob.A prob.b];   # <-- too expensive to form a and store this matrix at every iteration
    sa = hada(A,idx); # mat[:,1:end-1];  # S * A
    sb = hada(b,idx) # mat[:,end];   # S * b
    sas = hada(sa',idx);  # SAS^T
    # Implementation of
    # x =   x -S(S^TAS)^(-1)(S^TAx-S^Tb).
    vect = sa*x-sb; # 
    y = sas\vect;   # solving (S^TAS) y = (S^TAx-S^Tb)    
    S = hada(eye(prob.n),idx);  # NEED a more efficient function for calculating S
    STy = S'*y;     # calculating STy 
    x[:] = x[:] -STy; 
    
end
        
        
        
        
        
        
        
