include("hada.jl")

function boot_Hadamard(prob::Prob,options::MyOptions)

    name = "Hadamard";
    stepmethod = step_Hadamard
    #              cost of linear solve     #cost of calculating SA                                 #cost of adding x and subtracting Sb.
    flopsperiter = (options.sketchsize)^3 + prob.n*prob.options.sketchsize*log(options.sketchsize) +2*prob.n ; 


    method = Method(flopsperiter,name,step_Hadamard,boot_Hadamard)
    return method;
end

function step_Hadamard(prob::Prob, x::Array{Float64}, options::MyOptions )
    idx = sample(1:prob.n,options.sketchsize,replace=false);
    # SA(x_k+1 - x_k) = -S(Ax-b)
        
    
    M = [prob.A prob.b];
    
    mat = hada(M,idx); # that is S * M
    sa = mat[:,1:end-1];  # S * A
    sb = mat[:,end];   # S * b
    
   
    vect = sa*x-sb;
    x[:] = x[:] - sa\vect;
    # I suggest implementing this differently:
    # x =   x -S(S^TAS)^(-1)(S^TAx-S^Tb).
    # in other words
    # first solve a small positive definite linear system
    # (S^TAS) y = (S^TAx-S^Tb)    
    # then
    # x = x - S*y.
end
        
        
        
        
        
        
        
