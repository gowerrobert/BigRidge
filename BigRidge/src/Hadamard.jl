include("hada.jl")

function boot_Hadamard(prob::Prob,options::MyOptions)

    name = "Hadamard";
    stepmethod = step_Hadamard
    flopsperiter = (options.sketchsize)^3;


    method = Method(flopsperiter,name,step_Hadamard,boot_Hadamard)
    return method;
end

function step_Hadamard(prob::Prob, x::Array{Float64}, options::MyOptions )
idx = sample(1:prob.n,options.sketchsize,replace=true);
    
    # SA(x_k+1 - x_k) = -S(Ax-b)
        
    
    M = [prob.A prob.b];
    
    mat = hada(M,idx); # that is S * M
    sa = mat[:,1:end-1];  # S * A
    sb = mat[:,end];   # S * b
    
   
    vect = sa*x-sb;
x[:] = x[:] - sa\vect;
end
        
        
        
        
        
        
        
