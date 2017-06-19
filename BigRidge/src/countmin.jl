include("Countstream.jl");

function boot_countmin(prob::Prob,options::MyOptions)
    flopsperiter = (options.sketchsize)^3; # Re-think this
    name = "countmin";
    stepmethod = step_countmin
    method = Method(flopsperiter,name,step_countmin,boot_countmin)
    return method;
end

function step_countmin(prob::Prob, x::Array{Float64}, options::MyOptions )
s = options.sketchsize;
    
    C = Countstream(prob.A,s); # S^T
    sa = C*prob.A;
    sb = C*prob.b;
    
  vect = sa*x-sb;
    x[:] = x[:] - sa\vect;
end
        
        
        
        
        
        
        
