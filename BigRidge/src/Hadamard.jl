
function boot_Hadamard(prob::Prob,options::MyOptions)
    flopsperiter = (options.sketchsize)^3; # Re-think this
    name = "Hadamard";
    stepmethod = step_Hadamard
    method = Method(flopsperiter,name,step_Hadamard,boot_Hadamard)
    return method;
end

function step_Hadamard(prob::Prob, x::Array{Float64}, d::Array{Float64}, options::MyOptions )
idx = sample(1:prob.n,options.sketchsize,replace=true);
    
   S = zeros(prob.n,options.sketchsize);
    
x[idx] = x[idx] - S*(S'*prob.A*S)\S'*(prob.A*x -prob.b) ;
end
        
        
        
        
        
        
        
