

function boot_countmin(prob::Prob,options::MyOptions)
    flopsperiter = (options.sketchsize)^3; # Re-think this
    name = "countmin";
    stepmethod = step_countmin
    method = Method(flopsperiter,name,step_countmin,boot_countmin)
    return method;
end

function step_countmin(prob::Prob, x::Array{Float64}, d::Array{Float64}, options::MyOptions )
s = sample(1:prob.n,options.sketchsize,replace=false);
x[s] = x[s] -(prob.A[s,s]\(prob.A[s,:]*x -prob.b[s])) ;
end
        
        
        
        
        
        
        
