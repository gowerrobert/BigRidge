

function boot_grad(prob::Prob,options::MyOptions)
    flopsperiter = (prob.n)^2;
    name = "grad";
    stepmethod = step_grad
    method = Method(flopsperiter,name,step_grad,boot_grad)
    return method;
end

function step_grad(prob::Prob, x::Array{Float64}, d::Array{Float64}, options::MyOptions)
    d = - (prob.A*x-prob.b);
end
        
        
        
        
        
        
        
