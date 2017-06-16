

function boot_grad(prob::Prob,options::MyOptions)
    flopsperiter = (prob.n)^2;
    name = "grad";
    stepmethod = step_grad
    method = Method(flopsperiter,name,step_grad,boot_grad)
    return method;
end

function step_grad(prob::Prob, x::Array{Float64}, options::MyOptions)
     # Fix this later!
    #Make this efficient by storing the previous gradient and figuring out how to use it to calculate next grad
    # e.g. A(x_2-x_1) = \nabla f(x_1) - \nabla f(x_2)
    g= (prob.A*x-prob.b);
    alpha = norm(g)^2/dot(g,prob.A*g);
    d = - alpha*(prob.A*x-prob.b);
    x[:] = x[:]+d;
           # println("norm d:  ",norm(d))
end
        
        
        
        
        
        
        
