

function boot_CD(prob::Prob,options::MyOptions)
        #              inverting A[s,s]  +  calculating A[s,:]*x  +      adding x and subtracting prob.b[s]
    flopsperiter = (options.sketchsize)^3  + prob.n*options.sketchsize + 2*prob.n; # Re-think this
    name = string("CD-",s);
    stepmethod = step_CD
    grad = zeros(prob.n);
    grad[:] = prob.b;
    method = Method(flopsperiter,name,step_CD,boot_CD,grad)
    return method;
end

function step_CD(prob::Prob, x::Array{Float64}, options::MyOptions, method::Method )
s = sample(1:prob.n,options.sketchsize,replace=false);
x[s] =  x[s] -prob.A[s,s]\   #x[s] -(prob.A[s,s]\(prob.A[s,:]*x -prob.b[s])) ;
end
        
        
        
        
        
        
        
