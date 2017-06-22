

function boot_CD(prob::Prob,options::MyOptions)
        #              inverting A[s,s]  +  calculating A[s,:]*x  +      adding x[s] and subtracting prob.b[s]
    flopsperiter = (options.sketchsize)^3  + prob.n*options.sketchsize + 2*options.sketchsize; # Re-think this
    name = string("CD-",options.sketchsize);
    stepmethod = step_CD
    grad = zeros(prob.n);
    grad[:] = -prob.b[:];
    method = Method(flopsperiter,name,step_CD,boot_CD,[],[],grad,[],[])
    return method;
end

function step_CD(prob::Prob, x::Array{Float64}, options::MyOptions, method::SketchMethod )
s = sample(1:prob.n,options.sketchsize,replace=false);    
#d = -prob.A[s,s]\method.DATA[s];   # allocates s size vec
#x[s] =  x[s] +d;
#method.DATA =    method.DATA + prob.A[:,s]*d #update the gradient     
   
# Implementation with "brute force" update of gradient    
method.Sb[s] = prob.A[s,:]*x -prob.b[s];   #update the gradient        
x[s] =  x[s] -prob.A[s,s]\method.Sb[s];   #x[s] -(prob.A[s,s]\(prob.A[s,:]*x -prob.b[s])) ;
end
        
        
        
        
        
        
        
