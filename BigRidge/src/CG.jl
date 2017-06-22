

function boot_CG(prob::Prob,options::MyOptions)
   # flopsperiter = (prob.n)^2;
    name = "CG";
    stepmethod = step_CG
    method = CGMethod(0.0,name,step_CG,boot_CG,zeros(prob.n),zeros(prob.n),0.0 )

method.flopsperiter = 2*countnz(prob.A)+7*prob.n; #nnz(A) + 7n
method.r[:] = prob.b;
method.p[:] = method.r;# A'*(b-A*x);
method.rr = dot(method.r,method.r); 
    return method;
end

function step_CG(prob::Prob, x::Array{Float64}, options::MyOptions, method::CGMethod)
Ap = prob.A*method.p;
alpha = (method.rr)/(dot(method.p,Ap));
d= alpha*method.p;
method.r[:] = method.r-alpha*Ap;
rrnew = dot(method.r,method.r);
beta = (rrnew)/(method.rr);
method.rr = rrnew;
method.p[:] = method.r+beta*method.p;
x[:] = x[:]+d;
end
        
        
        
        
        
        
        
