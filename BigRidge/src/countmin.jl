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
ll = sample(1:s, prob.n); # sample n items from [s] with replacement
sgn = reshape(sample(1:2,prob.n) * 2 - 3, 1,prob.n); # one half are +1 and the rest are −1
    
prob.b = reshape(prob.b,1,prob.n);

as = Countstream(prob.A,s,ll,sgn); # A * S
as = Countstream(prob.b,s,ll,sgn); # b * S

    
S = prob.A \ as; # It has to be something else;
vect = as'*x-(bs)'; 
x[:] = x[:] -S*(S'*as\vect); #solving (S^TAS) y = (S^TAx-S^Tb)


end







