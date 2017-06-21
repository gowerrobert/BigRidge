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



b = reshape(prob.b,1,prob.n);




as = Countstream(prob.A,s,ll); # A * S


S = prob.A \ as;



sas = S'*as;

bs = b*S;


vect = as'*x-bs'; #
y = sas\vect;   # solving (S^TAS) y = (S^TAx-S^Tb)

x[:] = x[:] -S*y;

end







