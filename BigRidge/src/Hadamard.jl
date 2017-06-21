include("hada.jl");
function boot_Hadamard(prob::Prob,options::MyOptions)

name = "Hadamard";
stepmethod = step_Hadamard
# cost of linear solve #cost of calculating SA  #cost of adding x and subtracting Sb.
flopsperiter = (options.sketchsize)^3 + prob.n*options.sketchsize*convert(Int64,ceil(log(options.sketchsize))) +2*prob.n ;
println("flopsperiter: ", flopsperiter)

method = Method(flopsperiter,name,step_Hadamard,boot_Hadamard,[])
return method;
end


function step_Hadamard(prob::Prob, x::Array{Float64}, options::MyOptions, method::Method )

s = options.sketchsize;


idx = sample(1:prob.n,s,replace = false);

sgn = reshape(sample(1:2,prob.n) * 2 - 3,prob.n,1); # one half are +1 and the rest are −1




sa = hada(prob.A,idx,sgn);   # S * A

S = hada(eye(prob.n),idx,sgn);   # S


sb = hada(reshape(prob.b,prob.n,1),idx,sgn);
#sb = S*prob.b; # S * b
sas = hada(sa',idx,sgn);
#sas = S * sa'; # SAS^T

vect = sa*x-sb;
y = sas\vect;   # solving (S^TAS) y = (S^TAx-S^Tb)

#y_n = zeros(prob.n,1);
#y_n[idx] = y;

#x[:] = x[:] -hada(y_n,1:prob.n,sgn);
x[:] = x[:] - S'*y;
    
end









