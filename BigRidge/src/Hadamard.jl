include("hada.jl");
function boot_Hadamard(prob::Prob,options::MyOptions)

name = "Hadamard";
stepmethod = step_Hadamard
# cost of linear solve #cost of calculating SA  #cost of adding x and subtracting Sb.
flopsperiter = (options.sketchsize)^3 + prob.n*options.sketchsize*convert(Int64,ceil(log(options.sketchsize))) +2*prob.n ;
println("flopsperiter: ", flopsperiter)

method = Method(flopsperiter,name,step_Hadamard,boot_Hadamard)
return method;
end


function step_Hadamard(prob::Prob, x::Array{Float64}, options::MyOptions )

s = options.sketchsize;


idx = sample(1:prob.n,s,replace = false);
prob.b = reshape(prob.b,prob.n,1);



sa = hada(prob.A,idx);   # S * A
    
    S = sa/prob.A;
    
sb = hada(prob.b,idx)   # S * b
  #  sb = S*prob.b;

   # sas = hada(sa',idx);  # SAS^T
    sas = S * sa';
    
vect = sa*x-sb;
y = sas\vect;   # solving (S^TAS) y = (S^TAx-S^Tb)

#y_n = zeros(prob.n,1);
#y_n[idx] = y;
#STy =  hada(y_n,1:prob.n);
#x[:] = x[:] -STy;

    x[:] = x[:] -S'*y;



end








