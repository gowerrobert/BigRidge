
include("hada.jl");
function boot_Hadamard(prob::Prob,options::MyOptions)

    s = options.sketchsize;

name = "Hadamard";
stepmethod = step_Hadamard
# cost of linear solve #cost of calculating SA  #cost of adding x and subtracting Sb.
flopsperiter = (options.sketchsize)^3 + prob.n*options.sketchsize*convert(Int64,ceil(log(options.sketchsize))) +2*prob.n ;
println("flopsperiter: ", flopsperiter)
    
    SA = zeros(s,prob.n); # saving space for sketched matrix SA
    SAS = zeros(s,s);
    Sb = zeros(s);
    ind = sample(1:prob.n,s,replace=false); 
    sigs = zeros(s);

    method = SketchMethod(flopsperiter,name,step_Hadamard,boot_Hadamard,SA,SAS,Sb,ind,sigs)
return method;
end




function step_Hadamard(prob::Prob, x::Array{Float64}, options::MyOptions, method::SketchMethod )

s = options.sketchsize;


 sample!(1:prob.n,method.ind, replace = false);

    
method.SA[:] = hada(prob.A,method.ind);   # S * A

    method.Sb[:] = hada(reshape(prob.b,prob.n,1),method.ind);

    method.SAS[:] = hada(method.SA',method.ind);


     y = method.SAS\(method.SA*x-method.Sb);   # solving (S^TAS) y = (S^TAx-S^Tb)  


y_n = zeros(prob.n,1);
y_n[method.ind] = y;

x[:] = x[:] -hada(y_n,1:prob.n);
    
end











