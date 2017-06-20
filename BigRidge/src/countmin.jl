include("Countstream.jl");

function boot_countmin(prob::Prob,options::MyOptions)
    flopsperiter = (options.sketchsize)^3; # Re-think this
    name = "countmin";
    stepmethod = step_countmin
    method = Method(flopsperiter,name,step_countmin,boot_countmin)
    return method;
end

#function step_countmin(as::Array{Float64,2},bs::Array{Float64,2},sas::Array{Float64,2},S::Array{Float64,2},prob::Prob, x::Array{Float64}, options::MyOptions )

function step_countmin(prob::Prob, x::Array{Float64}, options::MyOptions )

s = options.sketchsize;

         b = reshape(prob.b,1,prob.n);

    mat = [prob.A;b];
    M = Countstream(mat,s);
    
    as = M[1:end-1,:];
    bs = M[end,:]';
    S = Countstream(eye(prob.n),s); # S
    sas = S'*as;

   # as = Countstream(prob.A,s); # A * S
  #  bs = Countstream(prob.b',s); # b * S
  #  sas = Countstream(as',s); # S^T * A * S
    
    
   
     vect = as'*x-bs'; # 
    y = pinv(sas)*vect;   # solving (S^TAS) y = (S^TAx-S^Tb)   

STy = S*y;     # calculating STy 
    x[:] = x[:] -STy; 
end
        
        
        
        
        
        
        