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

   A = broadcast(*, prob.A, sgn); # flip the signs of each column w.p. 50%
         b = reshape(prob.b,1,prob.n);
    
  # b = broadcast(*, b, sgn); # flip the signs of each column w.p. 50%

   # Id = broadcast(*, eye(prob.n), sgn); # flip the signs of each column w.p. 50%


  

    as = Countstream(A,s,ll); # A * S
    
  # S = Countstream(Id,s,ll); # S

        S = prob.A \ as;


    
 #  sas = Countstream(as',s,ll); # S^T * A * S
    sas = S'*as;
  #  bs = Countstream(b,s,ll); # b * S

bs = b*S;

    
     vect = as'*x-bs'; # 
    y = sas\vect;   # solving (S^TAS) y = (S^TAx-S^Tb)   

    STy = S*y;     # calculating STy 
    x[:] = x[:] -STy; 

end
        
        
        
        
        
        
        
