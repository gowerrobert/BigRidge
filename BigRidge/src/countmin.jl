#include("Countstream.jl");

function boot_countmin(prob::Prob,options::MyOptions)
flopsperiter = (options.sketchsize)^3; # Re-think this
name = "countmin";
stepmethod = step_countmin
method = Method(flopsperiter,name,step_countmin,boot_countmin,[])
return method;
end


function step_countmin(prob::Prob, x::Array{Float64}, options::MyOptions, method::Method )

s = options.sketchsize;
ll = sample(1:s, prob.n); # sample n items from [s] with replacement
sgn = reshape(sample(1:2,prob.n) * 2 - 3,prob.n,1); # one half are +1 and the rest are −1
    
    
    
    m,n = size(prob.A);


#A_sign = broadcast(*, prob.A, sgn); # flip the signs of each column w.p. 50%
#b_sign = broadcast(*, prob.b, sgn);
#Id_sign = broadcast(*, eye(prob.n), sgn);

as = zeros(m, s); # initialize A*S
bs = zeros(1, s);
S = zeros(n, s);
for j=1:n
        sgn = sample(1:2,1)*2-3;
    
    as[:, ll[j]] = as[:, ll[j]] + sgn[1]*prob.A[:, j];
        bs[ll[j]] +=  sgn[1]*prob.b[j];

    S[j, ll[j]] += sgn[1];


end
    
    
    
    
#as = Countstream(prob.A,s,ll,sgn); # A * S
#bs = Countstream(reshape(prob.b,1,prob.n),s,ll,sgn); # b * S
#S = Countstream(eye(prob.n),s,ll,sgn); # A * S


    
#S = prob.A \ as; # It has to be something else;
    
vect = as'*x-bs'; 
x[:] = x[:] -S*(S'*as\vect); #solving (S^TAS) y = (S^TAx-S^Tb)


end







