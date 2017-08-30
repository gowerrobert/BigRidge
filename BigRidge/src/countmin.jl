
function boot_countmin(prob::Prob,options::MyOptions)
flopsperiter = (options.sketchsize)^3; # Re-think this, I believe SA costs nnz(A), so the total number of flops should be around  nnz(A) + s^3.
name = "countmin";
stepmethod = step_countmin;
    s = options.sketchsize;

     SA = zeros(prob.n,s); # saving space for sketched matrix AS (that is AS here and not SA)
    SAS = zeros(prob.n+s,s); # here we use this yield to store S instead of SAS
    # S = SAS[1:prob.n,:]   and 'SAS' = SAS[prob.n+1:end,:]
    Sb = zeros(1,s);
    ind = sample(1:s, prob.n); # sample n items from [s] with replacement


    sigs = sample(1:2,1);


    
    
method = SketchMethod(flopsperiter,name,step_countmin,boot_countmin,SA,SAS,Sb,ind,sigs);
return method;
end


function step_countmin(prob::Prob, x::Array{Float64}, options::MyOptions, method::SketchMethod )

s = options.sketchsize;
sample!(1:s, method.ind); # sample n items from [s] with replacement
#sgn = reshape(sample(1:2,prob.n) * 2 - 3,prob.n,1); # one half are +1 and the rest are −1
    

for j=1:prob.n
        
    sample!(1:2,method.sigs);
    method.sigs[:] = 2.*method.sigs .- 3
    method.SA[:, method.ind[j]] +=  method.sigs[1].*prob.A[:, j];
    method.Sb[method.ind[j]] +=  method.sigs[1].*prob.b[j];

    method.SAS[j, method.ind[j]] += method.sigs[1]; #that's S in fact that we store in the yield SAS

end

    for i =1:s
        for j = 1:s
        method.SAS[prob.n+i,j] = sum(method.SAS[1:prob.n,i].*method.SA[:,j]);
        end
    end
    
    x[:] -= method.SAS[1:prob.n,:]*(method.SAS[prob.n+1:end,:]\(method.SA'*x-method.Sb')); #solving (S^TAS) y = (S^TAx-S^Tb) and multiplying y by S
    
  #  x[:] -= method.SAS*((method.SAS'*method.SA)\(method.SA'*x-method.Sb')); 
end







