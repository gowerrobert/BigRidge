function boot_method(method_name::AbstractString, prob,options::MyOptions)
@match method_name begin
    "grad"  => method = boot_grad(prob,options);
    "Hadamard"     => method = boot_Hadamard(prob,options);
    "countmin" => method = boot_countmin(prob,options);
    "CG" => method = boot_CG(prob,options);
    "CD" => method = boot_CD(prob,options);   
    "rademacher" => method = boot_rademacher(prob,options);  
                  _ => method = "METHOD DOES NOT EXIST"
    end
    return method;
end

